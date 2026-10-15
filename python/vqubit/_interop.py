"""Interoperability helpers.

Currently this exposes :func:`quirk_url`, which renders a circuit as a link to
the `Quirk <https://algassert.com/quirk>`_ circuit simulator. Only the gates
that Quirk understands natively are supported.
"""

from __future__ import annotations

import base64
import json
from typing import Any

from ._circuits import Circuit
from ._gates import (
	CNOTGate,
	CZGate,
	HPowGate,
	SWAPGate,
	XPowGate,
	YPowGate,
	ZPowGate,
)

_QUIRK_URL = "https://algassert.com/quirk#circuit="


def _single_qubit_column(gate: Any) -> str | None:
	if isinstance(gate, HPowGate) and gate.exponent == 1:
		return "H"
	if isinstance(gate, XPowGate) and gate.exponent == 1:
		return "X"
	if isinstance(gate, YPowGate) and gate.exponent == 1:
		return "Y"
	if isinstance(gate, ZPowGate):
		if gate.exponent == 1:
			return "Z"
		if gate.exponent == 0.5:
			return "Z^½"
		if gate.exponent == 0.25:
			return "Z^¼"
	return None


def _two_qubit_column(gate: Any) -> tuple[str, str] | None:
	if isinstance(gate, CNOTGate):
		return ("•", "X")
	if isinstance(gate, CZGate):
		return ("•", "•")
	if isinstance(gate, SWAPGate):
		return ("Swap", "Swap")
	return None


def quirk_url(circuit: Circuit) -> str:
	"""Return a Quirk link that reproduces ``circuit``.

	Raises ``TypeError`` if the circuit uses a gate Quirk cannot represent.
	"""
	qubits = sorted(circuit.all_qubits())
	index = {q: i for i, q in enumerate(qubits)}
	size = len(qubits)

	cols: list[list[Any]] = []
	for moment in circuit.moments:
		column: list[Any] = [1] * size
		for operation in moment.operations:
			gate = operation.gate
			if len(operation.qubits) == 1:
				entry = _single_qubit_column(gate)
				if entry is None:
					raise TypeError(f"{gate!r} cannot be written as a Quirk URL")
				column[index[operation.qubits[0]]] = entry
			elif len(operation.qubits) == 2:
				pair = _two_qubit_column(gate)
				if pair is None:
					raise TypeError(f"{gate!r} cannot be written as a Quirk URL")
				column[index[operation.qubits[0]]] = pair[0]
				column[index[operation.qubits[1]]] = pair[1]
			else:
				raise TypeError(f"{gate!r} cannot be written as a Quirk URL")
		cols.append(column)

	payload = {"cols": cols, "qubits": [f"q{i}" for i in range(size)]}
	raw = json.dumps(payload, separators=(",", ":")).encode()
	encoded = base64.urlsafe_b64encode(raw).decode().rstrip("=")
	return _QUIRK_URL + encoded


__all__ = ["quirk_url"]


_QUIRK_SINGLE = {
	"H": "H",
	"X": "X",
	"Y": "Y",
	"Z": "Z",
	"S": "S",
	"T": "T",
	"Z^½": "S",
	"Z^¼": "T",
}


def quirk_json_to_circuit(json_data: Any) -> Circuit:
	"""Build a circuit from a Quirk circuit JSON document."""
	data = json.loads(json_data) if isinstance(json_data, str) else json_data
	cols = data["cols"]
	names = data.get("qubits") or [f"q{i}" for i in range(len(cols[0]))]
	from ._gates import CNOT, CZ, SWAP, H, S, T, X, Y, Z
	from ._qubits import NamedQubit

	qubits = [NamedQubit(name) for name in names]
	gate_map = {"H": H, "X": X, "Y": Y, "Z": Z, "S": S, "T": T}
	operations: list[Any] = []
	for column in cols:
		for index, entry in enumerate(column):
			if entry in (1, "-"):
				continue
			if isinstance(entry, str) and entry in gate_map:
				operations.append(gate_map[entry](qubits[index]))
				continue
			if entry == "•":
				# A control: find the target in the same column.
				target = None
				for other_index, other in enumerate(column):
					if other_index != index and other in ("X", "•", "Swap"):
						target = other_index
						break
				if target is None:
					raise ValueError("Quirk control has no target")
				target_entry = column[target]
				if target_entry == "X":
					operations.append(CNOT(qubits[index], qubits[target]))
				elif target_entry == "•":
					operations.append(CZ(qubits[index], qubits[target]))
				else:
					raise ValueError(f"unsupported Quirk control target {target_entry!r}")
				column[target] = 1
				continue
			if entry == "Swap":
				# Find the other Swap in the column.
				for other_index, other in enumerate(column):
					if other_index != index and other == "Swap":
						operations.append(SWAP(qubits[index], qubits[other_index]))
						column[other_index] = 1
						break
				continue
			raise ValueError(f"unsupported Quirk gate {entry!r}")
	return Circuit(operations)


def quirk_url_to_circuit(url: str) -> Circuit:
	"""Build a circuit from a Quirk URL."""
	payload = url.split("circuit=", 1)[1]
	payload += "=" * (-len(payload) % 4)
	decoded = base64.urlsafe_b64decode(payload).decode()
	return quirk_json_to_circuit(decoded)


__all__ = ["quirk_json_to_circuit", "quirk_url", "quirk_url_to_circuit"]
