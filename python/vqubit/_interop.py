"""Quirk import and export.

Quirk (https://algassert.com/quirk) draws circuits column by column and stores
them as a JSON payload inside the URL. This module writes that payload (so a
vqubit circuit can be opened there) and reads it back.

The gate set is Quirk's: ``H``, the powers ``X^t``/``Y^t``/``Z^t`` (which is how
rotations are written), the controlled versions of those, ``Swap`` and
``Measure``.
"""

from __future__ import annotations

import base64
import json
import math
from typing import Any

from ._circuits import Circuit
from ._gates import (
	CNOTGate,
	ControlledGate,
	CZGate,
	HPowGate,
	IdentityGate,
	MeasurementGate,
	RxGate,
	RyGate,
	RzGate,
	SWAPGate,
	XPowGate,
	YPowGate,
	ZPowGate,
	is_plain_gate,
)
from ._ops import Operation
from ._qubits import NamedQubit, Qid

_QUIRK_URL = "https://algassert.com/quirk#circuit="

_AXES = {"X": XPowGate, "Y": YPowGate, "Z": ZPowGate}
_FRACTIONS = {"½": 0.5, "¼": 0.25, "⅛": 0.125}


def _format_exponent(exponent: float) -> str:
	"""Quirk's exponent text: ``0.25``, or ``-1`` for a whole number."""
	number = float(exponent)
	if number == int(number):
		return str(int(number))
	return repr(round(number, 12))


def _parse_exponent(text: str) -> float:
	"""Read the exponent Quirk writes, including its ``½``/``¼`` forms."""
	if text in _FRACTIONS:
		return _FRACTIONS[text]
	sign = -1.0 if text.startswith("-") else 1.0
	stripped = text[1:] if sign < 0 else text
	if stripped in _FRACTIONS:
		return sign * _FRACTIONS[stripped]
	return sign * float(stripped)


def _single_qubit_column(gate: Any) -> str | None:
	"""The Quirk column entry for a single-qubit gate, or ``None``."""
	if isinstance(gate, MeasurementGate):
		return "Measure"
	if isinstance(gate, HPowGate) and gate.exponent == 1:
		return "H"
	if isinstance(gate, IdentityGate):
		return None
	for axis, pow_gate in _AXES.items():
		if isinstance(gate, pow_gate):
			if gate.exponent == 1:
				return axis
			return f"{axis}^{_format_exponent(gate.exponent)}"
	# A rotation about an axis is the corresponding power, up to a global phase.
	rotations = ((RxGate, "X"), (RyGate, "Y"), (RzGate, "Z"))
	for rotation, axis in rotations:
		if isinstance(gate, rotation):
			return f"{axis}^{_format_exponent(gate.radians / math.pi)}"
	return None


def _two_qubit_column(gate: Any) -> tuple[str, str] | None:
	if isinstance(gate, CNOTGate) and is_plain_gate(gate):
		return ("•", "X")
	if isinstance(gate, CZGate) and is_plain_gate(gate):
		return ("•", "•")
	if isinstance(gate, SWAPGate) and is_plain_gate(gate):
		return ("Swap", "Swap")
	if isinstance(gate, ControlledGate) and gate.num_controls == 1 and all(gate.control_values):
		# A controlled power gate is a ``•`` over the power's own symbol.
		symbol = _single_qubit_column(gate.sub_gate)
		if symbol is not None and symbol != "Measure":
			return ("•", symbol)
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
			if isinstance(gate, MeasurementGate):
				# A measurement can act on several qubits at once, so it is written
				# before the one- and two-qubit dispatch.
				for qubit in operation.qubits:
					column[index[qubit]] = "Measure"
				continue
			if len(operation.qubits) == 1:
				entry = _single_qubit_column(gate)
				if entry is None:
					if isinstance(gate, IdentityGate):
						continue
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


def _gate_from_entry(entry: str) -> Any:
	"""The gate a Quirk column entry names, or ``None`` for a wire."""
	if entry in (1, "-"):
		return None
	if entry == "H":
		return HPowGate(exponent=1)
	axis, _, exponent = entry.partition("^")
	if axis in _AXES:
		return _AXES[axis](exponent=_parse_exponent(exponent) if _ else 1)
	return None


def quirk_json_to_circuit(json_data: Any) -> Circuit:
	"""Build a circuit from a Quirk circuit JSON document."""
	data = json.loads(json_data) if isinstance(json_data, str) else json_data
	cols = data["cols"]
	names = data.get("qubits") or [f"q{i}" for i in range(len(cols[0]))]
	qubits = [NamedQubit(name) for name in names]

	operations: list[Operation] = []
	for column in list(cols):
		measured: list[Qid] = []
		for position, entry in enumerate(column):
			if entry == "Measure":
				measured.append(qubits[position])
				column[position] = 1
		if measured:
			operations.append(_measure(measured))
		position = 0
		while position < len(column):
			entry = column[position]
			if entry in (1, "-", "Measure"):
				position += 1
				continue
			qubit = qubits[position]
			if entry == "•":
				control = qubit
				target = _find_controlled_target(column, position)
				if target is None:
					raise ValueError("Quirk control has no target")
				if column[target] == "•":
					# ``•`` over ``•`` is what Quirk writes for CZ.
					operations.append(CZGate().on(control, qubits[target]))
				else:
					gate = _gate_from_entry(column[target])
					if gate is None:
						raise ValueError(f"unsupported Quirk control target {column[target]!r}")
					if isinstance(gate, XPowGate) and gate.exponent == 1:
						operations.append(CNOTGate().on(control, qubits[target]))
					elif isinstance(gate, ZPowGate):
						operations.append(CZGate(exponent=gate.exponent).on(control, qubits[target]))
					else:
						operations.append(ControlledGate(gate, num_controls=1).on(control, qubits[target]))
				column[target] = 1
				column[position] = 1
				position += 1
				continue
			if entry == "Swap":
				for other in range(len(column)):
					if other != position and column[other] == "Swap":
						operations.append(SWAPGate().on(qubit, qubits[other]))
						column[other] = 1
						column[position] = 1
						break
				position += 1
				continue
			gate = _gate_from_entry(entry)
			if gate is None:
				raise ValueError(f"unsupported Quirk gate {entry!r}")
			operations.append(gate.on(qubit))
			column[position] = 1
			position += 1
	return Circuit(operations)


def _find_controlled_target(column: list[Any], control_position: int) -> int | None:
	"""The position of the gate a ``•`` controls, if there is one."""
	for position, entry in enumerate(column):
		if position == control_position:
			continue
		# A bare axis gate (``X``, ``Y``, ``Z``), a power of one, ``H`` and the
		# ``•`` that stands for the target of a CZ are all valid targets.
		if entry == "•" or entry == "H":
			return position
		if isinstance(entry, str) and entry.partition("^")[0] in _AXES:
			return position
	return None


def _measure(qubits: list[Qid]) -> Operation:
	from ._gates import measure

	return measure(*qubits, key="m")


def quirk_url_to_circuit(url: str) -> Circuit:
	"""Build a circuit from a Quirk URL."""
	payload = url.split("circuit=", 1)[1]
	payload += "=" * (-len(payload) % 4)
	decoded = base64.urlsafe_b64decode(payload).decode()
	return quirk_json_to_circuit(decoded)


def _import_cirq() -> Any:
	"""Cirq, imported on demand: it is optional and not a dependency.

	The import is dynamic so that a type checker run without Cirq installed does
	not report it as a missing module (the same trick ``vqubit._compat`` uses for
	numpy and sympy).
	"""
	import importlib

	try:
		return importlib.import_module("cirq")
	except ImportError as exc:  # pragma: no cover - only without cirq installed
		raise ImportError("the cirq bridge needs cirq installed") from exc


def to_cirq(circuit: Circuit) -> Any:
	"""Convert a vqubit circuit to a Cirq one (Cirq must be installed).

	The conversion goes through the JSON format the two libraries share, so it
	covers the gates that format carries (the same ones Cirq serializes).
	"""
	from ._json import to_json

	cirq = _import_cirq()
	return cirq.read_json(json_text=to_json(circuit))


def from_cirq(circuit: Any) -> Circuit:
	"""Convert a Cirq circuit to a vqubit one (Cirq must be installed)."""
	from ._json import read_json

	cirq = _import_cirq()
	return read_json(cirq.to_json(circuit))


__all__ = ["from_cirq", "quirk_json_to_circuit", "quirk_url", "quirk_url_to_circuit", "to_cirq"]
