"""Devices.

A device describes where a circuit may run: which qubits exist and which
connections between them are allowed. vqubit ships local versions of the vendor
devices -- Google's Sycamore, an ion trap, a neutral-atom register -- so code
written for those APIs runs against the local simulators with no network and no
dependencies.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from typing import Any

from ._gates import Gate
from ._linalg import Matrix
from ._ops import Operation
from ._qubits import GridQubit, LineQubit, Qid
from ._value import value_equality


class DeviceMetadata:
	"""What a device says about itself."""

	def __init__(
		self,
		*,
		name: str = "device",
		qubit_set: Iterable[Qid] | None = None,
		gate_durations: dict[Any, Any] | None = None,
	) -> None:
		self.name = name
		self.qubit_set = frozenset(qubit_set) if qubit_set is not None else None
		self.gate_durations = dict(gate_durations or {})

	def __repr__(self) -> str:
		return f"vqubit.DeviceMetadata(name={self.name!r})"

	def duration_of(self, operation: Operation) -> Any:
		"""Return the duration recorded for ``operation``'s gate, or ``None``."""
		return self.gate_durations.get(type(operation.gate).__name__)


class Device:
	"""A set of qubits a circuit may act on."""

	def __init__(self, qubits: Iterable[Qid], *, metadata: DeviceMetadata | None = None) -> None:
		self._qubits = frozenset(qubits)
		if metadata is None:
			metadata = DeviceMetadata(qubit_set=self._qubits)
		self.metadata = metadata

	@property
	def qubits(self) -> frozenset[Qid]:
		return self._qubits

	@property
	def qubit_set(self) -> frozenset[Qid]:
		return self._qubits

	def neighbors(self, qubit: Qid) -> list[Qid]:
		"""Return the qubits reachable from ``qubit`` in one operation."""
		raise NotImplementedError

	def is_adjacent(self, first: Qid, second: Qid) -> bool:
		"""Report whether an operation on ``first`` and ``second`` is allowed."""
		return second in self.neighbors(first)

	def validate_operation(self, operation: Operation) -> None:
		"""Raise ``ValueError`` if ``operation`` cannot run on this device."""
		for qubit in operation.qubits:
			if qubit not in self._qubits:
				raise ValueError(f"qubit {qubit!r} is not on the device")
		if len(operation.qubits) > 1:
			for i in range(len(operation.qubits)):
				for j in range(i + 1, len(operation.qubits)):
					if not self.is_adjacent(operation.qubits[i], operation.qubits[j]):
						raise ValueError(f"{operation.qubits[i]!r} and {operation.qubits[j]!r} are not adjacent")

	def validate_circuit(self, circuit: Any) -> None:
		"""Raise ``ValueError`` if any operation of ``circuit`` is not allowed."""
		for operation in circuit.all_operations():
			self.validate_operation(operation)

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}({len(self._qubits)} qubits)"


class GridDevice(Device):
	"""A device whose connectivity follows a grid."""

	def __init__(self, qubits: Iterable[GridQubit], *, metadata: DeviceMetadata | None = None) -> None:
		super().__init__(qubits, metadata=metadata)

	@classmethod
	def rect(cls, rows: int, cols: int, *, top: int = 0, left: int = 0) -> GridDevice:
		return cls(GridQubit.rect(range(top, top + rows), range(left, left + cols)))

	def neighbors(self, qubit: Qid) -> list[Qid]:
		if not isinstance(qubit, GridQubit):
			return []
		return [q for q in qubit.neighbors() if q in self._qubits]


class LineDevice(Device):
	"""A device whose qubits sit on a line with nearest-neighbour coupling."""

	def __init__(self, num_qubits: int) -> None:
		super().__init__(LineQubit.range(num_qubits))

	def neighbors(self, qubit: Qid) -> list[Qid]:
		if not isinstance(qubit, LineQubit):
			return []
		result: list[Qid] = []
		for offset in (-1, 1):
			position = qubit.x + offset
			if position < 0:
				continue
			candidate = LineQubit(position)
			if candidate in self._qubits:
				result.append(candidate)
		return result


@value_equality
class SycamoreGate(Gate):
	"""Google's fSim gate, the entangler of the Sycamore processor.

	It is ``fSim(pi/2, pi/6)``: a full iSWAP between the two qubits plus a
	``pi/6`` conditional phase.
	"""

	def _num_qubits_(self) -> int:
		return 2

	def _unitary_(self) -> Matrix:
		conditional_phase = complex(math.cos(math.pi / 6), -math.sin(math.pi / 6))
		return [
			[1 + 0j, 0j, 0j, 0j],
			[0j, 0j, -1j, 0j],
			[0j, -1j, 0j, 0j],
			[0j, 0j, 0j, conditional_phase],
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Syc", "Syc")

	def _value_equality_values_(self) -> Any:
		return ()

	def __str__(self) -> str:
		return "SycamoreGate"


def _sycamore_qubits() -> frozenset[GridQubit]:
	# A representative 54-qubit grid standing in for the real layout, which is
	# specific to the processor. Local emulation only.
	return frozenset(GridQubit.rect(range(6), range(9)))


class Sycamore(GridDevice):
	"""A local emulation of Google's Sycamore device."""

	def __init__(self) -> None:
		qubits = _sycamore_qubits()
		super().__init__(qubits, metadata=DeviceMetadata(name="Sycamore", qubit_set=qubits))


class XmonDevice(GridDevice):
	"""A superconducting-grid device, with optional gate durations."""

	def __init__(self, qubits: Iterable[GridQubit], *, gate_durations: dict[str, Any] | None = None) -> None:
		qubit_set = frozenset(qubits)
		super().__init__(
			qubit_set,
			metadata=DeviceMetadata(name="XmonDevice", qubit_set=qubit_set, gate_durations=gate_durations),
		)


class IonDevice(Device):
	"""An ion-trap device: all qubits are mutually connected."""

	def __init__(self, num_qubits: int) -> None:
		super().__init__(LineQubit.range(num_qubits))
		self._all_to_all = True

	def neighbors(self, qubit: Qid) -> list[Qid]:
		return [q for q in self._qubits if q != qubit]


class NeutralAtomDevice(Device):
	"""A neutral-atom register: qubits are coupled along a line, all-to-all too."""

	def __init__(self, qubits: Iterable[Qid] | None = None, *, num_qubits: int = 10) -> None:
		if qubits is None:
			qubits = LineQubit.range(num_qubits)
		super().__init__(qubits, metadata=DeviceMetadata(name="NeutralAtom"))

	def neighbors(self, qubit: Qid) -> list[Qid]:
		return [q for q in self._qubits if q != qubit]


class AQTDevice(LineDevice):
	"""A local emulation of an AQT ion-trap device (a line of qubits)."""


__all__ = [
	"AQTDevice",
	"Device",
	"DeviceMetadata",
	"GridDevice",
	"IonDevice",
	"LineDevice",
	"NeutralAtomDevice",
	"Sycamore",
	"SycamoreGate",
	"XmonDevice",
]
