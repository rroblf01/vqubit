"""Operations and moments.

An :class:`Operation` is a gate applied to specific qubits. A :class:`Moment`
is a group of operations that act on disjoint qubits, so they can happen at the
same time. A circuit is a sequence of moments.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Any

from ._qubits import Qid
from ._value import value_equality

if TYPE_CHECKING:
	from ._gates import Gate


class Operation:
	"""A quantum operation: something that can be placed in a circuit."""

	@property
	def qubits(self) -> tuple[Qid, ...]:
		raise NotImplementedError

	@property
	def gate(self) -> Gate:
		raise NotImplementedError

	def _circuit_diagram_info_(self) -> tuple[str, ...]:
		return (str(self),) * len(self.qubits)


@value_equality
class GateOperation(Operation):
	"""The application of a :class:`~vqubit.Gate` to specific qubits."""

	def __init__(self, gate: Gate, qubits: Iterable[Qid]) -> None:
		qubits = tuple(qubits)
		if len(qubits) != gate.num_qubits:
			raise ValueError(f"{gate!r} acts on {gate.num_qubits} qubit(s), got {len(qubits)}")
		if any(not isinstance(q, Qid) for q in qubits):
			raise TypeError("GateOperation qubits must all be Qid instances")
		self._gate = gate
		self._qubits = qubits

	@property
	def gate(self) -> Gate:
		return self._gate

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits

	def _value_equality_values_(self) -> Any:
		return (self._gate, self._qubits)

	def _circuit_diagram_info_(self) -> tuple[str, ...]:
		return self._gate._circuit_diagram_info_(self._qubits)

	def __str__(self) -> str:
		return f"{self._gate}{self._qubits if len(self._qubits) != 1 else ''}"

	def __repr__(self) -> str:
		return f"vqubit.GateOperation(gate={self._gate!r}, qubits={self._qubits!r})"


class Moment:
	"""A set of operations that act on disjoint qubits at the same time."""

	def __init__(self, operations: Iterable[Operation] = ()) -> None:
		operations = tuple(operations)
		seen: set[Qid] = set()
		for op in operations:
			for q in op.qubits:
				if q in seen:
					raise ValueError(f"overlapping operations on qubit {q!r}")
				seen.add(q)
		self._operations = operations

	@property
	def operations(self) -> tuple[Operation, ...]:
		return self._operations

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return tuple(sorted({q for op in self._operations for q in op.qubits}))

	def __iter__(self) -> Iterator[Operation]:
		return iter(self._operations)

	def __len__(self) -> int:
		return len(self._operations)

	def __getitem__(self, index: int) -> Operation:
		return self._operations[index]

	def __bool__(self) -> bool:
		return bool(self._operations)

	def __eq__(self, other: object) -> bool:
		if not isinstance(other, Moment):
			return NotImplemented
		return frozenset(self._operations) == frozenset(other._operations)

	def __hash__(self) -> int:
		return hash(frozenset(self._operations))

	def __str__(self) -> str:
		return "[" + ", ".join(str(op) for op in self._operations) + "]"

	def __repr__(self) -> str:
		return f"vqubit.Moment({list(self._operations)!r})"

	def with_operation(self, operation: Operation) -> Moment:
		"""Return a new moment with ``operation`` added."""
		return Moment((*self._operations, operation))

	def without_operations_touching(self, qubits: Iterable[Qid]) -> Moment:
		"""Return a new moment with every operation on ``qubits`` removed."""
		qubits = set(qubits)
		return Moment(op for op in self._operations if not qubits.intersection(op.qubits))
