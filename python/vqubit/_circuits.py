"""Circuits: ordered sequences of moments.

A circuit is built by handing operations to :class:`Circuit`, in the order they
should run. vqubit packs them into the earliest moments that fit, so gates on
different qubits share a moment and gates on the same qubit are separated in
time, exactly as Cirq does with its ``EARLIEST`` insertion strategy.
"""

from __future__ import annotations

import enum
from collections.abc import Callable, Iterable, Iterator
from typing import Any

from ._ops import Moment, Operation
from ._qubits import Qid
from ._value import value_equality


class InsertStrategy(enum.Enum):
	"""How a new operation is placed among the existing moments."""

	EARLIEST = 0
	NEW = 1

	def __str__(self) -> str:
		return self.name


def _coerce_ops(content: Any) -> list[Operation]:
	"""Turn the argument of ``append``/``insert`` into a flat list of operations."""
	if isinstance(content, Operation):
		return [content]
	if isinstance(content, Moment):
		return list(content.operations)
	if isinstance(content, (str, bytes)):
		raise TypeError(f"cannot append {type(content).__name__} to a circuit")
	if isinstance(content, Iterable):
		operations: list[Operation] = []
		for item in content:
			if isinstance(item, Operation):
				operations.append(item)
			elif isinstance(item, Moment):
				# A moment (or a circuit, which iterates as moments) unfolds into
				# its operations instead of being silently dropped.
				operations.extend(item.operations)
			else:
				raise TypeError(f"cannot append {type(item).__name__} to a circuit")
		return operations
	raise TypeError(f"cannot append {type(content).__name__} to a circuit")


class AbstractCircuit:
	"""The base class of :class:`Circuit` and :class:`FrozenCircuit`."""

	@property
	def moments(self) -> tuple[Any, ...]:
		"""The moments of the circuit, in time order."""
		raise NotImplementedError


class Circuit(AbstractCircuit):
	"""A mutable, time-ordered sequence of moments.

	>>> from vqubit import CNOT, H, LineQubit
	>>> q0, q1 = LineQubit.range(2)
	>>> circuit = Circuit(H(q0), CNOT(q0, q1))
	>>> len(circuit)
	2
	>>> len(circuit.all_qubits())
	2
	"""

	def __init__(
		self,
		*contents: Any,
		strategy: InsertStrategy = InsertStrategy.EARLIEST,
		tags: Any = (),
	) -> None:
		self._moments: list[Moment] = []
		self._tags = tuple(tags)
		for content in contents:
			self.append(content, strategy=strategy)

	# -- construction ---------------------------------------------------------

	def _split_index(self, operation: Operation) -> int:
		"""Return the earliest moment index at which ``operation`` may be placed."""
		qubits = operation.qubits
		for i in range(len(self._moments) - 1, -1, -1):
			moment_qubits = self._moments[i].qubits_set
			if any(q in moment_qubits for q in qubits):
				return i + 1
		return 0

	def _place(self, operation: Operation, split: int, strategy: InsertStrategy) -> None:
		if strategy is InsertStrategy.NEW:
			self._moments.insert(split, Moment([operation]))
			return
		if split < len(self._moments):
			self._moments[split] = self._moments[split].with_operation(operation)
		else:
			self._moments.append(Moment([operation]))

	def append(self, content: Any, strategy: InsertStrategy = InsertStrategy.EARLIEST) -> Circuit:
		"""Append operations (or a moment) to the circuit."""
		if isinstance(content, Moment):
			for op in content.operations:
				self.append(op, strategy=strategy)
			return self
		for operation in _coerce_ops(content):
			split = len(self._moments) if strategy is InsertStrategy.NEW else self._split_index(operation)
			self._place(operation, split, strategy)
		return self

	def insert(self, index: int, content: Any, strategy: InsertStrategy = InsertStrategy.EARLIEST) -> Circuit:
		"""Insert operations at ``index``, respecting dependencies."""
		if index < 0:
			index += len(self._moments) + 1
		if not 0 <= index <= len(self._moments):
			raise IndexError(f"insert index {index} out of range")
		if isinstance(content, Moment):
			for op in content.operations:
				self.insert(index, op, strategy=strategy)
			return self
		for operation in _coerce_ops(content):
			split = max(index, self._split_index(operation))
			if strategy is InsertStrategy.NEW and split == index:
				self._place(operation, split, InsertStrategy.NEW)
			else:
				self._place(operation, split, strategy)
		return self

	@classmethod
	def from_moments(cls, *moments: Any, tags: Any = ()) -> Circuit:
		"""Build a circuit from explicit moments, keeping their boundaries.

		``Circuit.from_moments(m1, m2)`` is Cirq's calling form; a single
		iterable of moments is accepted too, since that is what the transformers
		in this library build.

		>>> from vqubit import H, LineQubit
		>>> q0, q1 = LineQubit.range(2)
		>>> len(Circuit.from_moments(Moment([H(q0)]), Moment([H(q1)])))
		2
		"""
		if len(moments) == 1 and not isinstance(moments[0], Moment):
			moments = tuple(moments[0])
		for moment in moments:
			if not isinstance(moment, Moment):
				raise TypeError(f"from_moments expects Moment instances, got {type(moment).__name__}")
		circuit = cls(tags=tags)
		# Assign the list directly: appending operation by operation would let two
		# independent moments collapse into one.
		circuit._moments = [Moment(list(moment.operations)) for moment in moments]
		return circuit

	# -- access ---------------------------------------------------------------

	@property
	def moments(self) -> tuple[Moment, ...]:
		"""The moments of the circuit, in time order.

		>>> from vqubit import H, LineQubit
		>>> len(Circuit(H(LineQubit(0))).moments)
		1
		"""
		return tuple(self._moments)

	def all_operations(self) -> Iterator[Operation]:
		"""Every operation in the circuit, moment by moment."""
		for moment in self._moments:
			yield from moment.operations

	def all_qubits(self) -> tuple[Qid, ...]:
		"""The qubits the circuit uses, sorted."""
		return tuple(sorted({q for moment in self._moments for q in moment.qubits_set}))

	def findall_operations(self, predicate: Callable[[Operation], bool]) -> Iterator[tuple[int, Operation]]:
		"""Yield ``(moment_index, operation)`` for every operation matching ``predicate``."""
		for index, moment in enumerate(self._moments):
			for operation in moment.operations:
				if predicate(operation):
					yield index, operation

	def findall_operations_with_gate_type(self, gate_type: type) -> Iterator[tuple[int, Operation]]:
		"""Yield ``(moment_index, operation)`` for operations whose gate is ``gate_type``."""
		return self.findall_operations(lambda operation: isinstance(operation.gate, gate_type))

	def unitary(self) -> Any:
		"""Return the full unitary matrix over the sorted circuit qubits."""
		from ._protocols import unitary

		return unitary(self)

	def to_text_diagram(
		self,
		*,
		use_unicode_characters: bool = True,
		transpose: bool = False,
		include_tags: Any = True,
		precision: int | None = 3,
		qubit_order: Any = None,
	) -> str:
		"""Return an ASCII diagram of the circuit.

		``qubit_order`` decides the order of the lines (any iterable of qubits, or
		a :class:`~vqubit.QubitOrder`); the other arguments are accepted for
		compatibility and only ``use_unicode_characters`` and ``precision``
		change the output today: the renderer is always one line per qubit.
		"""
		from ._diagram import to_text_diagram

		if qubit_order is not None and not isinstance(qubit_order, (list, tuple)):
			order = getattr(qubit_order, "order_for", None)
			if order is not None:
				qubit_order = order(self.all_qubits())
		del transpose, include_tags
		return to_text_diagram(
			self,
			qubit_order=qubit_order,
			use_unicode_characters=use_unicode_characters,
			precision=precision,
		)

	@property
	def tags(self) -> tuple[Any, ...]:
		"""The tags carried by the circuit, empty when it has none."""
		return self._tags

	def with_tags(self, *tags: Any) -> Circuit:
		"""Return a copy of the circuit carrying ``tags``, as Cirq does."""
		return Circuit.from_moments(self._moments, tags=(*self._tags, *tags))

	def __add__(self, other: Any) -> Circuit:
		"""Concatenate a circuit, or append an operation as ``cirq.Circuit`` does."""
		if isinstance(other, Circuit):
			return Circuit.from_moments([*self._moments, *other._moments], tags=self._tags)
		if isinstance(other, (Moment, Operation, list, tuple)):
			result = Circuit.from_moments(self._moments, tags=self._tags)
			result.append(other)
			return result
		return NotImplemented

	def __radd__(self, other: Any) -> Circuit:
		if other == 0:  # ``sum([circuit, ...])``
			return self
		if isinstance(other, Circuit):
			return Circuit.from_moments([*other._moments, *self._moments], tags=self._tags)
		if isinstance(other, (Moment, Operation, list, tuple)):
			result = Circuit.from_moments(self._moments, tags=self._tags)
			result.append(other)
			return result
		return NotImplemented

	def __mul__(self, repetitions: int) -> Circuit:
		if not isinstance(repetitions, int) or isinstance(repetitions, bool):
			return NotImplemented
		return Circuit.from_moments([*self._moments] * repetitions, tags=self._tags)

	def __rmul__(self, repetitions: int) -> Circuit:
		return self.__mul__(repetitions)

	def __iter__(self) -> Iterator[Moment]:
		return iter(self._moments)

	def __len__(self) -> int:
		return len(self._moments)

	def __bool__(self) -> bool:
		return bool(self._moments)

	def __getitem__(self, index: int | slice) -> Moment | Circuit:
		if isinstance(index, slice):
			return Circuit.from_moments(self._moments[index])
		return self._moments[index]

	def __reversed__(self) -> Iterator[Moment]:
		return reversed(self._moments)

	def __eq__(self, other: object) -> bool:
		if not isinstance(other, Circuit):
			return NotImplemented
		if len(self._moments) != len(other._moments):
			return False
		return all(a == b for a, b in zip(self._moments, other._moments, strict=True))

	def __hash__(self) -> int:
		return hash(tuple(self._moments))

	def __str__(self) -> str:
		from ._diagram import to_text_diagram

		return to_text_diagram(self)

	def __repr__(self) -> str:
		return f"vqubit.Circuit({list(self._moments)!r})"


class FrozenCircuit(Circuit):
	"""An immutable, hashable circuit.

	The moments are fixed at construction; the mutation methods from
	:class:`Circuit` are disabled.
	"""

	def __init__(
		self,
		*contents: Any,
		strategy: InsertStrategy = InsertStrategy.EARLIEST,
		tags: Any = (),
	) -> None:
		built = Circuit(*contents, strategy=strategy, tags=tags)
		self._moments = list(built._moments)
		self._tags = built._tags

	@classmethod
	def from_moments(cls, *moments: Any, tags: Any = ()) -> FrozenCircuit:
		"""Build a frozen circuit from explicit moments, keeping their boundaries."""
		built = Circuit.from_moments(*moments, tags=tags)
		frozen = cls(tags=tags)
		frozen._moments = list(built._moments)
		return frozen

	def append(self, content: Any, strategy: InsertStrategy = InsertStrategy.EARLIEST) -> Circuit:
		"""Return a new frozen circuit with the operations appended."""
		raise TypeError("FrozenCircuit is immutable")

	def insert(self, index: int, content: Any, strategy: InsertStrategy = InsertStrategy.EARLIEST) -> Circuit:
		"""Return a new frozen circuit with the operations inserted at ``index``."""
		raise TypeError("FrozenCircuit is immutable")


@value_equality
class CircuitOperation(Operation):
	"""A circuit used as a single operation inside another circuit.

	The wrapped circuit acts on its own qubits; ``repetitions`` runs it several
	times in a row.
	"""

	def __init__(
		self,
		circuit: Circuit,
		*,
		repetitions: int = 1,
		param_resolver: Any = None,
	) -> None:
		if param_resolver is not None:
			from ._protocols import resolve_parameters

			circuit = resolve_parameters(circuit, param_resolver)
		if not isinstance(circuit, Circuit):
			raise TypeError("CircuitOperation needs a Circuit")
		if repetitions < 1:
			raise ValueError("repetitions must be at least 1")
		self._circuit = circuit
		self._repetitions = int(repetitions)

	@property
	def circuit(self) -> Circuit:
		"""The circuit this operation runs."""
		return self._circuit

	@property
	def repetitions(self) -> int:
		"""How many times the circuit is run."""
		return self._repetitions

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits the wrapped circuit uses."""
		return self._circuit.all_qubits()

	def _num_qubits_(self) -> int:
		return len(self.qubits)

	def _unitary_(self) -> Any:
		from ._linalg import matmul

		matrix = self._circuit.unitary()
		result = matrix
		for _ in range(self._repetitions - 1):
			result = matmul(matrix, result)
		return result

	def repeat(self, repetitions: int) -> CircuitOperation:
		"""Return a copy that runs the circuit ``repetitions`` more times."""
		return CircuitOperation(self._circuit, repetitions=self._repetitions * int(repetitions))

	def _circuit_diagram_info_(self) -> tuple[str, ...]:
		return ("CircuitOperation",) * len(self.qubits)

	def _is_measurement_(self) -> bool:
		return any(operation._is_measurement_() for operation in self._circuit.all_operations())

	def _measurement_key_obj_(self) -> Any:
		return None

	def _is_parameterized_(self) -> bool:
		from ._protocols import is_parameterized

		return is_parameterized(self._circuit)

	def _resolve_parameters_(self, resolver: Any) -> CircuitOperation:
		from ._protocols import resolve_parameters

		return CircuitOperation(resolve_parameters(self._circuit, resolver), repetitions=self._repetitions)

	def _value_equality_values_(self) -> Any:
		return (self._circuit, self._repetitions)

	def __str__(self) -> str:
		return f"CircuitOperation({self._circuit!r}, repetitions={self._repetitions})"

	def __repr__(self) -> str:
		return f"vqubit.CircuitOperation({self._circuit!r}, repetitions={self._repetitions})"
