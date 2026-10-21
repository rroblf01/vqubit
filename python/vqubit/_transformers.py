"""Circuit transformers.

A transformer takes a circuit and returns another circuit. These are the
housekeeping ones: remove the empty moments, drop operations that do nothing,
expand composite gates, and re-order operations into strata so that similar
work sits together.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import vqubit_native as _native

from ._circuits import Circuit
from ._gates import MatrixGate
from ._linalg import matmul, to_flat
from ._ops import Moment, Operation, TaggedOperation
from ._optimizers import single_qubit_matrix_to_phxz
from ._protocols import decompose
from ._qubits import Qid

Predicate = Callable[[Operation], bool]


def drop_empty_moments(circuit: Circuit) -> Circuit:
	"""Return ``circuit`` without its empty moments."""
	return Circuit.from_moments([moment for moment in circuit.moments if len(moment) > 0])


def _matrix_is_identity(matrix: list[list[complex]], atol: float) -> bool:
	return bool(_native.is_identity_matrix(to_flat(matrix), len(matrix), atol))


def drop_negligible_operations(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Return ``circuit`` without operations whose unitary is the identity."""
	kept: list[Operation] = []
	for operation in circuit.all_operations():
		if operation.gate._is_measurement_():
			kept.append(operation)
			continue
		matrix = operation.gate._unitary_()
		if matrix is None:
			kept.append(operation)
			continue
		if not _matrix_is_identity([list(row) for row in matrix], atol):
			kept.append(operation)
	return Circuit(kept)


def expand_composite(circuit: Circuit) -> Circuit:
	"""Return ``circuit`` with every decomposable gate expanded."""
	return Circuit(decompose(circuit))


def merge_single_qubit_moments_to_phxz(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Merge runs of single-qubit gates into one ``PhasedXZ`` gate per qubit.

	Consecutive moments that hold only single-qubit gates are combined; on each
	qubit the gates are multiplied together and replaced by the
	:class:`~vqubit.PhasedXZGate` that implements the same unitary, up to phase.
	"""
	result = Circuit()
	pending: dict[Qid, list[list[complex]]] = {}

	def flush() -> None:
		if not pending:
			return
		operations: list[Operation] = []
		for qubit, matrix in pending.items():
			if not _matrix_is_identity(matrix, atol):
				operations.append(single_qubit_matrix_to_phxz(matrix, atol=atol).on(qubit))
		if operations:
			result.append(operations)
		pending.clear()

	for moment in circuit.moments:
		mergeable = len(moment) > 0 and all(
			len(operation.qubits) == 1 and not operation.gate._is_measurement_() for operation in moment.operations
		)
		if not mergeable:
			flush()
			compiled: list[Operation] = []
			for operation in moment.operations:
				if len(operation.qubits) != 1 or operation.gate._is_measurement_():
					compiled.append(operation)
					continue
				matrix = operation.gate._unitary_()
				if matrix is None:
					compiled.append(operation)
				elif not _matrix_is_identity([list(row) for row in matrix], atol):
					compiled.append(
						single_qubit_matrix_to_phxz([list(row) for row in matrix], atol=atol).on(operation.qubits[0])
					)
			if compiled:
				result.append(compiled)
			continue
		for operation in moment.operations:
			matrix = operation.gate._unitary_()
			if matrix is None:
				flush()
				result.append(operation)
				continue
			qubit = operation.qubits[0]
			current = pending.get(qubit)
			pending[qubit] = (
				[list(row) for row in matrix] if current is None else matmul([list(row) for row in matrix], current)
			)
	flush()
	return result


def _as_predicate(category: Any) -> Predicate:
	if isinstance(category, type):
		return lambda operation: isinstance(operation.gate, category)
	if callable(category):
		return category
	raise TypeError(f"cannot use {category!r} as a category")


def stratified_circuit(circuit: Circuit, *positional: Any, categories: Any = ()) -> Circuit:
	"""Return ``circuit`` with its operations ordered by category.

	The operations matching the first category come first, then those matching
	the second, and so on; the original order is kept within each group, and
	anything matching no category is placed last. Moments are not preserved.
	The categories are taken positionally or by name, and a single iterable is
	unpacked, so all of Cirq's calling forms work.
	"""
	category_list: tuple[Any, ...] = tuple(positional) + tuple(categories)
	if len(category_list) == 1 and isinstance(category_list[0], (list, tuple)):
		category_list = tuple(category_list[0])
	predicates = [_as_predicate(category) for category in category_list]
	operations = list(circuit.all_operations())
	result = Circuit()
	for predicate in predicates:
		for operation in operations:
			if predicate(operation):
				result.append(operation)
	for operation in operations:
		if not any(predicate(operation) for predicate in predicates):
			result.append(operation)
	return result


def align_left(circuit: Circuit) -> Circuit:
	"""Return a circuit with operations packed as early as possible.

	This is exactly how vqubit builds circuits, so it is the identity here and
	exists to mirror Cirq's API.
	"""
	return Circuit.from_moments(circuit.moments)


def align_right(circuit: Circuit) -> Circuit:
	"""Return a circuit with operations packed as late as possible."""
	reversed_operations = list(reversed(list(circuit.all_operations())))
	packed = Circuit(reversed_operations)
	moments = [Moment(list(reversed(moment.operations))) for moment in reversed(packed.moments)]
	return Circuit.from_moments(moments)


__all__ = [
	"align_left",
	"align_right",
	"drop_empty_moments",
	"drop_negligible_operations",
	"expand_composite",
	"map_moments",
	"map_operations",
	"merge_k_qubit_unitaries",
	"merge_single_qubit_moments_to_phxz",
	"stratified_circuit",
	"synchronize_terminal_measurements",
	"toggle_tags",
]


def map_operations(circuit: Circuit, map_func: Callable[[Operation], Any]) -> Circuit:
	"""Return a circuit with every operation replaced by ``map_func``'s result.

	``map_func`` may return an operation or an iterable of operations.
	"""
	result = Circuit()
	for operation in circuit.all_operations():
		mapped = map_func(operation)
		if isinstance(mapped, Operation):
			result.append(mapped)
		else:
			for item in mapped:
				result.append(item)
	return result


def map_moments(circuit: Circuit, map_func: Callable[[Moment], Moment]) -> Circuit:
	"""Return a circuit with every moment replaced by ``map_func``'s result."""
	return Circuit.from_moments(map_func(moment) for moment in circuit.moments)


def toggle_tags(circuit: Circuit, *tags: Any) -> Circuit:
	"""Return ``circuit`` with ``tags`` added to every operation."""

	def add_tags(operation: Operation) -> Operation:
		if isinstance(operation, TaggedOperation):
			return TaggedOperation(operation.sub_operation, *operation.tags, *tags)
		return TaggedOperation(operation, *tags)

	return map_operations(circuit, add_tags)


def merge_k_qubit_unitaries(circuit: Circuit, k: int = 1, *, atol: float = 1e-8) -> Circuit:
	"""Merge runs of unitary gates on at most ``k`` qubits into ``MatrixGate``s.

	Consecutive unitary operations on the same qubits are composed into a single
	``MatrixGate``. A measurement or an operation on more than ``k`` qubits ends
	the runs it touches.
	"""
	result = Circuit()
	pending: dict[tuple[Qid, ...], list[list[complex]]] = {}

	def flush() -> None:
		if not pending:
			return
		operations: list[Operation] = []
		for qubits, matrix in pending.items():
			if not _matrix_is_identity(matrix, atol):
				operations.append(MatrixGate(matrix).on(*qubits))
		if operations:
			result.append(operations)
		pending.clear()

	for operation in circuit.all_operations():
		gate = operation.gate
		matrix = None if gate._is_measurement_() else gate._unitary_()
		if len(operation.qubits) <= k and matrix is not None:
			key = operation.qubits
			composed = [list(row) for row in matrix]
			pending[key] = composed if key not in pending else matmul(composed, pending[key])
			continue
		flush()
		result.append(operation)
	flush()
	return result


def synchronize_terminal_measurements(circuit: Circuit) -> Circuit:
	"""Move all terminal measurements into a single final moment.

	A measurement is terminal when nothing but measurements happen after it on
	its qubits. If any measurement is not terminal, the circuit is returned
	unchanged.
	"""
	operations = list(circuit.all_operations())
	measurements = [operation for operation in operations if operation.gate._is_measurement_()]
	if not measurements:
		return Circuit.from_moments(circuit.moments)
	for index, operation in enumerate(operations):
		if not operation.gate._is_measurement_():
			continue
		qubits = set(operation.qubits)
		for later in operations[index + 1 :]:
			if not later.gate._is_measurement_() and qubits.intersection(later.qubits):
				return Circuit.from_moments(circuit.moments)
	result = Circuit(operation for operation in operations if not operation.gate._is_measurement_())
	try:
		result.append(Moment(measurements))
	except ValueError:
		return Circuit.from_moments(circuit.moments)
	return result
