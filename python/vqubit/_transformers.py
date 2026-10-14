"""Circuit transformers.

A transformer takes a circuit and returns another circuit. These are the
housekeeping ones: remove the empty moments, drop operations that do nothing,
expand composite gates, and re-order operations into strata so that similar
work sits together.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from ._circuits import Circuit
from ._ops import Moment, Operation
from ._protocols import decompose

Predicate = Callable[[Operation], bool]


def drop_empty_moments(circuit: Circuit) -> Circuit:
	"""Return ``circuit`` without its empty moments."""
	return Circuit.from_moments([moment for moment in circuit.moments if len(moment) > 0])


def _matrix_is_identity(matrix: list[list[complex]], atol: float) -> bool:
	dim = len(matrix)
	for i in range(dim):
		for j in range(dim):
			expected = 1.0 + 0j if i == j else 0j
			if abs(matrix[i][j] - expected) > atol:
				return False
	return True


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


def _as_predicate(category: Any) -> Predicate:
	if isinstance(category, type):
		return lambda operation: isinstance(operation.gate, category)
	if callable(category):
		return category
	raise TypeError(f"cannot use {category!r} as a category")


def stratified_circuit(circuit: Circuit, *categories: Any) -> Circuit:
	"""Return ``circuit`` with its operations ordered by category.

	The operations matching the first category come first, then those matching
	the second, and so on; the original order is kept within each group, and
	anything matching no category is placed last. Moments are not preserved.
	"""
	if len(categories) == 1 and isinstance(categories[0], (list, tuple)):
		categories = tuple(categories[0])
	predicates = [_as_predicate(category) for category in categories]
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
	"stratified_circuit",
]
