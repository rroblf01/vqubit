"""Protocols: functions that ask a value what it can do.

These mirror Cirq's protocols. ``unitary(gate)`` returns the gate's matrix,
``num_qubits(value)`` its width, ``inverse(value)`` the undoing of it, and so
on. They accept gates, operations and circuits wherever that makes sense.
"""

from __future__ import annotations

from typing import Any

from ._circuits import Circuit
from ._gates import Gate, MeasurementGate
from ._linalg import Matrix, apply_gate, identity, matmul
from ._ops import GateOperation, Moment, Operation
from ._qubits import Qid


def _as_matrix(unitary_like: Any) -> Matrix | None:
	if unitary_like is NotImplemented or unitary_like is None:
		return None
	return [list(row) for row in unitary_like]


def unitary(value: Any) -> Matrix:
	"""Return the unitary matrix of ``value``, or raise ``TypeError``.

	For a gate, the matrix is in the gate's own qubit order. For an operation,
	it is in the order of ``operation.qubits``. For a circuit, it is the full
	matrix over the sorted circuit qubits.
	"""
	if isinstance(value, Circuit):
		return _circuit_unitary(value)
	if isinstance(value, (GateOperation, Operation)):
		matrix = _as_matrix(value.gate._unitary_())
		if matrix is None:
			raise TypeError(f"{value!r} has no unitary")
		return matrix
	if isinstance(value, Gate):
		matrix = _as_matrix(value._unitary_())
		if matrix is None:
			raise TypeError(f"{value!r} has no unitary")
		return matrix
	raise TypeError(f"{value!r} has no unitary")


def has_unitary(value: Any) -> bool:
	"""Report whether ``value`` has a unitary matrix."""
	try:
		unitary(value)
	except TypeError:
		return False
	return True


def _circuit_unitary(circuit: Circuit) -> Matrix:
	qubits = circuit.all_qubits()
	index = {q: i for i, q in enumerate(qubits)}
	dim = 1 << len(qubits)
	result = identity(dim)
	for moment in circuit.moments:
		result = matmul(_moment_unitary(moment, index, dim), result)
	return result


def _moment_unitary(moment: Moment, index: dict[Qid, int], dim: int) -> Matrix:
	columns: list[list[complex]] = []
	for j in range(dim):
		column: list[complex] = [0j] * dim
		column[j] = 1.0 + 0j
		for operation in moment.operations:
			if isinstance(operation.gate, MeasurementGate):
				raise TypeError("a circuit with measurements has no unitary")
			gate_matrix = _as_matrix(operation.gate._unitary_())
			if gate_matrix is None:
				raise TypeError(f"{operation.gate!r} has no unitary")
			targets = [index[q] for q in operation.qubits]
			column = apply_gate(column, gate_matrix, targets)
		columns.append(column)
	return [[columns[j][i] for j in range(dim)] for i in range(dim)]


def num_qubits(value: Any) -> int:
	"""Return the number of qubits ``value`` acts on."""
	if isinstance(value, Circuit):
		return len(value.all_qubits())
	if isinstance(value, (GateOperation, Operation)):
		return len(value.qubits)
	if isinstance(value, Gate):
		return value.num_qubits
	raise TypeError(f"cannot compute num_qubits of {value!r}")


def qid_shape(value: Any) -> tuple[int, ...]:
	"""Return the qid shape of ``value`` (all twos for qubits)."""
	return (2,) * num_qubits(value)


def inverse(value: Any) -> Any:
	"""Return the inverse of ``value``."""
	if isinstance(value, Gate):
		gate_inverse = value._inverse_() if hasattr(value, "_inverse_") else NotImplemented
		if gate_inverse is NotImplemented:
			raise TypeError(f"{value!r} has no inverse")
		return gate_inverse
	if isinstance(value, (GateOperation, Operation)):
		return inverse(value.gate).on(*value.qubits)
	if isinstance(value, Circuit):
		reversed_moments = [
			Moment([inverse(op) for op in reversed(moment.operations)]) for moment in reversed(value.moments)
		]
		return Circuit.from_moments(reversed_moments)
	raise TypeError(f"{value!r} has no inverse")


def pow(value: Any, exponent: float) -> Any:
	"""Return ``value`` raised to ``exponent``."""
	if isinstance(value, Gate):
		return value**exponent
	if isinstance(value, (GateOperation, Operation)):
		return (value.gate**exponent).on(*value.qubits)
	raise TypeError(f"{value!r} does not support powers")


def decompose(value: Any) -> list[Any]:
	"""Return a decomposition of ``value`` into simpler operations."""
	if isinstance(value, Circuit):
		return list(value.all_operations())
	if isinstance(value, (GateOperation, Operation)):
		return [value]
	if isinstance(value, Gate):
		return [value]
	raise TypeError(f"cannot decompose {value!r}")
