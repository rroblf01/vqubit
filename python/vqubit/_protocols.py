"""Protocols: functions that ask a value what it can do.

These mirror Cirq's protocols. ``unitary(gate)`` returns the gate's matrix,
``num_qubits(value)`` its width, ``inverse(value)`` the undoing of it, and so
on. They accept gates, operations and circuits wherever that makes sense.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from ._circuits import Circuit
from ._gates import Gate, MeasurementGate
from ._linalg import Matrix, apply_gate, identity, matmul
from ._ops import GateOperation, Moment, Operation
from ._parameter import ParamResolver, is_symbolic, resolve_value
from ._pauli import PauliString
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
	if isinstance(value, PauliString):
		return value.matrix()
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
		gate_inverse = value._inverse_()
		if gate_inverse is None:
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


def is_measurement(value: Any) -> bool:
	"""Report whether ``value`` contains a measurement."""
	if isinstance(value, Circuit):
		return any(is_measurement(operation) for operation in value.all_operations())
	if isinstance(value, Moment):
		return any(is_measurement(operation) for operation in value.operations)
	if isinstance(value, (Gate, GateOperation, Operation)):
		return bool(value._is_measurement_())
	return False


def measurement_key_names(value: Any) -> set[str]:
	"""Return the classical keys whose values ``value`` depends on."""
	if isinstance(value, Circuit):
		operations: Iterable[Any] = value.all_operations()
	elif isinstance(value, Moment):
		operations = value.operations
	elif isinstance(value, (GateOperation, Operation, Gate)):
		operations = [value]
	else:
		return set()
	keys: set[str] = set()
	for operation in operations:
		if is_measurement(operation):
			key = operation._measurement_key_name_()
			if key is not None:
				keys.add(key)
	return keys


def measurement_keys_touched(value: Any) -> set[str]:
	"""Alias of :func:`measurement_key_names`, matching Cirq's name."""
	return measurement_key_names(value)


def _commutator_is_zero(a: Matrix, b: Matrix, atol: float) -> bool:
	ab = matmul(a, b)
	ba = matmul(b, a)
	dim = len(ab)
	return all(abs(ab[i][j] - ba[i][j]) <= atol for i in range(dim) for j in range(dim))


def commutes(a: Any, b: Any, atol: float = 1e-8) -> bool:
	"""Report whether ``a`` and ``b`` commute.

	Operations on disjoint qubits always commute. Otherwise the unitaries are
	compared, and anything that cannot be decided is reported as not commuting.
	"""
	qubits_a = a.qubits if isinstance(a, (GateOperation, Operation)) else None
	qubits_b = b.qubits if isinstance(b, (GateOperation, Operation)) else None

	if qubits_a is not None and qubits_b is not None:
		if set(qubits_a).isdisjoint(qubits_b):
			return True
		if tuple(qubits_a) == tuple(qubits_b):
			matrix_a = _as_matrix(a.gate._unitary_())
			matrix_b = _as_matrix(b.gate._unitary_())
			if matrix_a is not None and matrix_b is not None and len(matrix_a) == len(matrix_b):
				return _commutator_is_zero(matrix_a, matrix_b, atol)
		return False

	if isinstance(a, Gate) and isinstance(b, Gate) and a.num_qubits == b.num_qubits:
		matrix_a = _as_matrix(a._unitary_())
		matrix_b = _as_matrix(b._unitary_())
		if matrix_a is not None and matrix_b is not None and len(matrix_a) == len(matrix_b):
			return _commutator_is_zero(matrix_a, matrix_b, atol)

	return False


def definitely_commutes(a: Any, b: Any, atol: float = 1e-8) -> bool:
	"""Report whether ``a`` and ``b`` are known to commute."""
	return commutes(a, b, atol=atol)


def is_parameterized(value: Any) -> bool:
	"""Report whether ``value`` depends on a symbolic parameter."""
	if isinstance(value, Circuit):
		return any(is_parameterized(operation) for operation in value.all_operations())
	if isinstance(value, Moment):
		return any(is_parameterized(operation) for operation in value.operations)
	if isinstance(value, (GateOperation, Operation)):
		return bool(value.gate._is_parameterized_())
	if isinstance(value, Gate):
		return bool(value._is_parameterized_())
	return is_symbolic(value)


def resolve_parameters(value: Any, resolver: ParamResolver) -> Any:
	"""Substitute parameters in ``value`` using ``resolver``."""
	if isinstance(value, Circuit):
		moments = [Moment([resolve_parameters(op, resolver) for op in moment.operations]) for moment in value.moments]
		return Circuit.from_moments(moments)
	if isinstance(value, Moment):
		return Moment([resolve_parameters(op, resolver) for op in value.operations])
	if isinstance(value, (GateOperation, Operation)):
		return resolve_parameters(value.gate, resolver).on(*value.qubits)
	if isinstance(value, Gate):
		return value._resolve_parameters_(resolver)
	return resolve_value(value, resolver)
