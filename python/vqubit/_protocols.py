"""Protocols: functions that ask a value what it can do.

These mirror Cirq's protocols. ``unitary(gate)`` returns the gate's matrix,
``num_qubits(value)`` its width, ``inverse(value)`` the undoing of it, and so
on. They accept gates, operations and circuits wherever that makes sense.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from itertools import product
from typing import Any

from ._circuits import Circuit
from ._gates import Gate, MatrixGate, MeasurementGate, ZPowGate
from ._linalg import Matrix, apply_gate, dagger, identity, kron, matmul
from ._ops import GateOperation, Moment, Operation
from ._parameter import ParamResolver, is_symbolic, resolve_value
from ._pauli import Pauli, PauliString, PauliSum
from ._qubits import LineQubit, Qid


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
	"""Return a decomposition of ``value`` into simpler operations.

	Gates that can be expanded (the Toffoli, the Fredkin, ...) are recursively
	replaced by their decomposition; everything else is returned unchanged.
	"""
	if isinstance(value, Circuit):
		result: list[Any] = []
		for operation in value.all_operations():
			result.extend(_decompose_operation(operation))
		return result
	if isinstance(value, (GateOperation, Operation)):
		return _decompose_operation(value)
	if isinstance(value, Gate):
		return [value]
	raise TypeError(f"cannot decompose {value!r}")


def decompose_once(value: Any) -> list[Any]:
	"""Return one level of decomposition of ``value``."""
	if isinstance(value, (GateOperation, Operation)):
		return _decompose_level(value)
	if isinstance(value, Gate):
		return [value]
	if isinstance(value, Circuit):
		result: list[Any] = []
		for operation in value.all_operations():
			result.extend(_decompose_level(operation))
		return result
	raise TypeError(f"cannot decompose {value!r}")


def _decompose_level(operation: Operation) -> list[Any]:
	gate = operation.gate
	expanded = gate._decompose_(operation.qubits)
	if expanded is None:
		return [operation]
	return list(expanded)


def _decompose_operation(operation: Operation, depth: int = 0) -> list[Any]:
	if depth >= 10:
		return [operation]
	gate = operation.gate
	expanded = gate._decompose_(operation.qubits)
	if expanded is None:
		return [operation]
	result: list[Any] = []
	for child in expanded:
		result.extend(_decompose_operation(child, depth + 1))
	return result


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


def pauli_expansion(value: Any, *, atol: float = 1e-8) -> PauliSum:
	"""Expand ``value`` as a weighted sum of Pauli strings.

	Every Pauli string over the value's qubits is tested; only those whose
	coefficient is larger than ``atol`` are kept.
	"""
	matrix = unitary(value)
	dim = len(matrix)
	num_qubits = dim.bit_length() - 1
	if isinstance(value, (GateOperation, Operation)):
		qubits: tuple[Qid, ...] = value.qubits
	else:
		qubits = tuple(LineQubit(i) for i in range(num_qubits))

	terms: dict[PauliString, complex] = {}
	for assignment in product((Pauli.I, Pauli.X, Pauli.Y, Pauli.Z), repeat=num_qubits):
		probe: Matrix = [[1 + 0j]]
		for pauli in assignment:
			probe = kron(probe, pauli.matrix())
		trace = 0j
		for i in range(dim):
			for j in range(dim):
				trace += probe[i][j] * matrix[j][i]
		coefficient = trace / dim
		if abs(coefficient) > atol:
			terms[PauliString(dict(zip(qubits, assignment, strict=True)))] = coefficient
	return PauliSum(terms)


def _embed_single(gate_matrix: Matrix, target: int, num_qubits: int) -> Matrix:
	dim = 1 << num_qubits
	columns: list[list[complex]] = []
	for column_index in range(dim):
		column = [0j] * dim
		column[column_index] = 1.0 + 0j
		columns.append(apply_gate(column, gate_matrix, [target]))
	return [[columns[j][i] for j in range(dim)] for i in range(dim)]


def phase_by(value: Any, phase_exponent: float, qubit_index: int) -> MatrixGate:
	"""Conjugate ``value`` by a phase on one of its qubits.

	The result is ``Z**phase_exponent`` applied to ``qubit_index``, ``value``,
	and the inverse phase, returned as a :class:`~vqubit.MatrixGate`.
	"""
	matrix = unitary(value)
	dim = len(matrix)
	num_qubits = dim.bit_length() - 1
	if not 0 <= qubit_index < num_qubits:
		raise ValueError(f"qubit_index {qubit_index} out of range for {num_qubits} qubits")
	rotation = ZPowGate(exponent=phase_exponent)._unitary_()
	if rotation is None:
		raise TypeError("phase_by needs a numeric phase exponent")
	embedding = _embed_single(rotation, qubit_index, num_qubits)
	return MatrixGate(matmul(matmul(embedding, matrix), dagger(embedding)))


def apply_unitary(value: Any, state_vector: Sequence[complex]) -> list[complex]:
	"""Apply the unitary of ``value`` to ``state_vector``.

	The vector is expressed over the value's own qubits, so its length must be
	the dimension of the unitary.
	"""
	matrix = unitary(value)
	dim = len(matrix)
	if len(state_vector) != dim:
		raise ValueError(f"state_vector has length {len(state_vector)}, expected {dim}")
	column = [complex(x) for x in state_vector]
	return [sum(matrix[i][j] * column[j] for j in range(dim)) for i in range(dim)]
