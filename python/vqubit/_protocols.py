"""Protocols: functions that ask a value what it can do.

These mirror Cirq's protocols. ``unitary(gate)`` returns the gate's matrix,
``num_qubits(value)`` its width, ``inverse(value)`` the undoing of it, and so
on. They accept gates, operations and circuits wherever that makes sense.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from typing import Any

import vqubit_native as _native

from . import _compat
from ._circuits import Circuit, CircuitOperation
from ._compat import asarray
from ._gates import Gate, MatrixGate, MeasurementGate, ZPowGate
from ._linalg import Matrix, dagger, from_flat, identity, matmul, to_flat, vector_from_flat, vector_to_flat
from ._measurement_key import MeasurementKey
from ._ops import GateOperation, Moment, Operation
from ._parameter import ParamResolver, is_symbolic, resolve_value
from ._pauli import Pauli, PauliString, PauliSum
from ._qubits import LineQubit, Qid

_PAULI_BY_DIGIT = (Pauli.I, Pauli.X, Pauli.Y, Pauli.Z)


def _as_matrix(unitary_like: Any) -> Matrix | None:
	if unitary_like is NotImplemented or unitary_like is None:
		return None
	return [list(row) for row in unitary_like]


def unitary(value: Any) -> Matrix:
	"""Return the unitary matrix of ``value``, or raise ``TypeError``.

	For a gate, the matrix is in the gate's own qubit order. For an operation,
	it is in the order of ``operation.qubits``. For a circuit, it is the full
	matrix over the sorted circuit qubits. When numpy is installed the result is
	a ``numpy.ndarray``.
	"""
	if isinstance(value, Circuit):
		return asarray(_circuit_unitary(value))
	if isinstance(value, CircuitOperation):
		return asarray(value._unitary_())
	if isinstance(value, (GateOperation, Operation)):
		matrix = _as_matrix(value.gate._unitary_())
		if matrix is None:
			matrix = _decomposed_unitary(value)
		if matrix is None:
			raise TypeError(f"{value!r} has no unitary")
		return asarray(matrix)
	if isinstance(value, Gate):
		matrix = _as_matrix(value._unitary_())
		if matrix is None:
			matrix = _decomposed_unitary(value)
		if matrix is None:
			raise TypeError(f"{value!r} has no unitary")
		return asarray(matrix)
	if isinstance(value, PauliString):
		return asarray(value.matrix())
	raise TypeError(f"{value!r} has no unitary")


def _decomposed_unitary(value: Any) -> Matrix | None:
	"""Build the unitary of a gate/operation from its decomposition, if any."""
	subject = value.gate if isinstance(value, (GateOperation, Operation)) else value
	decompose_method = getattr(subject, "_decompose_", None)
	if decompose_method is None:
		return None
	qubits = (
		value.qubits if isinstance(value, (GateOperation, Operation)) else tuple(LineQubit.range(subject.num_qubits))
	)
	try:
		operations = decompose_method(qubits)
	except TypeError:
		return None
	if operations is None:
		return None
	return _circuit_unitary(Circuit(list(operations)))


def has_unitary(value: Any) -> bool:
	"""Report whether ``value`` has a unitary matrix."""
	try:
		unitary(value)
	except TypeError:
		return False
	return True


def _circuit_unitary(circuit: Circuit) -> Any:
	if _compat.HAS_NUMPY:
		return _numpy_circuit_unitary(circuit)
	qubits = circuit.all_qubits()
	index = {q: i for i, q in enumerate(qubits)}
	dim = 1 << len(qubits)
	result = identity(dim)
	for moment in circuit.moments:
		result = matmul(_moment_unitary(moment, index, dim), result)
	return result


def _numpy_circuit_unitary(circuit: Circuit) -> Any:
	"""Build the circuit unitary with numpy (BLAS) when it is installed."""
	import numpy as np

	qubits = circuit.all_qubits()
	index = {q: i for i, q in enumerate(qubits)}
	num_qubits = len(qubits)
	dim = 1 << num_qubits
	result = np.eye(dim, dtype=np.complex128)
	for moment in circuit.moments:
		for operation in moment.operations:
			if isinstance(operation.gate, MeasurementGate):
				raise TypeError("a circuit with measurements has no unitary")
			matrix = _as_matrix(operation.gate._unitary_())
			if matrix is None:
				raise TypeError(f"{operation.gate!r} has no unitary")
			targets = [index[q] for q in operation.qubits]
			result = _numpy_apply_to_unitary(result, np.asarray(matrix, dtype=np.complex128), num_qubits, targets)
	return result


def _numpy_apply_to_unitary(U: Any, gate_matrix: Any, num_qubits: int, targets: Sequence[int]) -> Any:
	"""Left-multiply the running unitary by a gate on ``targets``."""
	import numpy as np

	k = len(targets)
	gate_tensor = gate_matrix.reshape((2,) * k + (2,) * k)
	tensor = U.reshape((2,) * (2 * num_qubits))
	other = [axis for axis in range(2 * num_qubits) if axis not in targets]
	permuted = tensor.transpose(other + list(targets))
	contracted = np.tensordot(
		gate_tensor, permuted, axes=(list(range(k, 2 * k)), list(range(2 * num_qubits - k, 2 * num_qubits)))
	)
	moved = np.moveaxis(contracted, list(range(k)), list(targets))
	dim = 1 << num_qubits
	return moved.reshape(dim, dim)


def _moment_unitary(moment: Moment, index: dict[Qid, int], dim: int) -> Matrix:
	num_qubits = dim.bit_length() - 1
	result = identity(dim)
	for operation in moment.operations:
		if isinstance(operation.gate, MeasurementGate):
			raise TypeError("a circuit with measurements has no unitary")
		gate_matrix = _as_matrix(operation.gate._unitary_())
		if gate_matrix is None:
			raise TypeError(f"{operation.gate!r} has no unitary")
		targets = [index[q] for q in operation.qubits]
		embedded = from_flat(_native.embed_operator(to_flat(gate_matrix), len(gate_matrix), targets, num_qubits), dim)
		result = matmul(embedded, result)
	return result


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
	return bool(_native.matrix_allclose(to_flat(ab), to_flat(ba), len(ab), atol))


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

	coefficients = _native.pauli_expansion(to_flat(matrix), num_qubits)
	terms: dict[PauliString, complex] = {}
	for combo in range(1 << (2 * num_qubits)):
		coefficient = complex(coefficients[2 * combo], coefficients[2 * combo + 1])
		if abs(coefficient) <= atol:
			continue
		assignment = tuple(_PAULI_BY_DIGIT[(combo >> (2 * (num_qubits - 1 - q))) & 3] for q in range(num_qubits))
		terms[PauliString(dict(zip(qubits, assignment, strict=True)))] = coefficient
	return PauliSum(terms)


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
	embedding = from_flat(
		_native.embed_operator(to_flat([list(row) for row in rotation]), 2, [qubit_index], num_qubits), dim
	)
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
	return vector_from_flat(_native.matvec(to_flat(matrix), vector_to_flat(state_vector), dim))


def kraus(value: Any) -> list[Matrix]:
	"""Return the Kraus operators of ``value``.

	A unitary gate has exactly one Kraus operator, its unitary matrix; a channel
	has as many as it needs.
	"""
	if isinstance(value, (GateOperation, Operation)):
		return kraus(value.gate)
	if isinstance(value, Gate):
		operators = value._kraus_()
		if operators is not None:
			return [asarray([list(row) for row in operator]) for operator in operators]
		matrix = _as_matrix(value._unitary_())
		if matrix is not None:
			return [asarray(matrix)]
		raise TypeError(f"{value!r} has no Kraus operators")
	raise TypeError(f"{value!r} has no Kraus operators")


def has_kraus(value: Any) -> bool:
	"""Report whether ``value`` has Kraus operators."""
	try:
		kraus(value)
	except TypeError:
		return False
	return True


def _matrices_approx_equal(a: Matrix, b: Matrix, atol: float) -> bool:
	return len(a) == len(b) and bool(_native.matrix_allclose(to_flat(a), to_flat(b), len(a), atol))


def _matrices_phase_equal(a: Matrix, b: Matrix, atol: float) -> bool:
	if len(a) != len(b):
		return False
	phase = None
	for i in range(len(b)):
		for j in range(len(b)):
			if abs(b[i][j]) > 1e-9:
				phase = a[i][j] / b[i][j]
				break
		if phase is not None:
			break
	if phase is None:
		return True
	scaled = [[phase * entry for entry in row] for row in b]
	return bool(_native.matrix_allclose(to_flat(a), to_flat(scaled), len(a), atol))


def approx_eq(val1: Any, val2: Any, *, atol: float = 1e-8) -> bool:
	"""Report whether two values are equal, comparing unitaries approximately."""
	try:
		u1 = unitary(val1)
		u2 = unitary(val2)
	except TypeError:
		return bool(val1 == val2)
	return _matrices_approx_equal(u1, u2, atol)


def equal_up_to_global_phase(val1: Any, val2: Any, *, atol: float = 1e-8) -> bool:
	"""Report whether two values are equal up to a global phase."""
	try:
		u1 = unitary(val1)
		u2 = unitary(val2)
	except TypeError:
		return bool(val1 == val2)
	return _matrices_phase_equal(u1, u2, atol)


def measurement_key_obj(value: Any) -> MeasurementKey | None:
	"""Return the :class:`~vqubit.MeasurementKey` a measurement names, if any."""
	if isinstance(value, (GateOperation, Operation)):
		return value.gate._measurement_key_obj_()
	if isinstance(value, Gate):
		return value._measurement_key_obj_()
	return None


def measurement_key_objs(value: Any) -> frozenset[MeasurementKey]:
	"""Return the measurement keys ``value`` depends on."""
	if isinstance(value, Circuit):
		operations: Iterable[Any] = value.all_operations()
	elif isinstance(value, Moment):
		operations = value.operations
	elif isinstance(value, (GateOperation, Operation, Gate)):
		operations = [value]
	else:
		return frozenset()
	keys = (measurement_key_obj(operation) for operation in operations)
	return frozenset(key for key in keys if key is not None)


def control_keys(value: Any) -> frozenset[MeasurementKey]:
	"""Return the classical keys that control ``value`` (none yet)."""
	return frozenset()


def with_measurement_key_mapping(value: Any, mapping: Any) -> Any:
	"""Return ``value`` with its measurement keys renamed by ``mapping``."""

	def remap(key: str) -> str:
		for source, target in dict(mapping).items():
			name = source.name if isinstance(source, MeasurementKey) else str(source)
			if name == key:
				return target.name if isinstance(target, MeasurementKey) else str(target)
		return key

	if isinstance(value, MeasurementGate):
		return MeasurementGate(remap(value.key), num_qubits=value.num_qubits, invert_mask=value.invert_mask)
	if isinstance(value, (GateOperation, Operation)):
		return with_measurement_key_mapping(value.gate, mapping).on(*value.qubits)
	if isinstance(value, Circuit):
		return Circuit([with_measurement_key_mapping(op, mapping) for op in value.all_operations()])
	return value


def trace_distance_bound(value: Any) -> float:
	"""Return an upper bound on the trace distance of ``value`` from the identity."""
	hook = getattr(value, "_trace_distance_bound_", None)
	if hook is not None:
		return float(hook())
	try:
		matrix = unitary(value)
	except TypeError:
		return 1.0
	dim = len(matrix)
	trace = sum(matrix[i][i] for i in range(dim)) / dim
	return math.sqrt(max(0.0, 1.0 - abs(trace) ** 2))
