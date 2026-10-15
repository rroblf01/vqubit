"""Quantum states: state vectors, density matrices and their measures.

These mirror Cirq's ``qis`` helpers. A state vector is a list of ``complex`` of
length ``2**n``; a density matrix is a nested list of ``complex`` of shape
``dim`` by ``dim``. Everything here is dependency-free.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence
from typing import Any

import vqubit_native as _native

from ._circuits import Circuit
from ._compat import asarray, asvector, is_numpy_array, to_nested_list
from ._linalg import Matrix, dagger, matmul

StateVectorLike = Sequence[complex]


def _is_density(value: Any) -> bool:
	if is_numpy_array(value):
		return value.ndim == 2
	return isinstance(value, Sequence) and len(value) > 0 and isinstance(value[0], Sequence)


def _dim_of(qid_shape: Sequence[int] | None, length: int) -> int:
	if qid_shape is not None:
		dim = 1
		for size in qid_shape:
			dim *= size
		return dim
	return length


def to_valid_state_vector(
	state: Any,
	qid_shape: Sequence[int],
	*,
	atol: float = 1e-8,
) -> list[complex]:
	"""Return ``state`` as a normalized state vector, or raise ``ValueError``."""
	state = state.tolist() if is_numpy_array(state) else state
	dim = 1
	for size in qid_shape:
		dim *= size
	if isinstance(state, int) and not isinstance(state, bool):
		if not 0 <= state < dim:
			raise ValueError(f"state {state} out of range for dimension {dim}")
		vector = [0j] * dim
		vector[state] = 1.0 + 0j
		return vector
	vector = [complex(x) for x in state]
	if len(vector) != dim:
		raise ValueError(f"state has length {len(vector)}, expected {dim}")
	norm = sum(abs(x) ** 2 for x in vector)
	if abs(norm - 1.0) > atol:
		raise ValueError(f"state is not normalized (norm squared {norm})")
	return vector


def state_vector(state: Any, *, qid_shape: Sequence[int] | None = None) -> list[complex]:
	"""Return ``state`` as a state vector.

	An integer is a computational basis state; a sequence is taken as amplitudes
	and validated. A :class:`~vqubit.Circuit` is simulated instead.
	"""
	if isinstance(state, Circuit):
		from ._sim import Simulator

		return Simulator().simulate(state).state_vector()
	shape = qid_shape
	if shape is None:
		shape = (2,) * ((len(state)).bit_length() - 1) if not isinstance(state, int) else ()
	return to_valid_state_vector(state, shape)


def to_valid_density_matrix(
	matrix: Any,
	qid_shape: Sequence[int],
	*,
	atol: float = 1e-8,
) -> Matrix:
	"""Return ``matrix`` as a density matrix (a state vector becomes ``|psi><psi|``)."""
	matrix = to_nested_list(matrix) if is_numpy_array(matrix) else matrix
	if isinstance(matrix, int) and not isinstance(matrix, bool):
		vector = to_valid_state_vector(matrix, qid_shape, atol=atol)
		return asarray([[vector[i] * vector[j].conjugate() for j in range(len(vector))] for i in range(len(vector))])
	if _is_density(matrix):
		rows = [[complex(x) for x in row] for row in matrix]
		dim = len(rows)
		if any(len(row) != dim for row in rows):
			raise ValueError("density matrix must be square")
		trace = sum(rows[i][i] for i in range(dim))
		if abs(trace - 1.0) > atol:
			raise ValueError(f"density matrix trace is {trace}, expected 1")
		return asarray(rows)
	vector = to_valid_state_vector(matrix, qid_shape, atol=atol)
	return asarray([[vector[i] * vector[j].conjugate() for j in range(len(vector))] for i in range(len(vector))])


def density_matrix(matrix: Any, *, qid_shape: Sequence[int] | None = None) -> Matrix:
	"""Return ``matrix`` as a density matrix."""
	if isinstance(matrix, Circuit):
		from ._sim import DensityMatrixSimulator

		return DensityMatrixSimulator().simulate(matrix).density_matrix()
	shape = qid_shape
	if shape is None:
		if _is_density(matrix):
			shape = (2,) * (len(matrix).bit_length() - 1)
		elif isinstance(matrix, int):
			shape = ()
		else:
			shape = (2,) * (len(matrix).bit_length() - 1)
	return to_valid_density_matrix(matrix, shape)


def bloch_vector_from_state_vector(state: StateVectorLike, index: int) -> tuple[float, float, float]:
	"""Return the Bloch vector ``(x, y, z)`` of qubit ``index`` of ``state``."""
	if is_numpy_array(state):
		import numpy as _np

		array = _np.asarray(state, dtype=_np.complex128)
		num_qubits = array.size.bit_length() - 1
		if not 0 <= index < num_qubits:
			raise ValueError(f"index {index} out of range for {num_qubits} qubits")
		tensor = array.reshape((2,) * num_qubits)
		zero = _np.take(tensor, 0, axis=index)
		one = _np.take(tensor, 1, axis=index)
		product = _np.sum(zero.conj() * one)
		x = 2 * product.real
		y = 2 * product.imag
		z = _np.sum(_np.abs(zero) ** 2) - _np.sum(_np.abs(one) ** 2)
		return (float(x), float(y), float(z))
	vector = list(state)
	num_qubits = len(vector).bit_length() - 1
	if not 0 <= index < num_qubits:
		raise ValueError(f"index {index} out of range for {num_qubits} qubits")
	x = 0.0
	y = 0.0
	z = 0.0
	for i, amplitude in enumerate(vector):
		if (i >> (num_qubits - 1 - index)) & 1 == 0:
			j = i | (1 << (num_qubits - 1 - index))
			other = vector[j]
			product = amplitude.conjugate() * other
			x += 2 * product.real
			y += 2 * product.imag
			z += abs(amplitude) ** 2 - abs(other) ** 2
	return (x, y, z)


def partial_trace(tensor: Sequence[Sequence[complex]], keep_indices: Sequence[int]) -> Matrix:
	"""Trace out the qubits of a density matrix that are not in ``keep_indices``."""
	rows = [[complex(x) for x in row] for row in tensor]
	dim = len(rows)
	num_qubits = dim.bit_length() - 1
	keep = sorted(set(keep_indices))
	for index in keep:
		if not 0 <= index < num_qubits:
			raise ValueError(f"keep index {index} out of range for {num_qubits} qubits")
	kept_dim = 1 << len(keep)
	reduced = [[0j] * kept_dim for _ in range(kept_dim)]
	for i in range(dim):
		kept_i = 0
		for index in keep:
			kept_i = (kept_i << 1) | ((i >> (num_qubits - 1 - index)) & 1)
		for j in range(dim):
			kept_j = 0
			for index in keep:
				kept_j = (kept_j << 1) | ((j >> (num_qubits - 1 - index)) & 1)
			# The traced bits must agree for the term to contribute.
			if _traced_bits_equal(i, j, keep, num_qubits):
				reduced[kept_i][kept_j] += rows[i][j]
	return asarray(reduced)


def _traced_bits_equal(i: int, j: int, keep: Sequence[int], num_qubits: int) -> bool:
	for index in range(num_qubits):
		if index in keep:
			continue
		if ((i >> (num_qubits - 1 - index)) & 1) != ((j >> (num_qubits - 1 - index)) & 1):
			return False
	return True


def _hermitian_eig(matrix: Matrix) -> tuple[list[float], Matrix]:
	"""Eigenvalues and eigenvectors (columns) of a Hermitian matrix.

	Uses the real symmetric embedding ``[[A, -B], [B, A]]`` for ``H = A + iB``,
	whose eigenvalues are the eigenvalues of ``H`` each doubled, and the native
	Jacobi eigensolver.
	"""
	n = len(matrix)
	size = 2 * n
	big = [0.0] * (size * size)
	for i in range(n):
		for j in range(n):
			re = matrix[i][j].real
			im = matrix[i][j].imag
			big[i * size + j] = re
			big[i * size + n + j] = -im
			big[(n + i) * size + j] = im
			big[(n + i) * size + n + j] = re
	result = _native.jacobi_eig(big, size)
	values = [result[k] for k in range(size)]
	vectors = result[size:]
	# Each eigenvalue of the embedding appears twice; sort so the pair is
	# adjacent, then take one of each pair.
	order = sorted(range(size), key=lambda index: values[index])
	eigenvalues: list[float] = []
	complex_vectors: Matrix = []
	for k in range(0, size, 2):
		index = order[k]
		eigenvalues.append(values[index])
		column = [complex(vectors[i * size + index], vectors[(n + i) * size + index]) for i in range(n)]
		complex_vectors.append(column)
	return eigenvalues, complex_vectors


def fidelity(state1: Any, state2: Any) -> float:
	"""Return the fidelity between two states (state vectors or density matrices)."""
	if _is_density(state1) and _is_density(state2):
		return _mixed_fidelity(list(state1), list(state2))
	return _mixed_fidelity(density_matrix(state1), density_matrix(state2))


def _mixed_fidelity(rho: Matrix, sigma: Matrix) -> float:
	dim = len(rho)
	root = _hermitian_sqrt(rho, dim)
	product = matmul(matmul(root, sigma), root)
	values, _ = _hermitian_eig(product)
	total = sum(math.sqrt(max(value, 0.0)) for value in values)
	return max(min(total * total, 1.0), 0.0)


def _hermitian_sqrt(matrix: Matrix, dim: int) -> Matrix:
	values, vectors = _hermitian_eig(matrix)
	result = [[0j] * dim for _ in range(dim)]
	for k, value in enumerate(values):
		root = math.sqrt(max(value, 0.0))
		vector = vectors[k]
		for i in range(dim):
			for j in range(dim):
				result[i][j] += root * vector[i] * vector[j].conjugate()
	return asarray(result)


def von_neumann_entropy(density: Matrix, *, base: float = 2.0) -> float:
	"""Return the von Neumann entropy ``-Tr(rho log rho)`` of a density matrix."""
	values, _ = _hermitian_eig(density)
	entropy = 0.0
	for value in values:
		if value > 1e-12:
			entropy -= value * (math.log(value) / math.log(base))
	return entropy


def random_state_vector(num_qubits: int, *, seed: int | None = None) -> list[complex]:
	"""Return a Haar-random state vector of ``num_qubits`` qubits."""
	rng = random.Random(seed)
	dim = 1 << num_qubits
	vector = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(dim)]
	norm = math.sqrt(sum(abs(x) ** 2 for x in vector))
	return asvector([x / norm for x in vector])


def random_density_matrix(num_qubits: int, *, seed: int | None = None) -> Matrix:
	"""Return a random density matrix (the Ginibre ensemble)."""
	rng = random.Random(seed)
	dim = 1 << num_qubits
	g = [[complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(dim)] for _ in range(dim)]
	product = matmul(g, dagger(g))
	trace = sum(product[i][i].real for i in range(dim))
	return asarray([[product[i][j] / trace for j in range(dim)] for i in range(dim)])


def kraus_channel(kraus_operators: Sequence[Sequence[Sequence[complex]]]) -> list[Matrix]:
	"""Validate and return the Kraus operators of a channel."""
	operators = [[[complex(x) for x in row] for row in op] for op in kraus_operators]
	if not operators:
		raise ValueError("a channel needs at least one Kraus operator")
	dim = len(operators[0])
	total = [[0j] * dim for _ in range(dim)]
	for operator in operators:
		if len(operator) != dim:
			raise ValueError("Kraus operators must share a shape")
		product = matmul(dagger(operator), operator)
		for i in range(dim):
			for j in range(dim):
				total[i][j] += product[i][j]
	for i in range(dim):
		for j in range(dim):
			expected = 1.0 if i == j else 0.0
			if abs(total[i][j] - expected) > 1e-8:
				raise ValueError("Kraus operators are not trace-preserving")
	return operators


__all__ = [
	"bloch_vector_from_state_vector",
	"density_matrix",
	"fidelity",
	"kraus_channel",
	"partial_trace",
	"random_density_matrix",
	"random_state_vector",
	"state_vector",
	"to_valid_density_matrix",
	"to_valid_state_vector",
	"von_neumann_entropy",
]
