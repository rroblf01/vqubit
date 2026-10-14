"""Complex linear algebra helpers.

Matrices are represented as nested Python lists of :class:`complex`, which is
what a Cirq user expects to see from ``vqubit.unitary``. The heavy operations
(matrix multiply, tensor product, apply and conjugate transpose) are delegated
to the native ``vqubit_native`` kernels, which work on flat, interleaved arrays
of ``float``: ``[re0, im0, re1, im1, ...]``.
"""

from __future__ import annotations

from collections.abc import Sequence

import vqubit_native as _native

Matrix = list[list[complex]]
Vector = list[complex]


def to_flat(matrix: Sequence[Sequence[complex]]) -> list[float]:
	"""Flatten a complex matrix into the native interleaved representation."""
	flat: list[float] = []
	for row in matrix:
		for z in row:
			flat.append(float(z.real))
			flat.append(float(z.imag))
	return flat


def from_flat(flat: Sequence[float], dim: int) -> Matrix:
	"""Rebuild a ``dim``-by-``dim`` complex matrix from the flat buffer."""
	return [[complex(flat[2 * (i * dim + j)], flat[2 * (i * dim + j) + 1]) for j in range(dim)] for i in range(dim)]


def vector_to_flat(state: Sequence[complex]) -> list[float]:
	flat: list[float] = []
	for z in state:
		flat.append(float(z.real))
		flat.append(float(z.imag))
	return flat


def vector_from_flat(flat: Sequence[float]) -> Vector:
	return [complex(flat[2 * i], flat[2 * i + 1]) for i in range(len(flat) // 2)]


def identity(dim: int) -> Matrix:
	return [[1.0 + 0j if i == j else 0j for j in range(dim)] for i in range(dim)]


def matmul(a: Sequence[Sequence[complex]], b: Sequence[Sequence[complex]]) -> Matrix:
	dim = len(a)
	return from_flat(_native.matmul(to_flat(a), to_flat(b), dim), dim)


def kron(a: Sequence[Sequence[complex]], b: Sequence[Sequence[complex]]) -> Matrix:
	a_dim = len(a)
	b_dim = len(b)
	return from_flat(_native.kron(to_flat(a), a_dim, to_flat(b), b_dim), a_dim * b_dim)


def dagger(a: Sequence[Sequence[complex]]) -> Matrix:
	dim = len(a)
	return from_flat(_native.dagger(to_flat(a), dim), dim)


def apply_gate(state: Sequence[complex], gate: Sequence[Sequence[complex]], targets: Sequence[int]) -> Vector:
	"""Apply ``gate`` to ``targets`` of the state vector, in big-endian order."""
	num_qubits = (len(state)).bit_length() - 1
	num_targets = len(targets)
	dim = len(gate)
	if dim != 1 << num_targets:
		raise ValueError(f"gate dimension {dim} does not match {num_targets} targets")
	flat_gate = to_flat(gate)
	return vector_from_flat(_native.apply_gate(vector_to_flat(state), flat_gate, num_qubits, list(targets)))


def is_unitary(matrix: Sequence[Sequence[complex]], tol: float = 1e-9) -> bool:
	dim = len(matrix)
	product = from_flat(_native.matmul(to_flat(matrix), to_flat(dagger(matrix)), dim), dim)
	for i in range(dim):
		for j in range(dim):
			expected = 1.0 + 0j if i == j else 0j
			if abs(product[i][j] - expected) > tol:
				return False
	return True


def transpose(matrix: Sequence[Sequence[complex]]) -> Matrix:
	"""Return the transpose (no conjugation)."""
	rows = len(matrix)
	cols = len(matrix[0]) if rows else 0
	return [[complex(matrix[i][j]) for i in range(rows)] for j in range(cols)]


def conjugate(matrix: Sequence[Sequence[complex]]) -> Matrix:
	"""Return the elementwise complex conjugate."""
	return [[complex(z).conjugate() for z in row] for row in matrix]


def det(matrix: Sequence[Sequence[complex]]) -> complex:
	"""Return the determinant, by Gaussian elimination with partial pivoting."""
	a = [[complex(z) for z in row] for row in matrix]
	n = len(a)
	value = 1 + 0j
	for col in range(n):
		pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
		if abs(a[pivot][col]) < 1e-14:
			return 0j
		if pivot != col:
			a[col], a[pivot] = a[pivot], a[col]
			value = -value
		value *= a[col][col]
		inv = 1.0 / a[col][col]
		for r in range(col + 1, n):
			factor = a[r][col] * inv
			if factor != 0:
				for c in range(col, n):
					a[r][c] -= factor * a[col][c]
	return value


def inverse(matrix: Sequence[Sequence[complex]]) -> Matrix:
	"""Return the inverse, by Gauss-Jordan elimination."""
	n = len(matrix)
	a = [[complex(z) for z in row] for row in matrix]
	m = identity(n)
	for col in range(n):
		pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
		if abs(a[pivot][col]) < 1e-14:
			raise ValueError("matrix is singular")
		a[col], a[pivot] = a[pivot], a[col]
		m[col], m[pivot] = m[pivot], m[col]
		inv = 1.0 / a[col][col]
		a[col] = [z * inv for z in a[col]]
		m[col] = [z * inv for z in m[col]]
		for r in range(n):
			if r == col:
				continue
			factor = a[r][col]
			if factor != 0:
				for c in range(n):
					a[r][c] -= factor * a[col][c]
					m[r][c] -= factor * m[col][c]
	return m


def real_symmetric_eig(
	matrix: Sequence[Sequence[float]],
	*,
	max_sweeps: int = 100,
	tol: float = 1e-14,
) -> tuple[list[float], list[list[float]]]:
	"""Diagonalize a real symmetric matrix with the cyclic Jacobi method.

	Returns the eigenvalues and a real orthogonal matrix whose *columns* are
	the corresponding eigenvectors, so ``A = V diag(values) V^T``.
	"""
	n = len(matrix)
	a = [[float(matrix[i][j]) for j in range(n)] for i in range(n)]
	v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
	for _ in range(max_sweeps):
		off = sum(a[i][j] * a[i][j] for i in range(n) for j in range(n) if i != j)
		if off < tol:
			break
		for p in range(n):
			for q in range(p + 1, n):
				if abs(a[p][q]) < tol:
					continue
				theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
				t = (
					1.0 / (abs(theta) + (theta * theta + 1) ** 0.5)
					if theta >= 0
					else -1.0 / (abs(theta) + (theta * theta + 1) ** 0.5)
				)
				c = 1.0 / (t * t + 1) ** 0.5
				s = t * c
				for k in range(n):
					akp, akq = a[k][p], a[k][q]
					a[k][p] = c * akp - s * akq
					a[k][q] = s * akp + c * akq
				for k in range(n):
					apk, aqk = a[p][k], a[q][k]
					a[p][k] = c * apk - s * aqk
					a[q][k] = s * apk + c * aqk
				for k in range(n):
					vkp, vkq = v[k][p], v[k][q]
					v[k][p] = c * vkp - s * vkq
					v[k][q] = s * vkp + c * vkq
	values = [a[i][i] for i in range(n)]
	return values, v
