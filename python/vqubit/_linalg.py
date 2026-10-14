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
	product = _native.matmul(to_flat(matrix), to_flat(dagger(matrix)), dim)
	return bool(_native.is_identity_matrix(product, dim, tol))


def transpose(matrix: Sequence[Sequence[complex]]) -> Matrix:
	"""Return the transpose (no conjugation)."""
	rows = len(matrix)
	cols = len(matrix[0]) if rows else 0
	return [[complex(matrix[i][j]) for i in range(rows)] for j in range(cols)]


def conjugate(matrix: Sequence[Sequence[complex]]) -> Matrix:
	"""Return the elementwise complex conjugate."""
	return [[complex(z).conjugate() for z in row] for row in matrix]


def det(matrix: Sequence[Sequence[complex]]) -> complex:
	"""Return the determinant, in the native kernel."""
	value = _native.matrix_det(to_flat(matrix), len(matrix))
	return complex(value[0], value[1])
