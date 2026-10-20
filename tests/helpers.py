"""Shared helpers for the test suite."""

from __future__ import annotations

import random
from collections.abc import Sequence
from typing import Any


def assert_matrix_close(actual: Any, expected: Any, tol: float = 1e-9) -> None:
	assert len(actual) == len(expected), f"row count {len(actual)} != {len(expected)}"
	for i, (row_a, row_e) in enumerate(zip(actual, expected, strict=True)):
		assert len(row_a) == len(row_e), f"column count in row {i}: {len(row_a)} != {len(row_e)}"
		for j, (a, e) in enumerate(zip(row_a, row_e, strict=True)):
			assert abs(a - e) <= tol, f"matrix[{i}][{j}]: {a!r} != {e!r}"


def assert_vector_close(actual: Sequence[complex], expected: Sequence[complex], tol: float = 1e-9) -> None:
	assert len(actual) == len(expected)
	for i, (a, e) in enumerate(zip(actual, expected, strict=True)):
		assert abs(a - e) <= tol, f"vector[{i}]: {a!r} != {e!r}"


def identity_matrix() -> list[list[complex]]:
	return [[1.0 + 0j, 0j], [0j, 1.0 + 0j]]


def random_unitary(dim: int, rng: random.Random) -> list[list[complex]]:
	"""Return a Haar-random ``dim``-by-``dim`` unitary, by Gram-Schmidt."""
	columns: list[list[complex]] = []
	for _ in range(dim):
		vector = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(dim)]
		for column in columns:
			overlap = sum(column[k].conjugate() * vector[k] for k in range(dim))
			vector = [vector[k] - overlap * column[k] for k in range(dim)]
		norm = sum(abs(z) ** 2 for z in vector) ** 0.5
		columns.append([z / norm for z in vector])
	return [[columns[j][i] for j in range(dim)] for i in range(dim)]


def up_to_phase(actual: Sequence[Sequence[complex]], expected: Sequence[Sequence[complex]], tol: float = 1e-6) -> bool:
	"""Report whether two matrices are equal up to a single global phase."""
	phase = None
	for i in range(len(expected)):
		for j in range(len(expected)):
			if abs(expected[i][j]) > 1e-9:
				phase = actual[i][j] / expected[i][j]
				break
		if phase is not None:
			break
	if phase is None:
		return False
	return all(
		abs(actual[i][j] - phase * expected[i][j]) < tol for i in range(len(expected)) for j in range(len(expected))
	)


__all__ = [
	"assert_matrix_close",
	"assert_vector_close",
	"identity_matrix",
	"random_unitary",
	"up_to_phase",
]
