"""Shared helpers for the test suite."""

from __future__ import annotations

from collections.abc import Sequence


def assert_matrix_close(
	actual: Sequence[Sequence[complex]], expected: Sequence[Sequence[complex]], tol: float = 1e-9
) -> None:
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


__all__ = ["assert_matrix_close", "assert_vector_close", "identity_matrix"]
