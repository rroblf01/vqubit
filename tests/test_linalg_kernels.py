"""Tests for the complex linear algebra kernels compiled from V.

The native module works on flat, interleaved complex arrays:
``[re0, im0, re1, im1, ...]``. These tests exercise the kernels directly, so a
regression in the V code fails here rather than in a simulator test.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import vqubit_native as native

SQRT_HALF = math.sqrt(0.5)

IDENTITY_2 = [
	1.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	1.0,
	0.0,
]
X = [
	0.0,
	0.0,
	1.0,
	0.0,
	1.0,
	0.0,
	0.0,
	0.0,
]
H = [
	SQRT_HALF,
	0.0,
	SQRT_HALF,
	0.0,
	SQRT_HALF,
	0.0,
	-SQRT_HALF,
	0.0,
]
CNOT_4 = [
	1.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	1.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	1.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	1.0,
	0.0,
	0.0,
	0.0,
]


def close(actual: Sequence[float], expected: Sequence[float], tol: float = 1e-9) -> bool:
	return len(actual) == len(expected) and all(abs(a - e) <= tol for a, e in zip(actual, expected, strict=True))


def test_matmul_identity() -> None:
	assert close(native.matmul(X, IDENTITY_2, 2), X)


def test_matmul_x_squared_is_identity() -> None:
	assert close(native.matmul(X, X, 2), IDENTITY_2)


def test_matmul_h_squared_is_identity() -> None:
	assert close(native.matmul(H, H, 2), IDENTITY_2)


def test_matmul_empty_dimension() -> None:
	assert native.matmul(X, X, 0) == []


def test_kron_dimension() -> None:
	assert len(native.kron(X, 2, H, 2)) == 2 * 4 * 4


def test_kron_empty_dimension() -> None:
	assert native.kron(X, 0, H, 2) == []
	assert native.kron(X, 2, H, 0) == []


def test_apply_gate_hadamard_on_first_qubit() -> None:
	state = [0.0] * 8
	state[0] = 1.0
	expected = [0.0] * 8
	expected[0] = SQRT_HALF
	expected[4] = SQRT_HALF
	assert close(native.apply_gate(state, H, 2, [0]), expected)


def test_apply_gate_hadamard_on_second_qubit() -> None:
	state = [0.0] * 8
	state[0] = 1.0
	expected = [0.0] * 8
	expected[0] = SQRT_HALF
	expected[2] = SQRT_HALF
	assert close(native.apply_gate(state, H, 2, [1]), expected)


def test_apply_gate_bell_state() -> None:
	state = [0.0] * 8
	state[0] = 1.0
	after_h = native.apply_gate(state, H, 2, [0])
	bell = native.apply_gate(after_h, CNOT_4, 2, [0, 1])
	expected = [0.0] * 8
	expected[0] = SQRT_HALF
	expected[6] = SQRT_HALF
	assert close(bell, expected)


def test_apply_gate_no_targets_is_identity() -> None:
	state = [0.0] * 8
	state[0] = 1.0
	assert native.apply_gate(state, [], 3, []) == state


def test_dagger_is_conjugate_transpose() -> None:
	x_dagger = native.dagger(X, 2)
	assert close(x_dagger, X)
	assert close(native.dagger(H, 2), H)


def test_add_scaled() -> None:
	assert close(native.add_scaled(X, X, 2, 1.0), [2.0 * v for v in X])


def test_add_scaled_empty_dimension() -> None:
	assert native.add_scaled(X, X, 0, 1.0) == []
