"""Tests for the public ``vqubit.linalg`` helpers."""

from __future__ import annotations

import math
import random

import pytest

import vqubit
from tests.helpers import assert_matrix_close, assert_vector_close, random_unitary, up_to_phase

X = [[0j, 1 + 0j], [1 + 0j, 0j]]
I2 = [[1 + 0j, 0j], [0j, 1 + 0j]]


def test_predicates() -> None:
	assert vqubit.is_unitary(X)
	assert vqubit.is_hermitian(X)
	assert vqubit.is_normal(X)
	assert vqubit.is_diagonal(I2)
	assert not vqubit.is_diagonal(X)
	assert vqubit.is_orthogonal([[0, 1], [-1, 0]])
	assert vqubit.is_special_orthogonal([[0, 1], [-1, 0]])
	assert vqubit.is_special_unitary(X) is False
	assert vqubit.is_special_unitary(vqubit.to_special(X))
	assert vqubit.matrix_commutes(X, X)
	assert not vqubit.matrix_commutes(X, [[1, 0], [0, -1]])
	assert vqubit.is_cptp(kraus_ops=[X])


def test_all_near_zero() -> None:
	assert vqubit.all_near_zero([0, 1e-12])
	assert not vqubit.all_near_zero([1e-3])
	assert vqubit.all_near_zero_mod([0, 2 * math.pi], 2 * math.pi)


def test_global_phase() -> None:
	phased = [[z * (1j) for z in row] for row in X]
	assert vqubit.allclose_up_to_global_phase(X, phased)
	assert not vqubit.allclose_up_to_global_phase(X, I2)
	a2, b2 = vqubit.match_global_phase(X, phased)
	assert_matrix_close(a2, b2)


def test_combinators() -> None:
	assert_matrix_close(vqubit.dot(X, X), I2)
	assert_matrix_close(vqubit.kron(X, X), vqubit.unitary(vqubit.XX))
	assert_matrix_close(vqubit.block_diag(X, I2), [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
	controlled = vqubit.kron_with_controls(vqubit.CONTROL_TAG, X)
	assert_matrix_close(controlled, vqubit.unitary(vqubit.CNOT))


def test_operator_spaces() -> None:
	assert vqubit.hilbert_schmidt_inner_product(X, X) == pytest.approx(2)
	coefficients = vqubit.expand_matrix_in_orthogonal_basis(X, vqubit.PAULI_BASIS)
	assert coefficients["X"] == pytest.approx(1)
	assert coefficients["I"] == pytest.approx(0)
	assert_matrix_close(vqubit.matrix_from_basis_coefficients({"X": 1}, vqubit.PAULI_BASIS), X)
	assert vqubit.pow_pauli_combination(0, 1, 0, 0, 2) == pytest.approx((1, 0, 0, 0))
	assert set(vqubit.kron_bases(vqubit.PAULI_BASIS, repeat=2)) == {
		"I" + "I", "I" + "X", "I" + "Y", "I" + "Z",
		"X" + "I", "X" + "X", "X" + "Y", "X" + "Z",
		"Y" + "I", "Y" + "X", "Y" + "Y", "Y" + "Z",
		"Z" + "I", "Z" + "X", "Z" + "Y", "Z" + "Z",
	}  # fmt: skip


def test_transformations() -> None:
	assert_matrix_close(vqubit.reflection_matrix_pow(X, 0.5), vqubit.unitary(vqubit.X**0.5))
	assert vqubit.is_special_unitary(vqubit.to_special(X))
	assert up_to_phase(vqubit.to_special(X), X)
	state = vqubit.state_vector_kronecker_product([1, 0], [0, 1])
	assert list(state) == pytest.approx([0, 1, 0, 0])
	assert_vector_close(vqubit.apply_matrix_to_slices([0, 1, 0, 0], X, [1, 2]), [0, 0, 1, 0])
	identity = vqubit.targeted_left_multiply(X, [[1, 0], [0, 1]], [0])
	assert_matrix_close(identity, X)


def test_sub_state_vector() -> None:
	factored = vqubit.sub_state_vector([1, 0, 0, 0], [1])
	assert abs(factored[0]) == pytest.approx(1)
	assert vqubit.sub_state_vector([1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)], [0], default="x") == "x"


def test_diagonalize_real_symmetric() -> None:
	matrix = [[2, 1], [1, 2]]
	p = vqubit.diagonalize_real_symmetric_matrix(matrix)
	pt = [[p[j][i] for j in range(2)] for i in range(2)]
	assert vqubit.is_diagonal(vqubit.dot(pt, matrix, p))


def test_unitary_eig_and_map_eigenvalues() -> None:
	values, _basis = vqubit.unitary_eig(X)
	assert sorted(complex(v).real for v in values) == pytest.approx([-1, 1])
	assert_matrix_close(vqubit.map_eigenvalues(X, lambda v: v), X)


def test_axis_angle_roundtrip() -> None:
	rng = random.Random(0)
	for _ in range(10):
		unitary = random_unitary(2, rng)
		decomposition = vqubit.axis_angle(unitary)
		assert up_to_phase(decomposition.matrix(), unitary, tol=1e-6)


def test_kak_vector_and_decomposition() -> None:
	vector = vqubit.kak_vector(vqubit.unitary(vqubit.CNOT))
	assert abs(complex(vector[0])) == pytest.approx(math.pi / 4, abs=1e-6)
	# canonical constraints: 0 <= |z| <= y <= x <= pi/4
	abs_vector = sorted((abs(complex(c)) for c in vector), reverse=True)
	assert abs_vector[2] <= abs_vector[1] + 1e-9 <= abs_vector[0] + 1e-9 <= math.pi / 4 + 1e-9
	rng = random.Random(1)
	for _ in range(10):
		unitary = random_unitary(4, rng)
		decomposition = vqubit.kak_decomposition(unitary)
		assert up_to_phase(decomposition.matrix(), unitary, tol=1e-5)
		assert vqubit.num_cnots_required(unitary) in (0, 1, 2, 3)


def test_slice_for_qubits_equal_to() -> None:
	assert vqubit.slice_for_qubits_equal_to([1, 3], 0b11) == (slice(None), 1, slice(None), 1, Ellipsis)
	with pytest.raises(ValueError):
		vqubit.slice_for_qubits_equal_to([0], 1, big_endian_qureg_value=1)


def test_bidiagonalize_unitary() -> None:
	rng = random.Random(3)
	for _ in range(5):
		unitary = random_unitary(4, rng)
		left, diagonal, right = vqubit.bidiagonalize_unitary_with_special_orthogonals(unitary)
		assert vqubit.is_special_orthogonal(left, atol=1e-6)
		assert vqubit.is_special_orthogonal(right, atol=1e-6)
		assert up_to_phase(vqubit.dot(left, unitary, right), diagonal, tol=1e-5)


def test_kak_canonicalize_vector() -> None:
	decomposition = vqubit.kak_canonicalize_vector(0.1, 0.2, 0.3)
	x, y, z = decomposition.interaction_coefficients
	assert 0 <= abs(z) <= y + 1e-9 <= x + 1e-9 <= math.pi / 4 + 1e-9


def test_partial_trace_of_state_vector_as_mixture() -> None:
	bell = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	mixture = vqubit.partial_trace_of_state_vector_as_mixture(bell, [0])
	total = sum(probability for probability, _ in mixture)
	assert total == pytest.approx(1)


def test_kron_factor_and_so4() -> None:
	matrix = vqubit.kron(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.H))
	phase, a, b = vqubit.kron_factor_4x4_to_2x2s(matrix)
	assert up_to_phase(vqubit.kron(a, b), matrix, tol=1e-6)
	assert abs(phase) == pytest.approx(1, abs=1e-6)
