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


def _identity(n: int) -> list[list[complex]]:
	return [[1 + 0j if i == j else 0j for j in range(n)] for i in range(n)]


def _transpose(matrix: list[list[complex]]) -> list[list[complex]]:
	return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]


def test_all_near_zero_nested_and_tolerance() -> None:
	assert vqubit.all_near_zero([[0, 0j], [1e-12, -2e-9]])
	assert not vqubit.all_near_zero([0, 1e-4])
	assert vqubit.all_near_zero([1e-9], atol=1e-8)
	assert not vqubit.all_near_zero([1e-7], atol=1e-8)


def test_all_near_zero_mod_period_and_negative() -> None:
	period = 2 * math.pi
	assert vqubit.all_near_zero_mod([0, period, -period, 3 * period], period)
	assert not vqubit.all_near_zero_mod(period / 2, period)
	assert vqubit.all_near_zero_mod([1.0, -1.0], 1.0)
	assert not vqubit.all_near_zero_mod(0.25, 1.0)


def test_predicates_reject_wrong_shapes() -> None:
	assert not vqubit.is_hermitian([[1, 0, 0], [0, 1, 0]])
	assert not vqubit.is_unitary([[1, 0, 0], [0, 1, 0]])
	assert not vqubit.is_normal([[1, 0, 0], [0, 1, 0]])
	assert not vqubit.matrix_commutes(_identity(2), [[1, 0, 0], [0, 1, 0]])
	assert not vqubit.is_orthogonal([[1j, 0], [0, 1j]])
	assert vqubit.is_diagonal([[1, 0], [0, 2]])
	assert not vqubit.is_diagonal([[1, 1e-6], [0, 2]])


def test_special_predicates_use_determinant() -> None:
	rotation = [[0, -1], [1, 0]]
	reflection = [[1, 0], [0, -1]]
	assert vqubit.is_orthogonal(rotation)
	assert vqubit.is_special_orthogonal(rotation)
	assert vqubit.is_orthogonal(reflection)
	assert not vqubit.is_special_orthogonal(reflection)
	special_unitary = [[1j, 0], [0, -1j]]
	assert vqubit.is_unitary(special_unitary)
	assert vqubit.is_special_unitary(special_unitary)
	assert not vqubit.is_special_unitary(vqubit.unitary(vqubit.X))


def test_is_cptp_accepts_complete_kraus_sets_only() -> None:
	assert vqubit.is_cptp(kraus_ops=[])
	decay = [[1, 0], [0, math.sqrt(0.5)]]
	flip = [[0, math.sqrt(0.5)], [0, 0]]
	assert vqubit.is_cptp(kraus_ops=[decay, flip])
	assert not vqubit.is_cptp(kraus_ops=[[[2, 0], [0, 1]]])


def test_match_global_phase_makes_matrices_equal() -> None:
	a = [[1j, 0], [0, -1j]]
	b = [[-1, 0], [0, 1]]
	matched_a, matched_b = vqubit.match_global_phase(a, b)
	assert_matrix_close(matched_a, matched_b)
	assert abs(matched_a[0][0]) == pytest.approx(1)
	assert vqubit.allclose_up_to_global_phase(a, b)
	# Mismatched shapes are returned unchanged.
	x, y = vqubit.match_global_phase(_identity(2), [[1, 0, 0], [0, 1, 0]])
	assert len(x[0]) == 2
	assert len(y[0]) == 3


def test_allclose_up_to_global_phase_rejects_different_matrices() -> None:
	assert not vqubit.allclose_up_to_global_phase(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.Z))
	assert not vqubit.allclose_up_to_global_phase(_identity(2), [[1, 0, 0], [0, 1, 0]])
	assert vqubit.allclose_up_to_global_phase([[0, 0], [0, 0]], [[0, 0], [0, 0]])


def test_dot_matrix_vector_and_single_argument() -> None:
	assert_vector_close(vqubit.dot(vqubit.unitary(vqubit.X), [1, 0]), [0, 1])
	assert_matrix_close(vqubit.dot(vqubit.unitary(vqubit.X)), vqubit.unitary(vqubit.X))
	assert_vector_close(vqubit.dot([1, 2]), [1, 2])
	with pytest.raises(ValueError):
		vqubit.dot()


def test_dot_chain_three_matrices() -> None:
	product = vqubit.dot(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.Y), vqubit.unitary(vqubit.Z))
	assert_matrix_close(product, [[1j, 0], [0, 1j]])


def test_kron_multiple_factors_and_scalars() -> None:
	x = vqubit.unitary(vqubit.X)
	product = vqubit.kron(x, vqubit.unitary(vqubit.Y), x)
	assert len(product) == 8
	assert len(product[0]) == 8
	assert sum(1 for row in product for z in row if abs(z) > 1e-9) == 8
	assert vqubit.kron(2, 3) == 6
	assert_matrix_close(vqubit.kron(), [[1]])


def test_kron_with_controls_replaces_nan_blocks_with_identity() -> None:
	controlled = vqubit.kron_with_controls(vqubit.CONTROL_TAG, vqubit.unitary(vqubit.Z))
	assert_matrix_close(controlled, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
	assert all(not (math.isnan(z.real) or math.isnan(z.imag)) for row in controlled for z in row)


def test_block_diag_empty_and_non_square() -> None:
	assert len(vqubit.block_diag()) == 0
	combined = vqubit.block_diag([[1]], [[2, 3], [4, 5]], [[6]])
	assert_matrix_close(
		combined,
		[[1, 0, 0, 0], [0, 2, 3, 0], [0, 4, 5, 0], [0, 0, 0, 6]],
	)
	with pytest.raises(ValueError, match="square"):
		vqubit.block_diag([[1, 2, 3], [4, 5, 6]])


def test_hilbert_schmidt_inner_product_of_paulis() -> None:
	assert vqubit.hilbert_schmidt_inner_product(_identity(2), _identity(2)) == pytest.approx(2)
	assert vqubit.hilbert_schmidt_inner_product(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.X)) == pytest.approx(2)
	assert abs(vqubit.hilbert_schmidt_inner_product(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.Y))) < 1e-9


def test_pauli_basis_expansion_round_trip() -> None:
	matrix = [[0.5, 0], [0, -0.5]]
	coefficients = vqubit.expand_matrix_in_orthogonal_basis(matrix, vqubit.PAULI_BASIS)
	assert coefficients["I"] == pytest.approx(0)
	assert coefficients["Z"] == pytest.approx(0.5)
	assert_matrix_close(vqubit.matrix_from_basis_coefficients(coefficients, vqubit.PAULI_BASIS), matrix)


def test_kron_bases_repeat_and_empty() -> None:
	assert vqubit.kron_bases() == {"": [[1 + 0j]]}
	assert set(vqubit.kron_bases(vqubit.PAULI_BASIS)) == {"I", "X", "Y", "Z"}
	pairs = vqubit.kron_bases(vqubit.PAULI_BASIS, repeat=2)
	assert len(pairs) == 16
	assert_matrix_close(
		pairs["XZ"],
		vqubit.kron(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.Z)),
	)


def test_pow_pauli_combination_known_cases() -> None:
	assert vqubit.pow_pauli_combination(0, 0, 0, 0, 0) == pytest.approx((1, 0, 0, 0))
	assert vqubit.pow_pauli_combination(0, 1, 0, 0, 1) == pytest.approx((0, 1, 0, 0))
	assert vqubit.pow_pauli_combination(0, 1, 0, 0, 2) == pytest.approx((1, 0, 0, 0))
	assert vqubit.pow_pauli_combination(0, 0, 0, 1, 3) == pytest.approx((0, 0, 0, 1))
	# (I + X)**2 = 2I + 2X
	assert vqubit.pow_pauli_combination(1, 1, 0, 0, 2) == pytest.approx((2, 2, 0, 0))


def test_reflection_matrix_pow_matches_gate_powers() -> None:
	x = vqubit.unitary(vqubit.X)
	for exponent in (0, 1, 2, 3):
		assert_matrix_close(vqubit.reflection_matrix_pow(x, exponent), vqubit.unitary(vqubit.X**exponent))
	assert_matrix_close(vqubit.reflection_matrix_pow(x, 0.5), vqubit.unitary(vqubit.X**0.5), tol=1e-12)
	assert_matrix_close(
		vqubit.reflection_matrix_pow(vqubit.unitary(vqubit.Z), 0.5),
		vqubit.unitary(vqubit.Z**0.5),
		tol=1e-12,
	)


def test_to_special_has_unit_determinant() -> None:
	special = vqubit.to_special(vqubit.unitary(vqubit.H))
	determinant = special[0][0] * special[1][1] - special[0][1] * special[1][0]
	assert complex(determinant) == pytest.approx(1, abs=1e-9)
	assert vqubit.is_special_unitary(special)


def test_state_and_density_matrix_kronecker_products() -> None:
	assert_vector_close(vqubit.state_vector_kronecker_product([1, 0], [0, 1]), [0, 1, 0, 0])
	density = vqubit.density_matrix_kronecker_product([[1, 0], [0, 0]], [[1, 0], [0, 0]])
	assert len(density) == 4
	assert len(density[0]) == 4
	assert_matrix_close(density, [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])


def test_targeted_left_multiply_on_each_axis() -> None:
	x = vqubit.unitary(vqubit.X)
	assert_matrix_close(vqubit.targeted_left_multiply(x, [[1, 2], [3, 4]], [0]), [[3, 4], [1, 2]])
	assert_matrix_close(vqubit.targeted_left_multiply(x, [[1, 2], [3, 4]], [1]), [[2, 1], [4, 3]])
	cnot = vqubit.unitary(vqubit.CNOT)
	assert_matrix_close(vqubit.targeted_left_multiply(cnot, [[0, 0], [1, 0]], [0, 1]), [[0, 0], [0, 1]])


def test_apply_matrix_to_non_adjacent_slices() -> None:
	result = vqubit.apply_matrix_to_slices([1, 2, 3, 4], [[0, 1], [1, 0]], [0, 3])
	assert_vector_close(result, [4, 2, 3, 1])


def test_sub_state_vector_factorizes_and_rejects_entanglement() -> None:
	# The extracted factor is defined up to a phase, so compare magnitudes.
	first = vqubit.sub_state_vector([1, 0, 0, 0], [0])
	assert [abs(z) for z in first] == pytest.approx([1, 0])
	second = vqubit.sub_state_vector([1, 0, 0, 0], [1])
	assert [abs(z) for z in second] == pytest.approx([1, 0])
	bell = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	assert vqubit.sub_state_vector(bell, [0], default="mixed") == "mixed"
	# A zero vector is not a state at all: the input is rejected, as in Cirq.
	with pytest.raises(ValueError, match="normalized"):
		vqubit.sub_state_vector([0, 0, 0, 0], [0], default=None)


def test_partial_trace_of_bell_state_as_mixture() -> None:
	bell = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	mixture = vqubit.partial_trace_of_state_vector_as_mixture(bell, [0])
	assert len(mixture) == 2
	assert [probability for probability, _ in mixture] == pytest.approx([0.5, 0.5])
	assert all(abs(sum(abs(z) ** 2 for z in vector) - 1) < 1e-9 for _, vector in mixture)
	pure = vqubit.partial_trace_of_state_vector_as_mixture([1, 0, 0, 0], [0])
	assert len(pure) == 1
	assert pure[0][0] == pytest.approx(1)


def test_diagonalize_real_symmetric_matrix_eigendecomposition() -> None:
	matrix = [[2, 1], [1, 2]]
	p = vqubit.diagonalize_real_symmetric_matrix(matrix)
	assert vqubit.is_orthogonal(p)
	assert vqubit.is_diagonal(vqubit.dot(_transpose(p), matrix, p))
	with pytest.raises(ValueError):
		vqubit.diagonalize_real_symmetric_matrix([[0, 1], [-1, 0]])


def test_diagonalize_commuting_pair() -> None:
	symmetric = [[2, 1, 0, 0], [1, 2, 0, 0], [0, 0, 3, 1], [0, 0, 1, 3]]
	diagonal = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]]
	p = vqubit.diagonalize_real_symmetric_and_sorted_diagonal_matrices(symmetric, diagonal)
	assert vqubit.is_orthogonal(p)
	assert vqubit.is_diagonal(vqubit.dot(_transpose(p), symmetric, p))
	assert vqubit.is_diagonal(vqubit.dot(_transpose(p), diagonal, p))
	with pytest.raises(ValueError):
		vqubit.diagonalize_real_symmetric_and_sorted_diagonal_matrices(
			[[2, 1], [1, 2]],
			[[1, 0], [0, -1]],
		)


def test_bidiagonalize_real_matrix_pair() -> None:
	# Commuting symmetric matrices have symmetric products, as the protocol requires.
	first = [[1, 2], [2, 1]]
	second = [[2, 1], [1, 2]]
	left, right = vqubit.bidiagonalize_real_matrix_pair_with_symmetric_products(first, second)
	assert vqubit.is_orthogonal(left)
	assert vqubit.is_orthogonal(right)
	assert vqubit.is_diagonal(vqubit.dot(left, first, right))
	assert vqubit.is_diagonal(vqubit.dot(left, second, right))
	with pytest.raises(ValueError, match="real"):
		vqubit.bidiagonalize_real_matrix_pair_with_symmetric_products([[1j, 0], [0, 1]], _identity(2))


def test_bidiagonalize_unitary_2x2() -> None:
	unitary = vqubit.unitary(vqubit.H)
	left, diagonal, right = vqubit.bidiagonalize_unitary_with_special_orthogonals(unitary)
	assert vqubit.is_special_orthogonal(left)
	assert vqubit.is_special_orthogonal(right)
	assert up_to_phase(vqubit.dot(left, unitary, right), diagonal, tol=1e-6)
	with pytest.raises(ValueError, match="unitary"):
		vqubit.bidiagonalize_unitary_with_special_orthogonals([[1, 0], [0, 2]])


def test_unitary_eig_reconstructs_the_matrix() -> None:
	unitary = vqubit.unitary(vqubit.S)
	values, vectors = vqubit.unitary_eig(unitary)
	assert sorted(abs(complex(v)) for v in values) == pytest.approx([1, 1])
	assert vqubit.is_unitary(vectors)
	diagonal = [[values[i] if i == j else 0j for j in range(2)] for i in range(2)]
	conjugate = [[vectors[j][i].conjugate() for j in range(2)] for i in range(2)]
	assert_matrix_close(vqubit.dot(vectors, diagonal, conjugate), unitary, tol=1e-6)
	with pytest.raises(ValueError, match="normal"):
		vqubit.unitary_eig([[1, 1], [0, 1]])


def test_map_eigenvalues_hermitian_and_normal() -> None:
	assert_matrix_close(vqubit.map_eigenvalues(vqubit.unitary(vqubit.X), lambda v: v * v), _identity(2))
	assert_matrix_close(
		vqubit.map_eigenvalues(vqubit.unitary(vqubit.S), lambda v: v),
		vqubit.unitary(vqubit.S),
		tol=1e-6,
	)


def test_so4_to_magic_su2s_returns_special_unitaries() -> None:
	theta = 0.3
	rotation = [[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]]
	matrix = [[rotation[i // 2][j // 2] * rotation[i % 2][j % 2] for j in range(4)] for i in range(4)]
	assert vqubit.is_special_orthogonal(matrix)
	a, b = vqubit.so4_to_magic_su2s(matrix)
	assert vqubit.is_special_unitary(a)
	assert vqubit.is_special_unitary(b)
	assert vqubit.is_unitary(vqubit.kron(a, b))
	with pytest.raises(ValueError):
		vqubit.so4_to_magic_su2s(vqubit.unitary(vqubit.SWAP))


def test_axis_angle_decomposition_of_known_unitaries() -> None:
	identity_decomposition = vqubit.axis_angle(_identity(2))
	assert identity_decomposition.angle == pytest.approx(0)
	assert identity_decomposition.global_phase == pytest.approx(1)
	assert tuple(identity_decomposition.axis) == pytest.approx((1, 0, 0))
	x_decomposition = vqubit.axis_angle(vqubit.unitary(vqubit.X))
	assert up_to_phase(x_decomposition.matrix(), vqubit.unitary(vqubit.X), tol=1e-6)
	assert_matrix_close(x_decomposition._unitary_(), x_decomposition.matrix())
	with pytest.raises(ValueError):
		vqubit.axis_angle([[1, 0], [0, 2]])
	with pytest.raises(ValueError, match="normalized"):
		vqubit.AxisAngleDecomposition(angle=1, axis=(1, 1, 0), global_phase=1)


def test_kak_decomposition_class_matrix() -> None:
	identity_decomposition = vqubit.KakDecomposition(interaction_coefficients=(0, 0, 0))
	assert_matrix_close(identity_decomposition.matrix(), _identity(4))
	decomposition = vqubit.KakDecomposition(interaction_coefficients=(math.pi / 4, 0, 0))
	s = 1 / math.sqrt(2)
	expected = [
		[s, 0, 0, 1j * s],
		[0, s, 1j * s, 0],
		[0, 1j * s, s, 0],
		[1j * s, 0, 0, s],
	]
	assert_matrix_close(decomposition.matrix(), expected)
	assert vqubit.is_unitary(decomposition.matrix())
	assert "KakDecomposition" in repr(decomposition)


def test_kak_canonicalize_vector_orders_coefficients() -> None:
	decomposition = vqubit.kak_canonicalize_vector(0.1, 0.2, 0.3)
	x, y, z = decomposition.interaction_coefficients
	assert (x, y, z) == pytest.approx((0.3, 0.2, 0.1))
	assert abs(z) <= y + 1e-9 <= x + 1e-9 <= math.pi / 4 + 1e-9
	assert vqubit.is_unitary(decomposition.matrix())
	assert vqubit.kak_canonicalize_vector(0, 0, 0).interaction_coefficients == pytest.approx((0, 0, 0))


def test_kak_vector_single_batch_and_errors() -> None:
	vector = vqubit.kak_vector(vqubit.unitary(vqubit.CNOT))
	assert list(vector) == pytest.approx([math.pi / 4, 0, 0], abs=1e-6)
	batch = vqubit.kak_vector([vqubit.unitary(vqubit.CNOT), _identity(4)])
	assert len(batch) == 2
	assert len(batch[0]) == 3
	assert list(batch[0]) == pytest.approx([math.pi / 4, 0, 0], abs=1e-6)
	assert list(batch[1]) == pytest.approx([0, 0, 0], abs=1e-6)
	assert len(vqubit.kak_vector([])) == 0
	with pytest.raises(ValueError):
		vqubit.kak_vector(_identity(2))
	with pytest.raises(ValueError):
		vqubit.kak_decomposition([[1, 0], [0, 1]])


def test_num_cnots_required_for_known_gates() -> None:
	assert vqubit.num_cnots_required(_identity(4)) == 0
	assert vqubit.num_cnots_required(vqubit.unitary(vqubit.CNOT)) == 1
	assert vqubit.num_cnots_required(vqubit.unitary(vqubit.ISWAP)) == 2
	assert vqubit.num_cnots_required(vqubit.unitary(vqubit.SWAP)) == 3
	with pytest.raises(ValueError, match="shape"):
		vqubit.num_cnots_required(_identity(2))


def test_match_global_phase_with_generic_complex_entry() -> None:
	a = [[1 + 1j, 0], [0, 1 - 1j]]
	# Multiply by a unit phase, so magnitudes are preserved.
	b = [[1j * z for z in row] for row in a]
	matched_a, matched_b = vqubit.match_global_phase(a, b)
	assert_matrix_close(matched_a, matched_b)
	assert vqubit.allclose_up_to_global_phase(a, b)


def test_targeted_left_multiply_on_a_vector() -> None:
	assert_vector_close(vqubit.targeted_left_multiply(vqubit.unitary(vqubit.X), [1, 2], [0]), [2, 1])


def test_diagonalize_pair_rejects_complex_symmetric_matrix() -> None:
	with pytest.raises(ValueError, match="symmetric_matrix"):
		vqubit.diagonalize_real_symmetric_and_sorted_diagonal_matrices(
			[[1, 1j], [1j, 2]],
			[[1, 0], [0, 2]],
		)


def test_axis_angle_canonicalization_and_repr() -> None:
	decomposition = vqubit.AxisAngleDecomposition(angle=3 * math.pi, axis=(1, 0, 0), global_phase=1)
	canonical = decomposition.canonicalize()
	assert -math.pi <= canonical.angle <= math.pi
	assert_matrix_close(canonical.matrix(), decomposition.matrix(), tol=1e-9)
	assert "AxisAngleDecomposition" in repr(canonical)


def test_slice_for_qubits_equal_to_requires_one_endian_value() -> None:
	with pytest.raises(ValueError, match="exactly one"):
		vqubit.slice_for_qubits_equal_to([0], 1, big_endian_qureg_value=1)


def test_slice_for_qubits_equal_to_little_and_big_endian() -> None:
	assert vqubit.slice_for_qubits_equal_to([0, 2], 1, num_qubits=3) == (1, slice(None), 0)
	assert vqubit.slice_for_qubits_equal_to([0, 2], big_endian_qureg_value=1, num_qubits=3) == (
		0,
		slice(None),
		1,
	)
	assert vqubit.slice_for_qubits_equal_to([0]) == (0, Ellipsis)
	assert vqubit.slice_for_qubits_equal_to([0], qid_shape=(2, 3)) == (0, slice(None))
	with pytest.raises(ValueError, match="len"):
		vqubit.slice_for_qubits_equal_to([0], num_qubits=2, qid_shape=(2, 2, 2))
