"""Regression tests for bugs found while auditing the numerical core.

Each test names the behaviour that was wrong, so a future refactor that brings
it back fails with a message that explains what is at stake. The cases come
from the coverage and parity audit: wrong eigenvalues for unitary matrices, a
reduced state that picked the wrong factor, a "Kronecker product" that was not
one, and a handful of protocol details.
"""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close

np = pytest.importorskip("numpy", exc_type=ImportError)

# --- the eigensolver ---------------------------------------------------------


@pytest.mark.parametrize(
	("matrix", "expected"),
	[
		# A rotation has eigenvalues i and -i: no real shift finds them, and the
		# old solver returned [0, 0].
		([[0, -1], [1, 0]], {1j, -1j}),
		# A permutation matrix has three eigenvalues of modulus one; the QR
		# iteration used to stall on it.
		([[0, 0, 1], [1, 0, 0], [0, 1, 0]], {1, complex(-0.5, math.sqrt(3) / 2), complex(-0.5, -math.sqrt(3) / 2)}),
		# The Hadamard and the identity are their own eigenproblems.
		([[1, 1], [1, -1]], {math.sqrt(2), -math.sqrt(2)}),
		([[1, 0], [0, 1]], {1}),
	],
)
def test_unitary_eig_finds_the_eigenvalues_of_a_normal_matrix(matrix: list, expected: set) -> None:
	values, vectors = vqubit.unitary_eig(matrix)
	# Match as a multiset of eigenvalues, up to numerical noise.
	remaining = list(values)
	for wanted in expected:
		closest = min(remaining, key=lambda value: abs(value - wanted))
		assert abs(closest - wanted) < 1e-7, (matrix, values)
		remaining.remove(closest)
	# Every eigenpair holds, and the basis is unitary (the columns are the
	# eigenvectors, as in Cirq).
	matrix_array = np.asarray(matrix, dtype=complex)
	basis = np.asarray(vectors)
	for index, value in enumerate(values):
		assert np.allclose(matrix_array @ basis[:, index], value * basis[:, index], atol=1e-7)
	assert np.allclose(basis.conj().T @ basis, np.eye(len(values)), atol=1e-7)


def test_map_eigenvalues_of_a_rotation() -> None:
	# squaring the eigenvalues of the 90-degree rotation gives -I.
	assert_matrix_close(np.asarray(vqubit.map_eigenvalues([[0, -1], [1, 0]], lambda v: v * v)), -np.eye(2))
	# ... and mapping to the identity gives the identity whatever the input.
	assert_matrix_close(np.asarray(vqubit.map_eigenvalues([[0, -1], [1, 0]], lambda v: 1 + 0j)), np.eye(2))


def test_unitary_eig_rejects_a_non_normal_matrix() -> None:
	with pytest.raises(ValueError, match="normal"):
		vqubit.unitary_eig([[1, 2], [0, 1]])


# --- the reduced state -------------------------------------------------------


def test_sub_state_vector_returns_the_kept_factor() -> None:
	# |10> on two qubits: qubit 0 is |1>, qubit 1 is |0>.
	assert [abs(z) for z in vqubit.sub_state_vector([0, 0, 1, 0], [0])] == pytest.approx([0, 1])
	assert [abs(z) for z in vqubit.sub_state_vector([0, 1, 0, 0], [0])] == pytest.approx([1, 0])
	assert [abs(z) for z in vqubit.sub_state_vector([0, 1, 0, 0], [1])] == pytest.approx([0, 1])
	# Keeping every qubit gives the state back, and keeping none the scalar 1.
	assert [abs(z) for z in vqubit.sub_state_vector([0, 0, 1, 0], [0, 1])] == pytest.approx([0, 0, 1, 0])
	assert list(vqubit.sub_state_vector([1, 0, 0, 0], [])) == [1]
	# A nested (row shaped) state vector is accepted too.
	assert [abs(z) for z in vqubit.sub_state_vector([[0, 0], [1, 0]], [0])] == pytest.approx([0, 1])


def test_sub_state_vector_rejects_what_it_cannot_factor() -> None:
	bell = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	assert vqubit.sub_state_vector(bell, [0], default="mixed") == "mixed"
	with pytest.raises(ValueError, match="normalized"):
		vqubit.sub_state_vector([0, 0, 0, 0], [0])
	with pytest.raises(ValueError, match="does not represent a state"):
		vqubit.sub_state_vector([1, 0, 0], [0])
	with pytest.raises(ValueError, match="unique"):
		vqubit.sub_state_vector([1, 0, 0, 0], [0, 0])
	with pytest.raises(ValueError, match="invalid subset"):
		vqubit.sub_state_vector([1, 0, 0, 0], [2])


# --- the linear-algebra helpers ---------------------------------------------


def test_density_matrix_kronecker_product_is_a_kronecker_product() -> None:
	result = np.asarray(vqubit.density_matrix_kronecker_product([[1, 0], [0, 0]], [[0, 0], [0, 1]]))
	# |0><0| ⊗ |1><1| is the basis matrix |01><01|, with trace one.
	assert_matrix_close(result, np.diag([0, 1, 0, 0]))
	assert np.trace(result) == pytest.approx(1)
	# Non-square factors keep the block structure of np.kron.
	tall = np.asarray(vqubit.density_matrix_kronecker_product(np.ones((2, 3)).tolist(), np.ones((2, 2)).tolist()))
	assert tall.shape == (4, 6)


def test_targeted_conjugate_about_works_with_numpy_installed() -> None:
	# A · X · A^† on a single qubit, with A the Hadamard: this is the X basis
	# change, so the conjugated X is Z.
	hadamard = np.asarray(vqubit.unitary(vqubit.H))
	conjugated = np.asarray(vqubit.targeted_conjugate_about(hadamard, vqubit.unitary(vqubit.X), [0]))
	assert_matrix_close(conjugated, vqubit.unitary(vqubit.Z))


def test_kron_scales_with_a_scalar() -> None:
	x = vqubit.unitary(vqubit.X)
	assert_matrix_close(np.asarray(vqubit.kron(2, x)), 2 * np.asarray(x))
	assert_matrix_close(np.asarray(vqubit.kron(x, 2)), 2 * np.asarray(x))
	assert_matrix_close(np.asarray(vqubit.kron(2, x, 3)), 6 * np.asarray(x))
	assert vqubit.kron(2) == 2


def test_diagonal_and_empty_matrices() -> None:
	# A rectangular matrix is never diagonal, even if it is "diagonal" as written.
	assert not vqubit.is_diagonal([[1, 0, 0], [0, 1, 0]])
	assert vqubit.is_diagonal([[1, 0], [0, 2]])
	# The 0x0 matrix is square and vacuously fine.
	assert vqubit.is_unitary([])
	assert vqubit.is_hermitian([])


def test_so4_to_magic_su2s_requires_a_four_by_four_matrix() -> None:
	with pytest.raises(ValueError, match="4x4"):
		vqubit.so4_to_magic_su2s([[1, 0], [0, 1]])
	first, second = vqubit.so4_to_magic_su2s([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
	assert_matrix_close(np.asarray(first), np.eye(2))
	assert_matrix_close(np.asarray(second), np.eye(2))


# --- protocols ---------------------------------------------------------------


def test_parameter_names_looks_inside_expressions_and_circuits() -> None:
	t = vqubit.Parameter("t")
	q0 = vqubit.LineQubit(0)
	assert vqubit.parameter_names(vqubit.X**t) == {"t"}  # ty: ignore[unsupported-operator]
	assert vqubit.parameter_names(vqubit.X ** (2 * t + 1)) == {"t"}  # ty: ignore[unsupported-operator]
	assert vqubit.parameter_names(vqubit.rz(t)) == {"t"}
	assert vqubit.parameter_names(vqubit.Circuit((vqubit.X ** (t / 2)).on(q0))) == {"t"}  # ty: ignore[unsupported-operator]
	assert vqubit.parameter_names(vqubit.X) == set()


def test_qid_shape_follows_the_declared_dimensions() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.qid_shape(vqubit.X) == (2,)
	assert vqubit.qid_shape(vqubit.CNOT) == (2, 2)
	assert vqubit.qid_shape(q0) == (2,)
	assert vqubit.qid_shape(vqubit.X(q0)) == (2,)
	assert vqubit.qid_shape(vqubit.Circuit(vqubit.CNOT(q0, q1))) == (2, 2)
	# A qutrit channel keeps its dimension instead of being reported as a qubit.
	assert vqubit.qid_shape(vqubit.ResetChannel(3)) == (3,)


def test_unitary_of_a_partly_idle_decomposition_spans_every_qubit() -> None:
	# The decomposition of a uniform superposition over [0, 2) only touches one
	# qubit, but the gate acts on two.
	assert np.asarray(vqubit.unitary(vqubit.UniformSuperpositionGate(2, 2))).shape == (4, 4)
	assert np.asarray(vqubit.unitary(vqubit.UniformSuperpositionGate(1, 2))).shape == (4, 4)


def test_mixed_kinds_of_qid_sort_deterministically() -> None:
	ordered = sorted([vqubit.NamedQubit("a"), vqubit.LineQubit(1), vqubit.GridQubit(0, 0), vqubit.LineQubit(0)])
	assert ordered == [vqubit.LineQubit(0), vqubit.LineQubit(1), vqubit.GridQubit(0, 0), vqubit.NamedQubit("a")]


# --- Pauli strings -----------------------------------------------------------


def test_pauli_string_is_a_mapping_of_its_qubits() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	string = vqubit.PauliString({q0: vqubit.Pauli.X, q1: vqubit.Pauli.Y})
	# Iterating yields the qubits (as in Cirq), and it terminates.
	assert list(string) == [q0, q1]
	assert len(string) == 2
	assert list(string.values()) == [vqubit.Pauli.X, vqubit.Pauli.Y]
	assert string.get(vqubit.LineQubit(9)) is vqubit.Pauli.I
	assert string.get(vqubit.LineQubit(9), vqubit.Pauli.Z) is vqubit.Pauli.Z


def test_pauli_string_powers_use_the_parity() -> None:
	q0 = vqubit.LineQubit(0)
	string = vqubit.PauliString({q0: vqubit.Pauli.X})
	assert string**-1 == string
	assert string**3 == string
	assert (string**2).qubits == ()


def test_pauli_sum_has_value_equality() -> None:
	q0 = vqubit.LineQubit(0)
	first = vqubit.PauliSum([(vqubit.PauliString({q0: vqubit.Pauli.X}), 0.5)])
	second = vqubit.PauliSum([(vqubit.PauliString({q0: vqubit.Pauli.X}), 0.5)])
	assert first == second
	assert first != vqubit.PauliSum([(vqubit.PauliString({q0: vqubit.Pauli.X}), 0.25)])
	assert hash(first) == hash(second)


# --- gates -------------------------------------------------------------------


def test_controlled_gate_resolves_parameters_of_its_sub_gate() -> None:
	t = vqubit.Parameter("t")
	gate = vqubit.ControlledGate(vqubit.ZPowGate(exponent=t), num_controls=1)
	assert vqubit.is_parameterized(gate)
	resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({"t": 1}))
	assert not vqubit.is_parameterized(resolved)
	assert resolved == vqubit.ControlledGate(vqubit.Z, num_controls=1)


def test_parameter_expressions_have_an_evaluable_repr() -> None:
	t = vqubit.Parameter("t")
	expression = t / math.pi
	assert eval(repr(expression), {"vqubit": vqubit, "math": math}) == expression
	assert eval(repr(-t), {"vqubit": vqubit}) == -t
