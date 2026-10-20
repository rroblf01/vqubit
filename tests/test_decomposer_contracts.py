"""Tests for the closing contracts of the 0.5 audit (ROADMAP section 9).

Each test pins a behaviour the audit found missing: the decomposers reject a
matrix they cannot handle instead of panicking in the native core, the
diagonal-plus-CZ decomposition returns the real diagonal, the QASM importer
honours the bit a measurement goes to, ``ApplyUnitaryArgs`` really transposes
its axes, and the bidiagonalization survives a rank-deficient pair.
"""

from __future__ import annotations

import numpy as np
import pytest

import vqubit
from tests.helpers import assert_matrix_close, up_to_phase

Q0, Q1 = vqubit.LineQubit.range(2)
IDENTITY_4 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]


@pytest.mark.parametrize(
	("call", "match"),
	[
		(lambda: vqubit.kron_factor_4x4_to_2x2s([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), "4x4"),
		(lambda: vqubit.two_qubit_matrix_to_cz_isometry(Q0, Q1, [[1, 0], [0, 1]]), "isometry"),
		(lambda: vqubit.single_qubit_matrix_to_gates([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), "2x2"),
		(lambda: vqubit.single_qubit_matrix_to_phxz([[1, 1], [0, 1]]), "unitary"),
		(lambda: vqubit.deconstruct_single_qubit_matrix_into_angles([[1, 1], [0, 1]]), "unitary"),
		(lambda: vqubit.single_qubit_op_to_framed_phase_form([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), "2x2"),
		(lambda: vqubit.prepare_two_qubit_state_using_cz([1, 0, 0], [Q0, Q1]), "4 amplitudes"),
		(lambda: vqubit.prepare_two_qubit_state_using_iswap([2, 0, 0, 0], [Q0, Q1]), "normalized"),
	],
)
def test_decomposers_raise_value_error_instead_of_panicking(call: object, match: str) -> None:
	# The audit found these panicking inside the V core (RuntimeError) or
	# silently accepting nonsense.
	with pytest.raises(ValueError, match=match):
		call()  # ty: ignore[call-non-callable]


def test_the_diagonal_part_really_is_diagonal() -> None:
	def rebuilt(matrix: object) -> list[list[complex]]:
		diagonal, operations = vqubit.two_qubit_matrix_to_diagonal_and_cz_operations(Q0, Q1, matrix)  # ty: ignore[invalid-argument-type]
		ops_matrix = vqubit.unitary(vqubit.Circuit(operations)) if operations else IDENTITY_4
		return vqubit.linalg.dot(ops_matrix, vqubit.unitary(diagonal))

	# A three-CNOT unitary gets a non-trivial diagonal; a two-CNOT one an
	# identity diagonal; a diagonal one no operations at all.
	for gate in (vqubit.SWAP, vqubit.SWAP**0.5, vqubit.CNOT, vqubit.ISWAP):
		assert up_to_phase(rebuilt(vqubit.unitary(gate)), vqubit.unitary(gate)), gate

	diagonal, _ = vqubit.two_qubit_matrix_to_diagonal_and_cz_operations(Q0, Q1, vqubit.unitary(vqubit.SWAP))
	assert any(abs(angle) > 1e-6 for angle in diagonal.diag_angles_radians)

	target = vqubit.DiagonalGate([0.1, 0.2, 0.3, 0.4])
	diagonal, operations = vqubit.two_qubit_matrix_to_diagonal_and_cz_operations(Q0, Q1, vqubit.unitary(target))
	assert operations == []
	assert up_to_phase(vqubit.unitary(diagonal), vqubit.unitary(target))


def test_qasm_measurement_bits_keep_their_position() -> None:
	# Measuring q[0] into c[1] and q[1] into c[0] swaps the bits: the measurement
	# gate must name the qubits in bit order.
	program = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
measure q[0] -> c[1];
measure q[1] -> c[0];
"""
	circuit = vqubit.circuit_from_qasm(program)
	operation = next(iter(circuit.all_operations()))
	assert vqubit.measurement_key_name(operation) == "c"
	assert operation.qubits == (vqubit.LineQubit(1), vqubit.LineQubit(0))
	# The classical-assignment form follows the same rule.
	assigned = vqubit.circuit_from_qasm(
		"""OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
"""
	)
	assert next(iter(assigned.all_operations())).qubits == (Q0, Q1)


def test_apply_args_transpose_the_target_axes_to_the_front() -> None:
	# A three-qubit tensor, acting on axis 1: the transposed tensor swaps axes
	# 1 and 0 and renumbers the axes to 0.
	tensor = np.arange(8, dtype=complex).reshape(2, 2, 2)
	moved = vqubit.ApplyUnitaryArgs(tensor, tensor.copy(), [1]).with_axes_transposed_to_start()
	assert moved.axes == [0]
	assert np.array_equal(moved.target_tensor, np.transpose(tensor, (1, 0, 2)))
	assert np.array_equal(moved.available_buffer, np.transpose(tensor, (1, 0, 2)))

	# The nested-list form has to agree with numpy.
	nested = [[[complex(i * 4 + j * 2 + k) for k in range(2)] for j in range(2)] for i in range(2)]
	plain = vqubit.ApplyUnitaryArgs(nested, [[[0j] * 2] * 2] * 2, [2]).with_axes_transposed_to_start()
	assert plain.axes == [0]
	for a in range(2):
		for b in range(2):
			for c in range(2):
				assert plain.target_tensor[a][b][c] == nested[b][c][a]

	# No axes is a no-op, and an out-of-range axis is an error.
	assert vqubit.ApplyUnitaryArgs(nested).with_axes_transposed_to_start().axes == []
	with pytest.raises(ValueError, match="not axes"):
		vqubit.ApplyUnitaryArgs(nested, nested, [3]).with_axes_transposed_to_start()


@pytest.mark.parametrize(
	("first", "second"),
	[
		# Rank one: the second singular value is zero and the left basis has to
		# be completed, which used to break the orthogonality.
		([[1, 1], [1, 1]], [[2, 2], [2, 2]]),
		([[1, 1], [1, 1]], [[2, 1], [1, 2]]),
		([[0, 0], [0, 0]], [[2, 1], [1, 2]]),
		([[1, 0], [0, 1]], [[2, 1], [1, 2]]),
	],
)
def test_bidiagonalization_survives_a_rank_deficient_pair(first: list, second: list) -> None:
	left, right = vqubit.bidiagonalize_real_matrix_pair_with_symmetric_products(first, second)
	left, right = np.asarray(left), np.asarray(right)
	# Both are orthogonal and both products come out diagonal.
	assert_matrix_close(left @ left.T, np.eye(2))
	assert_matrix_close(right @ right.T, np.eye(2))
	for matrix in (first, second):
		product = np.asarray(left @ np.asarray(matrix, dtype=float) @ right).real
		assert_matrix_close(product, np.diag(np.diag(product)))


def test_bidiagonalization_rejects_what_it_cannot_diagonalize() -> None:
	with pytest.raises(ValueError, match="real"):
		vqubit.bidiagonalize_real_matrix_pair_with_symmetric_products([[1, 0], [0, 1j]], [[1, 0], [0, 1]])
	# mat1 @ mat2.T must be symmetric for the construction to exist.
	with pytest.raises(ValueError, match=r"mat1 @ mat2\.T"):
		vqubit.bidiagonalize_real_matrix_pair_with_symmetric_products([[1, 2], [2, 3]], [[2, 1], [1, 2]])


def test_isometry_completion_is_unitary() -> None:
	# Two orthonormal columns of a unitary complete to a unitary that keeps them
	# intact.
	full = np.asarray(vqubit.unitary(vqubit.ISWAP**0.3))
	isometry = full[:, :2]
	operations = vqubit.two_qubit_matrix_to_cz_isometry(Q0, Q1, isometry.tolist())
	matrix = np.asarray(vqubit.unitary(vqubit.Circuit(operations)))
	assert np.isclose(abs(np.linalg.det(matrix)), 1, atol=1e-6)
	for column in range(2):
		vector, reference = matrix[:, column], isometry[:, column]
		phase = vector[np.argmax(np.abs(reference))] / reference[np.argmax(np.abs(reference))]
		assert np.allclose(vector, phase * reference)
