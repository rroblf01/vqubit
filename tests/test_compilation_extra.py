"""Behavioural tests for the analytical decompositions in ``_compilation.py``.

Each decomposition is exercised through the public API and checked by building a
circuit from the returned operations and comparing its unitary with the matrix
that was decomposed.  The point is the contract, not the exact gate list, so most
comparisons are up to a global phase.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence
from typing import Any

import pytest

import vqubit
from tests.helpers import random_unitary, up_to_phase

_PAULI_GATE = {vqubit.Pauli.X: vqubit.X, vqubit.Pauli.Y: vqubit.Y, vqubit.Pauli.Z: vqubit.Z}


def _operations_circuit(operations: Any) -> vqubit.Circuit:
	return operations if isinstance(operations, vqubit.Circuit) else vqubit.Circuit(operations)


def _assert_implements(operations: Any, matrix: Sequence[Sequence[complex]], atol: float = 1e-5) -> None:
	circuit = _operations_circuit(operations)
	assert vqubit.linalg.allclose_up_to_global_phase(vqubit.unitary(circuit), matrix, atol=atol)


def _block_with_ancillas_zero(circuit: vqubit.Circuit, active_qubits: int) -> Sequence[Sequence[complex]]:
	"""Extract the block where the trailing ancilla qubits are ``|0>``.

	The decompositions put fresh ancillas on higher-numbered ``LineQubit``s, so
	they are the last bits in the circuit's sorted qubit order.  Inputs and
	outputs are restricted to ancilla ``0`` to compare with the operation on the
	active qubits alone.
	"""
	qubits = circuit.all_qubits()
	extra = len(qubits) - active_qubits
	matrix = vqubit.unitary(circuit)
	if extra == 0:
		return matrix
	sub = 1 << extra
	dim = 1 << active_qubits
	return [[matrix[a * sub][b * sub] for b in range(dim)] for a in range(dim)]


def _assert_vectors_close_up_to_phase(
	actual: Sequence[complex], expected: Sequence[complex], atol: float = 1e-6
) -> None:
	phase = None
	for a, e in zip(actual, expected, strict=True):
		if abs(e) > 1e-9:
			phase = a / e
			break
	assert phase is not None
	for a, e in zip(actual, expected, strict=True):
		assert abs(a - phase * e) <= atol


def _assert_states_close(actual: Sequence[complex], expected: Sequence[complex], atol: float = 1e-7) -> None:
	_assert_vectors_close_up_to_phase(actual, expected, atol=atol)


def test_single_qubit_matrix_to_pauli_rotations_is_a_zyz_sequence() -> None:
	for gate in (vqubit.I, vqubit.X, vqubit.Y, vqubit.Z, vqubit.T):
		rotations = vqubit.single_qubit_matrix_to_pauli_rotations(vqubit.unitary(gate))
		assert [pauli for pauli, _ in rotations] == [vqubit.Pauli.Z, vqubit.Pauli.Y, vqubit.Pauli.Z]
		assert len(rotations) == 3


def test_single_qubit_matrix_to_pauli_rotations_reconstructs_the_unitary() -> None:
	q0 = vqubit.LineQubit(0)
	rng = random.Random(0)
	matrices = [vqubit.unitary(vqubit.H), vqubit.unitary(vqubit.T), *[random_unitary(2, rng) for _ in range(3)]]
	for matrix in matrices:
		rotations = vqubit.single_qubit_matrix_to_pauli_rotations(matrix)
		operations = [(_PAULI_GATE[pauli] ** turns).on(q0) for pauli, turns in rotations]
		assert up_to_phase(vqubit.unitary(vqubit.Circuit(operations)), matrix, tol=1e-6)


def test_single_qubit_op_to_framed_phase_form() -> None:
	# The contract is M = U^-1 · diag(1, r) · U · diag(g, g), as in Cirq: U
	# rotates the axis onto Z, r is the relative phase and g the global one.
	def conjugate_transpose(matrix: Any) -> Any:
		return [[matrix[j][i].conjugate() for j in range(len(matrix))] for i in range(len(matrix))]

	for gate in (
		vqubit.H,
		vqubit.X,
		vqubit.Y,
		vqubit.Z,
		vqubit.S,
		vqubit.T,
		vqubit.rx(0.3),
		vqubit.ry(1.1),
		vqubit.rz(0.9),
	):
		matrix = vqubit.unitary(gate)
		u, r, g = vqubit.single_qubit_op_to_framed_phase_form(matrix)
		assert abs(abs(r) - 1) < 1e-8
		assert abs(abs(g) - 1) < 1e-8
		rebuilt = vqubit.linalg.dot(conjugate_transpose(u), [[1, 0], [0, r]])
		rebuilt = vqubit.linalg.dot(rebuilt, u)
		rebuilt = vqubit.linalg.dot(rebuilt, [[g, 0], [0, g]])
		assert up_to_phase(rebuilt, matrix, tol=1e-6), gate


def test_single_qubit_op_to_framed_phase_form_of_a_z_rotation_is_diagonal() -> None:
	# A Z rotation already rotates about the Z axis, so U is a permutation and r
	# is the phase between the two eigenvalues.
	matrix = vqubit.unitary(vqubit.rz(0.7))
	_, r, _ = vqubit.single_qubit_op_to_framed_phase_form(matrix)
	assert abs(abs(r) - 1) < 1e-8
	assert abs(r - 1) > 1e-3  # the two eigenvalues differ


def test_cphase_intervals_handle_angles_past_a_quarter_turn() -> None:
	# theta = 2 rad puts |arcsin(sin theta)| on the sine's descending side, and
	# the mirror interval comes first, both of which are branches of the
	# interval computation.
	intervals = vqubit.compute_cphase_exponents_for_fsim_decomposition(vqubit.FSimGate(2.0, 0.0))
	assert all(0 <= start < stop <= 2 for start, stop in intervals)


def test_cphase_decomposition_rejects_a_parametrized_fsim() -> None:
	with pytest.raises(ValueError, match="parametrized gate"):
		vqubit.decompose_cphase_into_two_fsim(vqubit.CZ, fsim_gate=vqubit.FSimGate(vqubit.Parameter("t"), 0.0))


def test_cphase_decomposition_lists_the_exponents_when_none_fit() -> None:
	# An fSim with a large cphase cannot build this CZ: the error names the
	# intervals that would work.
	with pytest.raises(ValueError, match="Valid intervals"):
		vqubit.decompose_cphase_into_two_fsim(vqubit.CZ**0.25, fsim_gate=vqubit.FSimGate(0.1, 0.05))


def test_four_fsim_synthesis_rejects_bad_arguments() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError, match="parameterized"):
		vqubit.decompose_two_qubit_interaction_into_four_fsim_gates(
			vqubit.ISWAP, fsim_gate=vqubit.FSimGate(vqubit.Parameter("t"), 0.0)
		)
	with pytest.raises(ValueError, match="pair of qubits"):
		vqubit.decompose_two_qubit_interaction_into_four_fsim_gates(
			vqubit.ISWAP, fsim_gate=vqubit.FSimGate(math.pi / 2, 0.0), qubits=(q0,)
		)


def test_two_qubit_matrix_to_diagonal_and_cz_extracts_the_diagonal() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	identity = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

	def rebuilt(matrix: list[list[complex]]) -> list[list[complex]]:
		diagonal, operations = vqubit.two_qubit_matrix_to_diagonal_and_cz_operations(q0, q1, matrix)
		ops_matrix = vqubit.unitary(vqubit.Circuit(operations)) if operations else identity
		return vqubit.linalg.dot(ops_matrix, vqubit.unitary(diagonal))

	for gate in (vqubit.CNOT, vqubit.SWAP, vqubit.ISWAP, vqubit.SWAP**0.5):
		assert up_to_phase(rebuilt(vqubit.unitary(gate)), vqubit.unitary(gate)), gate

	# A three-CNOT unitary comes back with a non-trivial diagonal part, which is
	# the whole point of the decomposition.
	diagonal, _operations = vqubit.two_qubit_matrix_to_diagonal_and_cz_operations(q0, q1, vqubit.unitary(vqubit.SWAP))
	assert any(abs(angle) > 1e-6 for angle in diagonal.diag_angles_radians)

	# A diagonal unitary needs no operations at all: the diagonal is the gate.
	diagonal_gate = vqubit.DiagonalGate([0.1, 0.2, 0.3, 0.4])
	diagonal, operations = vqubit.two_qubit_matrix_to_diagonal_and_cz_operations(q0, q1, vqubit.unitary(diagonal_gate))
	assert operations == []
	assert up_to_phase(vqubit.unitary(diagonal), vqubit.unitary(diagonal_gate))


def test_decompose_multi_controlled_x_with_no_controls_is_x() -> None:
	target = vqubit.LineQubit(0)
	operations = vqubit.decompose_multi_controlled_x([], target)
	assert operations == [vqubit.X(target)]


@pytest.mark.parametrize("num_controls", [1, 2])
def test_decompose_multi_controlled_x_one_and_two_controls(num_controls: int) -> None:
	qubits = vqubit.LineQubit.range(num_controls + 1)
	controls, target = qubits[:num_controls], qubits[num_controls]
	operations = vqubit.decompose_multi_controlled_x(controls, target)
	reference = vqubit.unitary(vqubit.ControlledGate(vqubit.X, num_controls=num_controls))
	_assert_implements(operations, reference, atol=1e-6)


def test_decompose_multi_controlled_x_three_controls_uses_fresh_ancillas() -> None:
	qubits = vqubit.LineQubit.range(4)
	controls, target = qubits[:3], qubits[3]
	operations = vqubit.decompose_multi_controlled_x(controls, target)
	circuit = vqubit.Circuit(operations)
	assert set(circuit.all_qubits()) == set(vqubit.LineQubit.range(6))
	reference = vqubit.unitary(vqubit.ControlledGate(vqubit.X, num_controls=3))
	assert vqubit.linalg.allclose_up_to_global_phase(_block_with_ancillas_zero(circuit, 4), reference, atol=1e-6)


@pytest.mark.parametrize("num_controls", [1, 2, 3])
def test_decompose_multi_controlled_x_with_false_control_values(num_controls: int) -> None:
	qubits = vqubit.LineQubit.range(num_controls + 1)
	controls, target = qubits[:num_controls], qubits[num_controls]
	control_values = [False] * num_controls
	operations = vqubit.decompose_multi_controlled_x(controls, target, control_values=control_values)
	circuit = vqubit.Circuit(operations)
	block = _block_with_ancillas_zero(circuit, num_controls + 1)
	reference = vqubit.unitary(
		vqubit.ControlledGate(vqubit.X, num_controls=num_controls, control_values=control_values)
	)
	assert vqubit.linalg.allclose_up_to_global_phase(block, reference, atol=1e-6)


def test_decompose_multi_controlled_x_validates_control_values() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	with pytest.raises(ValueError):
		vqubit.decompose_multi_controlled_x([q0, q1], q2, control_values=[True])


def test_decompose_multi_controlled_rotation_without_controls_is_a_matrix_gate() -> None:
	q0 = vqubit.LineQubit(0)
	matrix = vqubit.unitary(vqubit.H)
	operations = vqubit.decompose_multi_controlled_rotation(matrix, [], q0)
	assert len(operations) == 1
	assert isinstance(operations[0].gate, vqubit.MatrixGate)
	assert up_to_phase(vqubit.unitary(vqubit.Circuit(operations)), matrix, tol=1e-6)


@pytest.mark.parametrize("control_value", [True, False])
def test_decompose_multi_controlled_rotation_single_control(control_value: bool) -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(1)
	for _ in range(3):
		matrix = random_unitary(2, rng)
		operations = vqubit.decompose_multi_controlled_rotation(matrix, [q0], q1, control_values=[control_value])
		reference = vqubit.unitary(
			vqubit.ControlledGate(vqubit.MatrixGate(matrix), num_controls=1, control_values=[control_value])
		)
		_assert_implements(operations, reference, atol=1e-6)


def test_decompose_multi_controlled_rotation_identity_has_no_phase_correction() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operations = vqubit.decompose_multi_controlled_rotation(vqubit.unitary(vqubit.I), [q0], q1)
	assert not any(isinstance(operation.gate, vqubit.ZPowGate) for operation in operations)
	reference = vqubit.unitary(vqubit.ControlledGate(vqubit.MatrixGate(vqubit.unitary(vqubit.I)), num_controls=1))
	_assert_implements(operations, reference, atol=1e-6)


@pytest.mark.parametrize("num_controls", [2, 3])
def test_decompose_multi_controlled_rotation_with_multiple_controls(num_controls: int) -> None:
	qubits = vqubit.LineQubit.range(num_controls + 1)
	controls, target = qubits[:num_controls], qubits[num_controls]
	rng = random.Random(2)
	for control_values in ([True] * num_controls, [False] * num_controls):
		matrix = random_unitary(2, rng)
		operations = vqubit.decompose_multi_controlled_rotation(matrix, controls, target, control_values=control_values)
		circuit = vqubit.Circuit(operations)
		block = _block_with_ancillas_zero(circuit, num_controls + 1)
		reference = vqubit.unitary(
			vqubit.ControlledGate(vqubit.MatrixGate(matrix), num_controls=num_controls, control_values=control_values)
		)
		assert vqubit.linalg.allclose_up_to_global_phase(block, reference, atol=1e-5)


def test_decompose_multi_controlled_rotation_validates_control_values() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	with pytest.raises(ValueError):
		vqubit.decompose_multi_controlled_rotation(vqubit.unitary(vqubit.H), [q0, q1], q2, control_values=[True])


@pytest.mark.parametrize("use_sqrt_iswap_inv", [False, True])
def test_two_qubit_matrix_to_sqrt_iswap_operations(use_sqrt_iswap_inv: bool) -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(3)
	for _ in range(3):
		matrix = random_unitary(4, rng)
		operations = vqubit.two_qubit_matrix_to_sqrt_iswap_operations(
			matrix, qubits=(q0, q1), use_sqrt_iswap_inv=use_sqrt_iswap_inv
		)
		_assert_implements(operations, matrix)


@pytest.mark.parametrize("use_sqrt_iswap_inv", [False, True])
def test_two_qubit_matrix_to_sqrt_iswap_operations_only_uses_target_gates(use_sqrt_iswap_inv: bool) -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	matrix = random_unitary(4, random.Random(4))
	operations = vqubit.two_qubit_matrix_to_sqrt_iswap_operations(
		matrix, qubits=(q0, q1), use_sqrt_iswap_inv=use_sqrt_iswap_inv
	)
	for operation in operations:
		if len(operation.qubits) == 2:
			assert isinstance(operation.gate, vqubit.ISwapPowGate)
			assert abs(abs(operation.gate.exponent) - 0.5) < 1e-8
			if use_sqrt_iswap_inv:
				assert operation.gate.exponent < 0


def test_two_qubit_matrix_to_sqrt_iswap_operations_rejects_wrong_shape() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	with pytest.raises((ValueError, RuntimeError)):
		vqubit.two_qubit_matrix_to_sqrt_iswap_operations(vqubit.unitary(vqubit.H), qubits=(q0, q1))


def test_two_qubit_matrix_to_ion_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(5)
	for _ in range(2):
		matrix = random_unitary(4, rng)
		operations = vqubit.two_qubit_matrix_to_ion_operations(matrix, qubits=(q0, q1))
		_assert_implements(operations, matrix)
		assert any(isinstance(operation.gate, vqubit.MSGate) for operation in operations)


def test_two_qubit_matrix_to_ion_operations_rejects_non_unitary_matrices() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	non_unitary = [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 3, 0], [0, 0, 0, 4]]
	with pytest.raises(ValueError, match="4x4 unitary"):
		vqubit.two_qubit_matrix_to_ion_operations(non_unitary, qubits=(q0, q1))


@pytest.mark.parametrize("allow_partial_czs", [False, True])
def test_two_qubit_matrix_to_diagonal_and_cz_operations(allow_partial_czs: bool) -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	matrix = random_unitary(4, random.Random(6))
	diagonal, operations = vqubit.two_qubit_matrix_to_diagonal_and_cz_operations(
		q0, q1, matrix, allow_partial_czs=allow_partial_czs
	)
	assert diagonal.num_qubits == 2
	# The contract is ``matrix ≈ Circuit(operations) @ diagonal``: the diagonal
	# is now the real one, so the operations alone do not implement the matrix.
	combined = vqubit.dot(vqubit.unitary(vqubit.Circuit(operations)), vqubit.unitary(diagonal))
	assert vqubit.linalg.allclose_up_to_global_phase(combined, matrix, atol=1e-5)


def test_two_qubit_matrix_to_cz_isometry_preserves_the_isometry_columns() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	full = random_unitary(4, random.Random(7))
	isometry = [[full[i][j] for j in range(2)] for i in range(4)]
	operations = vqubit.two_qubit_matrix_to_cz_isometry(q0, q1, isometry)
	matrix = vqubit.unitary(vqubit.Circuit(operations))
	phase = None
	for row in range(4):
		for column in range(2):
			if abs(isometry[row][column]) > 1e-9:
				phase = matrix[row][column] / isometry[row][column]
				break
		if phase is not None:
			break
	assert phase is not None
	for row in range(4):
		for column in range(2):
			assert abs(matrix[row][column] - phase * isometry[row][column]) < 1e-6


def test_decompose_cphase_into_two_fsim_implements_the_cz() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	fsim = vqubit.FSimGate(math.pi / 4, 0.0)
	operations = vqubit.decompose_cphase_into_two_fsim(vqubit.CZ, fsim_gate=fsim, qubits=(q0, q1))
	# Two copies of the fSim, with single-qubit rotations around them.
	assert [operation.gate for operation in operations if len(operation.qubits) == 2] == [fsim, fsim]
	_assert_implements(operations, vqubit.unitary(vqubit.CZ))
	# The default qubits are a fresh pair.
	assert vqubit.decompose_cphase_into_two_fsim(vqubit.CZ, fsim_gate=fsim)[0].qubits[0] == vqubit.LineQubit(0)


def test_decompose_cphase_into_two_fsim_rejects_what_it_cannot_do() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	# sin(theta)^2 == sin(phi/2)^2 leaves the construction undefined.
	with pytest.raises(ValueError, match="too close"):
		vqubit.decompose_cphase_into_two_fsim(
			vqubit.CZ, fsim_gate=vqubit.FSimGate(math.pi / 4, math.pi / 2), qubits=(q0, q1)
		)
	# A parametrized gate cannot be decomposed numerically.
	with pytest.raises(ValueError, match="parametrized"):
		vqubit.decompose_cphase_into_two_fsim(
			vqubit.CZ ** vqubit.Parameter("t"),  # ty: ignore[unsupported-operator]
			fsim_gate=vqubit.FSimGate(0.3, 0.4),
			qubits=(q0, q1),
		)
	with pytest.raises(ValueError, match="pair of qubits"):
		vqubit.decompose_cphase_into_two_fsim(vqubit.CZ, fsim_gate=vqubit.FSimGate(0.3, 0.4), qubits=(q0,))


def test_compute_cphase_exponents_for_fsim_decomposition() -> None:
	# The intervals come from the two angles of the fSim: each bounds the phase,
	# so the valid exponents are the sine's hump and its mirror.
	intervals = vqubit.compute_cphase_exponents_for_fsim_decomposition(vqubit.FSimGate(0.3, 0.4))
	assert isinstance(intervals, list)
	assert all(len(interval) == 2 and interval[0] < interval[1] for interval in intervals)
	assert all(start >= 0 and stop <= 2 for start, stop in intervals)
	# A degenerate fSim (theta = phi = 0) cannot decompose any CZ: both bounds
	# collapse and the intervals are empty, as in Cirq.
	assert vqubit.compute_cphase_exponents_for_fsim_decomposition(vqubit.FSimGate(0.0, 0.0)) == []
	with pytest.raises(ValueError, match="Symbolic"):
		vqubit.compute_cphase_exponents_for_fsim_decomposition(vqubit.FSimGate(vqubit.Parameter("t"), 0.0))


def test_decompose_two_qubit_interaction_into_four_fsim_gates() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	fsim = vqubit.FSimGate(math.pi / 2, 0.0)
	circuit = vqubit.decompose_two_qubit_interaction_into_four_fsim_gates(vqubit.ISWAP, fsim_gate=fsim, qubits=(q0, q1))
	operations = list(circuit.all_operations())
	# Exactly four copies of the requested gate, and the whole thing is an iSWAP.
	assert [operation.gate for operation in operations if operation.gate == fsim] == [fsim] * 4
	assert up_to_phase(vqubit.unitary(circuit), vqubit.unitary(vqubit.ISWAP))

	# An ISwapPowGate is mapped onto the fSim it parametrises, and the default
	# qubits are a fresh pair.
	by_power = vqubit.decompose_two_qubit_interaction_into_four_fsim_gates(vqubit.ISWAP, fsim_gate=vqubit.ISWAP)
	assert len([op for op in by_power.all_operations() if isinstance(op.gate, vqubit.ISwapPowGate)]) == 4
	assert up_to_phase(vqubit.unitary(by_power), vqubit.unitary(vqubit.ISWAP))

	# Gates outside the supported window are rejected.
	with pytest.raises(ValueError, match="3π/8"):
		vqubit.decompose_two_qubit_interaction_into_four_fsim_gates(vqubit.ISWAP, fsim_gate=vqubit.FSimGate(0.3, 0.0))
	with pytest.raises(ValueError, match="phi"):
		vqubit.decompose_two_qubit_interaction_into_four_fsim_gates(
			vqubit.ISWAP, fsim_gate=vqubit.FSimGate(math.pi / 2, math.pi / 2)
		)


def test_prepare_two_qubit_state_using_cz() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	state = [1 / math.sqrt(2), 0, 0, 1j / math.sqrt(2)]
	operations = vqubit.prepare_two_qubit_state_using_cz(state, [q0, q1])
	column = [vqubit.unitary(vqubit.Circuit(operations))[i][0] for i in range(4)]
	_assert_vectors_close_up_to_phase(column, state)


def test_prepare_two_qubit_state_using_iswap_and_sqrt_iswap_alias() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	state = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	iswap_operations = vqubit.prepare_two_qubit_state_using_iswap(state, [q0, q1])
	sqrt_iswap_operations = vqubit.prepare_two_qubit_state_using_sqrt_iswap(state, [q0, q1])
	assert sqrt_iswap_operations == iswap_operations
	column = [vqubit.unitary(vqubit.Circuit(iswap_operations))[i][0] for i in range(4)]
	_assert_vectors_close_up_to_phase(column, state)


@pytest.mark.parametrize("use_sqrt_iswap_inv", [False, True])
def test_parameterized_2q_op_to_sqrt_iswap_operations(use_sqrt_iswap_inv: bool) -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gates = [
		vqubit.CZPowGate(exponent=0.3),
		vqubit.SwapPowGate(exponent=0.4),
		vqubit.ISwapPowGate(exponent=0.6),
		vqubit.FSimGate(0.3, 0.4),
	]
	for gate in gates:
		operations = vqubit.parameterized_2q_op_to_sqrt_iswap_operations(
			gate.on(q0, q1), use_sqrt_iswap_inv=use_sqrt_iswap_inv
		)
		_assert_implements(operations, vqubit.unitary(gate))


def test_parameterized_2q_op_returns_not_implemented_for_other_gates() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.parameterized_2q_op_to_sqrt_iswap_operations(vqubit.CNOT(q0, q1)) is NotImplemented


def test_parameterized_2q_op_for_identity_fsim_is_empty() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operations = vqubit.parameterized_2q_op_to_sqrt_iswap_operations(vqubit.FSimGate(0.0, 0.0).on(q0, q1))
	assert operations == []


def test_sqrt_iswap_target_gateset_optimizes_two_qubit_gates() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	matrix = random_unitary(4, random.Random(8))
	gateset = vqubit.SqrtIswapTargetGateset()
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.MatrixGate(matrix).on(q0, q1))
	compiled = gateset.optimize_circuit(circuit)
	assert not any(isinstance(operation.gate, vqubit.MatrixGate) for operation in compiled.all_operations())
	assert vqubit.linalg.allclose_up_to_global_phase(vqubit.unitary(compiled), vqubit.unitary(circuit), atol=1e-5)


def test_sqrt_iswap_target_gateset_metadata_and_keeps_compatible_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gateset = vqubit.SqrtIswapTargetGateset(atol=1e-8)
	assert gateset.num_qubits == 2
	assert len(gateset.entanglers) == 1
	assert isinstance(gateset.entanglers[0], vqubit.ISwapPowGate)
	assert gateset.preprocess_transformers() == []
	assert gateset.postprocess_transformers() == []
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.SQRT_ISWAP(q0, q1), vqubit.measure(q1, key="m"))
	assert list(gateset.optimize_circuit(circuit).all_operations()) == list(circuit.all_operations())
	matrix = random_unitary(4, random.Random(9))
	_assert_implements(gateset.decompose_two_qubit(matrix, (q0, q1)), matrix)


def test_two_qubit_gate_product_tabulation_result_implements_the_gate() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	tabulation = vqubit.two_qubit_gate_product_tabulation(vqubit.SQRT_ISWAP, 0.01)
	assert tabulation.base_gate is vqubit.SQRT_ISWAP
	assert tabulation.max_infidelity == 0.01
	matrix = vqubit.unitary(vqubit.CNOT)
	result = tabulation(vqubit.MatrixGate(matrix), [q0, q1])
	assert result.success
	assert "success=True" in repr(result)
	_assert_implements(result.operations, matrix)


def _build_tableau(num_qubits: int, steps: list[tuple[str, tuple[int, ...]]]) -> vqubit.CliffordTableau:
	tableau = vqubit.CliffordTableau(num_qubits)
	for method, targets in steps:
		getattr(tableau, method)(*targets)
	return tableau


def _simulated_state(circuit: vqubit.Circuit, qubits: tuple[vqubit.LineQubit, ...]) -> Sequence[complex]:
	return vqubit.CliffordSimulator().simulate(circuit, qubit_order=qubits).final_state.state_vector()


def test_decompose_clifford_tableau_to_operations_reproduces_known_states() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	cases: list[tuple[int, list[tuple[str, tuple[int, ...]]], list[vqubit.Operation]]] = [
		(2, [], []),
		(2, [("apply_h", (0,))], [vqubit.H(q0)]),
		(2, [("apply_x", (0,))], [vqubit.X(q0)]),
		(2, [("apply_y", (1,))], [vqubit.Y(q1)]),
		(2, [("apply_z", (1,))], [vqubit.Z(q1)]),
		(2, [("apply_cx", (0, 1))], [vqubit.CNOT(q0, q1)]),
		# The two CNOTs make a SWAP, which exercises the swap branch.
		(
			2,
			[("apply_cx", (0, 1)), ("apply_cx", (1, 0))],
			[vqubit.CNOT(q0, q1), vqubit.CNOT(q1, q0)],
		),
		(2, [("apply_cz", (0, 1))], [vqubit.CZ(q0, q1)]),
		(2, [("apply_h", (0,)), ("apply_cx", (0, 1))], [vqubit.H(q0), vqubit.CNOT(q0, q1)]),
		(
			3,
			[("apply_h", (0,)), ("apply_cx", (0, 1)), ("apply_cz", (1, 2))],
			[vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.CZ(q1, q2)],
		),
	]
	for num_qubits, steps, reference_operations in cases:
		qubits = vqubit.LineQubit.range(num_qubits)
		tableau = _build_tableau(num_qubits, steps)
		operations = vqubit.decompose_clifford_tableau_to_operations(qubits, tableau)
		decomposed = vqubit.Circuit(operations)
		assert set(decomposed.all_qubits()) <= set(qubits)
		_assert_states_close(
			_simulated_state(decomposed, qubits), _simulated_state(vqubit.Circuit(reference_operations), qubits)
		)


def test_decompose_clifford_tableau_to_operations_validates_its_input() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError, match="CliffordTableau"):
		vqubit.decompose_clifford_tableau_to_operations([q0], object())
	with pytest.raises(ValueError, match="number of qubits"):
		vqubit.decompose_clifford_tableau_to_operations([q0], vqubit.CliffordTableau(2))


def test_decompose_multi_controlled_rotation_defaults_to_active_high_controls() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	matrix = vqubit.unitary(vqubit.H)
	operations = vqubit.decompose_multi_controlled_rotation(matrix, [q0, q1], q2)
	circuit = vqubit.Circuit(operations)
	reference = vqubit.unitary(vqubit.ControlledGate(vqubit.MatrixGate(matrix), num_controls=2))
	assert vqubit.linalg.allclose_up_to_global_phase(_block_with_ancillas_zero(circuit, 3), reference, atol=1e-6)


@pytest.mark.parametrize("num_qubits", [1, 2, 3])
def test_quantum_shannon_decomposition(num_qubits: int) -> None:
	qubits = vqubit.LineQubit.range(num_qubits)
	matrix = random_unitary(1 << num_qubits, random.Random(10 + num_qubits))
	operations = vqubit.quantum_shannon_decomposition(qubits, matrix)
	circuit = vqubit.Circuit(operations)
	block = _block_with_ancillas_zero(circuit, num_qubits)
	assert vqubit.linalg.allclose_up_to_global_phase(block, matrix, atol=1e-4)


def test_quantum_shannon_decomposition_of_the_identity_is_empty() -> None:
	qubits = vqubit.LineQubit.range(2)
	matrix = [[1.0 if row == column else 0.0 for column in range(4)] for row in range(4)]
	assert vqubit.quantum_shannon_decomposition(qubits, matrix) == []


def test_three_qubit_matrix_to_operations() -> None:
	qubits = vqubit.LineQubit.range(3)
	matrix = random_unitary(8, random.Random(15))
	operations = vqubit.three_qubit_matrix_to_operations(qubits, matrix)
	circuit = vqubit.Circuit(operations)
	block = _block_with_ancillas_zero(circuit, 3)
	assert vqubit.linalg.allclose_up_to_global_phase(block, matrix, atol=1e-4)


def test_two_qubit_compilation_target_gateset_is_abstract() -> None:
	gateset = vqubit.TwoQubitCompilationTargetGateset(vqubit.SQRT_ISWAP)
	q0, q1 = vqubit.LineQubit.range(2)
	with pytest.raises(NotImplementedError):
		gateset.decompose_two_qubit(vqubit.unitary(vqubit.CNOT), (q0, q1))
