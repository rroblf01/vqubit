"""Tests for the analytical decompositions and compilation gatesets."""

from __future__ import annotations

import math
import random

import vqubit
from tests.helpers import random_unitary, up_to_phase

_PAULI_GATE = {vqubit.Pauli.X: vqubit.X, vqubit.Pauli.Y: vqubit.Y, vqubit.Pauli.Z: vqubit.Z}


def _reconstruct_pauli_rotations(rotations, qubit):
	operations = [(_PAULI_GATE[pauli] ** turns).on(qubit) for pauli, turns in rotations]
	return vqubit.unitary(vqubit.Circuit(operations))


def test_single_qubit_matrix_to_pauli_rotations() -> None:
	q0 = vqubit.LineQubit(0)
	rng = random.Random(0)
	for _ in range(5):
		matrix = random_unitary(2, rng)
		rotations = vqubit.single_qubit_matrix_to_pauli_rotations(matrix)
		assert up_to_phase(_reconstruct_pauli_rotations(rotations, q0), matrix, tol=1e-6)


def test_single_qubit_op_to_framed_phase_form() -> None:
	q0 = vqubit.LineQubit(0)
	_phase, phased_x, z_gate = vqubit.single_qubit_op_to_framed_phase_form(vqubit.H(q0))
	matrix = vqubit.dot(vqubit.unitary(z_gate), vqubit.unitary(phased_x))
	assert up_to_phase(matrix, vqubit.unitary(vqubit.H), tol=1e-6)


def test_decompose_multi_controlled_x() -> None:
	for n in (1, 2):
		qubits = vqubit.LineQubit.range(n + 1)
		operations = vqubit.decompose_multi_controlled_x(qubits[:n], qubits[n])
		circuit = vqubit.Circuit(operations)
		reference = vqubit.unitary(vqubit.ControlledGate(vqubit.X, num_controls=n))
		assert up_to_phase(vqubit.unitary(circuit), reference, tol=1e-6)


def test_decompose_multi_controlled_rotation() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(1)
	for _ in range(5):
		matrix = random_unitary(2, rng)
		operations = vqubit.decompose_multi_controlled_rotation(matrix, [q0], q1)
		reference = vqubit.unitary(vqubit.ControlledGate(vqubit.MatrixGate(matrix), num_controls=1))
		assert up_to_phase(vqubit.unitary(vqubit.Circuit(operations)), reference, tol=1e-6)


def test_two_qubit_matrix_to_sqrt_iswap() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(2)
	for _ in range(5):
		matrix = random_unitary(4, rng)
		operations = vqubit.two_qubit_matrix_to_sqrt_iswap_operations(matrix, qubits=(q0, q1))
		assert up_to_phase(vqubit.unitary(vqubit.Circuit(operations)), matrix, tol=1e-5)


def test_two_qubit_matrix_to_ion() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(3)
	for _ in range(5):
		matrix = random_unitary(4, rng)
		operations = vqubit.two_qubit_matrix_to_ion_operations(matrix, qubits=(q0, q1))
		assert up_to_phase(vqubit.unitary(vqubit.Circuit(operations)), matrix, tol=1e-5)


def test_parameterized_2q_op_to_sqrt_iswap() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	for turns in (1.0, 0.5, -0.3):
		operation = vqubit.CZPowGate(exponent=turns).on(q0, q1)
		operations = vqubit.parameterized_2q_op_to_sqrt_iswap_operations(operation)
		assert up_to_phase(
			vqubit.unitary(vqubit.Circuit(operations)), vqubit.unitary(vqubit.CZPowGate(exponent=turns)), tol=1e-6
		)


def test_quantum_shannon_decomposition() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(4)
	for _ in range(3):
		matrix = random_unitary(4, rng)
		operations = vqubit.quantum_shannon_decomposition([q0, q1], matrix)
		assert up_to_phase(vqubit.unitary(vqubit.Circuit(operations)), matrix, tol=1e-5)


def test_three_qubit_matrix_to_operations() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	rng = random.Random(5)
	matrix = random_unitary(8, rng)
	circuit = vqubit.Circuit(vqubit.three_qubit_matrix_to_operations([q0, q1, q2], matrix))
	full = vqubit.unitary(circuit)
	extra = len(circuit.all_qubits()) - 3
	stride = 1 << extra
	block = [[full[a * stride][b * stride] for b in range(8)] for a in range(8)]
	assert up_to_phase(block, matrix, tol=1e-4)


def test_sqrt_iswap_target_gateset() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(6)
	matrix = random_unitary(4, rng)
	circuit = vqubit.Circuit(vqubit.MatrixGate(matrix).on(q0, q1))
	gateset = vqubit.SqrtIswapTargetGateset()
	compiled = gateset.optimize_circuit(circuit)
	assert up_to_phase(vqubit.unitary(compiled), matrix, tol=1e-5)


def test_two_qubit_gate_product_tabulation() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	tabulation = vqubit.two_qubit_gate_product_tabulation(vqubit.SQRT_ISWAP, 0.01)
	assert isinstance(tabulation, vqubit.TwoQubitGateTabulation)
	result = tabulation(vqubit.MatrixGate(vqubit.unitary(vqubit.CNOT)), [q0, q1])
	assert isinstance(result, vqubit.TwoQubitGateTabulationResult)
	assert result.success


def test_prepare_two_qubit_state_using_cz() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	state = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	operations = vqubit.prepare_two_qubit_state_using_cz(state, [q0, q1])
	prepared = vqubit.unitary(vqubit.Circuit(operations))
	column = [prepared[i][0] for i in range(4)]
	overlap = sum(column[i].conjugate() * state[i] for i in range(4))
	assert abs(abs(overlap) - 1) < 1e-6


def test_other_compilation_helpers_run() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert isinstance(vqubit.compute_cphase_exponents_for_fsim_decomposition(vqubit.FSimGate(0.3, 0.4)), tuple)
	assert len(vqubit.decompose_cphase_into_two_fsim(q0, q1)) > 0
	assert len(vqubit.decompose_two_qubit_interaction_into_four_fsim_gates(q0, q1)) == 4
	diagonal, operations = vqubit.two_qubit_matrix_to_diagonal_and_cz_operations(
		vqubit.unitary(vqubit.CNOT), qubits=(q0, q1)
	)
	assert diagonal.num_qubits == 2
	assert len(operations) > 0
	assert len(vqubit.two_qubit_matrix_to_cz_isometry(q0, q1, [[1, 0], [0, 1], [0, 0], [0, 0]])) > 0
