"""Extra tests for the gates in :mod:`vqubit._gates_more`.

The existing ``test_gates_extra.py`` checks the headline identities. These tests
pin down the remaining public surface: parameter resolution, validation,
diagram symbols, powers and the factory helpers.
"""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close, up_to_phase
from vqubit import testing


def test_phased_iswap_properties_and_powers() -> None:
	gate = vqubit.PhasedISwapPowGate(phase_exponent=0.1, exponent=0.5, global_shift=0.2)
	assert gate.num_qubits == 2
	assert gate.phase_exponent == pytest.approx(0.1)
	assert gate.exponent == pytest.approx(0.5)
	assert gate.global_shift == pytest.approx(0.2)
	squared = gate**2
	assert isinstance(squared, vqubit.PhasedISwapPowGate)
	assert squared.exponent == pytest.approx(1.0)
	assert squared.phase_exponent == pytest.approx(0.1)
	assert squared.global_shift == pytest.approx(0.2)
	# Squaring the gate squares its matrix.
	assert_matrix_close(vqubit.unitary(squared), vqubit.dot(vqubit.unitary(gate), vqubit.unitary(gate)))


def test_phased_iswap_inverse_undoes_the_gate() -> None:
	gate = vqubit.PhasedISwapPowGate(phase_exponent=0.3, exponent=1.5)
	inverse = vqubit.inverse(gate)
	assert isinstance(inverse, vqubit.PhasedISwapPowGate)
	assert inverse.exponent == pytest.approx(-1.5)
	assert_matrix_close(
		vqubit.dot(vqubit.unitary(gate), vqubit.unitary(inverse)), vqubit.unitary(vqubit.IdentityGate(2))
	)
	assert vqubit.inverse(inverse) == gate


def test_phased_iswap_is_an_iswap_for_zero_phase_exponent() -> None:
	# With no Z conjugation the matrix is the plain iSWAP matrix.
	assert_matrix_close(
		vqubit.unitary(vqubit.PhasedISwapPowGate(phase_exponent=0.0, exponent=1.0)),
		[
			[1 + 0j, 0j, 0j, 0j],
			[0j, 0j, 1j, 0j],
			[0j, 1j, 0j, 0j],
			[0j, 0j, 0j, 1 + 0j],
		],
	)
	# A global shift multiplies the whole matrix by a phase.
	shifted = vqubit.unitary(vqubit.PhasedISwapPowGate(phase_exponent=0.0, exponent=1.0, global_shift=0.5))
	assert_matrix_close(shifted, [[1j * entry for entry in row] for row in vqubit.unitary(vqubit.ISWAP)])


def test_phased_iswap_is_parameterized_and_resolves() -> None:
	parameter = vqubit.Parameter("t")
	gate = vqubit.PhasedISwapPowGate(phase_exponent=parameter, exponent=parameter)
	assert vqubit.is_parameterized(gate)
	assert not vqubit.has_unitary(gate)
	assert not vqubit.has_unitary(vqubit.PhasedISwapPowGate(exponent=parameter))
	resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({parameter: 0.25}))
	assert isinstance(resolved, vqubit.PhasedISwapPowGate)
	assert resolved.phase_exponent == pytest.approx(0.25)
	assert resolved.exponent == pytest.approx(0.25)
	assert not vqubit.is_parameterized(resolved)


def test_phased_iswap_symbols_and_repr() -> None:
	gate = vqubit.PhasedISwapPowGate(phase_exponent=0.25, exponent=0.5)
	assert vqubit.circuit_diagram_info(gate).wire_symbols == ("PhISwap", "PhISwap")
	assert "phase_exponent=0.25" in str(gate)
	testing.assert_equivalent_repr(gate)


def test_givens_is_zero_at_zero_angle() -> None:
	assert isinstance(vqubit.givens(0.0), vqubit.PhasedISwapPowGate)
	assert_matrix_close(vqubit.unitary(vqubit.givens(0.0)), vqubit.unitary(vqubit.IdentityGate(2)))
	# givens(angle) is the PhasedISwap power with exponent 2*angle/pi.
	assert_matrix_close(
		vqubit.unitary(vqubit.givens(0.4)), vqubit.unitary(vqubit.PhasedISwapPowGate(exponent=0.8 / math.pi))
	)


def test_diagonal_gate_applies_phases_on_basis_states() -> None:
	angles = [0.0, math.pi / 2, math.pi, 3 * math.pi / 2]
	gate = vqubit.DiagonalGate(angles)
	assert gate.num_qubits == 2
	assert gate.diag_angles_radians == tuple(angles)
	matrix = vqubit.unitary(gate)
	for i, angle in enumerate(angles):
		for j in range(4):
			expected = complex(math.cos(angle), math.sin(angle)) if i == j else 0j
			assert matrix[i][j] == pytest.approx(expected)
	assert vqubit.circuit_diagram_info(gate).wire_symbols == ("diag", "diag")
	assert "DiagonalGate" in str(gate)


def test_diagonal_gate_validates_its_size() -> None:
	with pytest.raises(ValueError):
		vqubit.DiagonalGate([])
	with pytest.raises(ValueError):
		vqubit.TwoQubitDiagonalGate([0.0, 0.0])
	with pytest.raises(ValueError):
		vqubit.ThreeQubitDiagonalGate([0.0] * 4)
	# The subclasses are strict about the number of angles they name.
	assert vqubit.TwoQubitDiagonalGate([0.0] * 4).num_qubits == 2
	assert vqubit.ThreeQubitDiagonalGate([0.0] * 8).num_qubits == 3


def test_three_qubit_diagonal_gate_can_build_ccz() -> None:
	assert_matrix_close(
		vqubit.unitary(vqubit.ThreeQubitDiagonalGate([0.0] * 7 + [math.pi])),
		vqubit.unitary(vqubit.CCZ),
	)


def test_diagonal_gate_is_parameterized_and_resolves() -> None:
	parameter = vqubit.Parameter("p")
	gate = vqubit.DiagonalGate([0.0, parameter])
	assert vqubit.is_parameterized(gate)
	with pytest.raises(TypeError):
		vqubit.unitary(gate)
	resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({parameter: 0.4}))
	assert resolved == vqubit.DiagonalGate([0.0, 0.4])
	assert not vqubit.is_parameterized(resolved)


def test_qubit_permutation_maps_basis_states_to_new_positions() -> None:
	# permutation[i] is the position qubit i moves to. A lone excitation on
	# qubit 0 (|100>, index 4) moves to position 1 (|010>, index 2).
	gate = vqubit.QubitPermutationGate([1, 2, 0])
	assert gate.num_qubits == 3
	assert gate.permutation == (1, 2, 0)
	matrix = vqubit.unitary(gate)
	assert matrix[2][4] == pytest.approx(1)
	# Excitations are preserved by the permutation.
	assert all(sum(1 for entry in row if abs(entry) > 1e-9) <= 1 for row in matrix)
	assert vqubit.is_unitary(matrix)
	assert vqubit.circuit_diagram_info(gate).wire_symbols == ("[0>1]", "[1>2]", "[2>0]")


def test_qubit_permutation_validates_and_supports_the_identity() -> None:
	with pytest.raises(ValueError):
		vqubit.QubitPermutationGate([])
	with pytest.raises(ValueError):
		vqubit.QubitPermutationGate([0, 0])
	with pytest.raises(ValueError):
		vqubit.QubitPermutationGate([0, 2])
	assert_matrix_close(vqubit.unitary(vqubit.QubitPermutationGate([0, 1])), vqubit.unitary(vqubit.IdentityGate(2)))
	assert vqubit.QubitPermutationGate([1, 0]) == vqubit.QubitPermutationGate([1, 0])
	assert vqubit.QubitPermutationGate([1, 0]) != vqubit.QubitPermutationGate([0, 1])
	testing.assert_equivalent_repr(vqubit.QubitPermutationGate([0, 2, 1, 3]))


def test_fourier_transform_of_one_qubit_is_hadamard() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.QuantumFourierTransformGate(1)), vqubit.unitary(vqubit.H))
	expected_two = [
		[0.5 + 0j, 0.5 + 0j, 0.5 + 0j, 0.5 + 0j],
		[0.5 + 0j, 0.5j, -0.5 + 0j, -0.5j],
		[0.5 + 0j, -0.5 + 0j, 0.5 + 0j, -0.5 + 0j],
		[0.5 + 0j, -0.5j, -0.5 + 0j, 0.5j],
	]
	assert_matrix_close(vqubit.unitary(vqubit.QuantumFourierTransformGate(2)), expected_two)


def test_fourier_transform_without_reverse_reorders_the_output_rows() -> None:
	normal = vqubit.unitary(vqubit.QuantumFourierTransformGate(2))
	without_reverse = vqubit.unitary(vqubit.QuantumFourierTransformGate(2, without_reverse=True))
	# Reversing two bits swaps output rows 1 and 2 and leaves 0 and 3 in place.
	assert_matrix_close([without_reverse[i] for i in (0, 2, 1, 3)], normal)
	assert vqubit.QuantumFourierTransformGate(2) != vqubit.QuantumFourierTransformGate(2, without_reverse=True)
	assert vqubit.circuit_diagram_info(vqubit.QuantumFourierTransformGate(3)).wire_symbols == ("QFT", "QFT", "QFT")


def test_qft_factory_inverse_is_the_dagger() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	forward = vqubit.qft(q0, q1)
	inverse = vqubit.qft(q0, q1, inverse=True)
	assert isinstance(forward.gate, vqubit.QuantumFourierTransformGate)
	assert isinstance(inverse.gate, vqubit.MatrixGate)
	assert vqubit.num_qubits(forward) == 2
	assert_matrix_close(
		vqubit.dot(vqubit.unitary(forward), vqubit.unitary(inverse)), vqubit.unitary(vqubit.IdentityGate(2))
	)
	# without_reverse is forwarded to the gate.
	reversed_qft = vqubit.qft(q0, q1, without_reverse=True)
	assert isinstance(reversed_qft.gate, vqubit.QuantumFourierTransformGate)
	assert_matrix_close(
		vqubit.unitary(reversed_qft),
		vqubit.unitary(vqubit.QuantumFourierTransformGate(2, without_reverse=True)),
	)


def test_factories_forward_their_parameters() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.cphase(0.3) == vqubit.CZPowGate(exponent=0.3 / math.pi)
	assert vqubit.cphase(0.3).exponent == pytest.approx(0.3 / math.pi)
	assert vqubit.ms(0.4).rads == pytest.approx(0.4)
	assert vqubit.riswap(0.4) == vqubit.ISwapPowGate(exponent=0.8 / math.pi)
	assert vqubit.q(2) == vqubit.LineQubit(2)
	with pytest.raises(ValueError):
		vqubit.q(-1)
	assert vqubit.wait(12.5, 2).duration_nanos == 12.5
	assert vqubit.wait(12.5, 2).num_qubits == 2
	assert vqubit.reset_each() == ()
	assert vqubit.identity_each() == ()
	assert vqubit.identity_each(q0)[0].qubits == (q0,)


def test_global_phase_operation_produces_a_1x1_unitary() -> None:
	for coefficient in (1j, -1 + 0j, 2 + 0j):
		operation = vqubit.global_phase_operation(coefficient)
		assert operation.qubits == ()
		assert_matrix_close(vqubit.unitary(operation), [[coefficient / abs(coefficient)]])


def test_parallel_gate_applies_the_sub_gate_to_two_groups() -> None:
	gate = vqubit.ParallelGate(vqubit.H)
	assert gate.num_qubits == 2
	assert gate.sub_gate == vqubit.H
	assert_matrix_close(vqubit.unitary(gate), vqubit.kron(vqubit.unitary(vqubit.H), vqubit.unitary(vqubit.H)))
	assert vqubit.circuit_diagram_info(gate).wire_symbols == ("H", "H")
	assert vqubit.ParallelGate(vqubit.H) == vqubit.ParallelGate(vqubit.H)
	assert vqubit.ParallelGate(vqubit.H) != vqubit.ParallelGate(vqubit.X)
	assert "ParallelGate" in str(gate)
	testing.assert_equivalent_repr(gate)


def test_parallel_gate_rejects_zero_qubit_sub_gates_and_non_unitary_ones() -> None:
	with pytest.raises(ValueError):
		vqubit.ParallelGate(vqubit.IdentityGate(0))
	channel_gate = vqubit.ParallelGate(vqubit.bit_flip(0.1))
	assert channel_gate.num_qubits == 2
	assert not vqubit.has_unitary(channel_gate)
	assert vqubit.circuit_diagram_info(channel_gate).wire_symbols == ("BF", "BF")


def test_parallel_gate_op_returns_one_operation_over_all_qubits() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.parallel_gate_op(vqubit.X, q0, q1)
	assert operation.qubits == (q0, q1)
	assert isinstance(operation.gate, vqubit.ParallelGate)
	assert vqubit.num_qubits(operation) == 2
	assert_matrix_close(vqubit.unitary(operation), vqubit.kron(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.X)))


def test_phase_gradient_gate_is_a_diagonal_phase_ladder() -> None:
	one_qubit = vqubit.PhaseGradientGate(num_qubits=1, exponent=1)
	assert_matrix_close(vqubit.unitary(one_qubit), vqubit.unitary(vqubit.Z))
	gate = vqubit.PhaseGradientGate(num_qubits=2, exponent=1)
	assert gate.num_qubits == 2
	assert gate.exponent == pytest.approx(1)
	diagonal = [vqubit.unitary(gate)[i][i] for i in range(4)]
	assert diagonal[0] == pytest.approx(1)
	assert diagonal[1] == pytest.approx(1j)
	assert diagonal[2] == pytest.approx(-1)
	assert diagonal[3] == pytest.approx(-1j)
	assert_matrix_close(
		vqubit.unitary(vqubit.PhaseGradientGate(num_qubits=2, exponent=0)), vqubit.unitary(vqubit.IdentityGate(2))
	)
	assert vqubit.circuit_diagram_info(vqubit.PhaseGradientGate(num_qubits=3, exponent=1)).wire_symbols == (
		"Grad",
		"#2",
		"#3",
	)


def test_phase_gradient_gate_powers_and_parameterization() -> None:
	gate = vqubit.PhaseGradientGate(num_qubits=2, exponent=0.5)
	assert (gate**2).exponent == pytest.approx(1.0)
	assert (gate**0).exponent == pytest.approx(0.0)
	assert_matrix_close(vqubit.unitary(gate**0), vqubit.unitary(vqubit.IdentityGate(2)))

	parameter = vqubit.Parameter("t")
	symbolic = vqubit.PhaseGradientGate(num_qubits=2, exponent=parameter)
	assert vqubit.is_parameterized(symbolic)
	assert not vqubit.has_unitary(symbolic)
	resolved = vqubit.resolve_parameters(symbolic, vqubit.ParamResolver({parameter: 0.5}))
	assert resolved == gate
	assert not vqubit.is_parameterized(resolved)
	assert gate != vqubit.PhaseGradientGate(num_qubits=3, exponent=0.5)
	testing.assert_equivalent_repr(gate)


def test_phased_fsim_reduces_to_fsim_and_identity() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.PhasedFSimGate(0.0)), vqubit.unitary(vqubit.IdentityGate(2)))
	# With zeta = chi = gamma = 0 the gate is FSimGate(theta, phi).
	assert_matrix_close(
		vqubit.unitary(vqubit.PhasedFSimGate(0.5, phi=0.3)),
		vqubit.unitary(vqubit.FSimGate(0.5, 0.3)),
	)
	# phi only phases the |11> amplitude.
	matrix = vqubit.unitary(vqubit.PhasedFSimGate(0.0, phi=0.5))
	assert matrix[3][3] == pytest.approx(complex(math.cos(0.5), -math.sin(0.5)))


def test_phased_fsim_properties_and_angle_canonicalization() -> None:
	gate = vqubit.PhasedFSimGate(0.5, 0.1, 0.2, 0.3, 0.4)
	assert gate.num_qubits == 2
	assert (gate.theta, gate.zeta, gate.chi, gate.gamma, gate.phi) == (
		pytest.approx(0.5),
		pytest.approx(0.1),
		pytest.approx(0.2),
		pytest.approx(0.3),
		pytest.approx(0.4),
	)
	# Angles are canonicalized into [-pi, pi).
	assert vqubit.PhasedFSimGate(math.pi + 0.1).theta == pytest.approx(-math.pi + 0.1)
	assert vqubit.PhasedFSimGate(2 * math.pi + 0.1).theta == pytest.approx(0.1)
	assert vqubit.circuit_diagram_info(gate).wire_symbols == ("PhFSim", "PhFSim")
	testing.assert_equivalent_repr(gate)


def test_phased_fsim_from_fsim_rz_roundtrips_the_rz_angles() -> None:
	before = (0.2, -0.3)
	after = (0.4, 0.1)
	gate = vqubit.PhasedFSimGate.from_fsim_rz(0.5, 0.2, before, after)
	assert isinstance(gate, vqubit.PhasedFSimGate)
	assert gate.theta == pytest.approx(0.5)
	assert gate.phi == pytest.approx(0.2)
	# The angles are defined up to a common shift: only their differences are
	# physical. The gate is Rz(after) FSim Rz(before), up to a global phase.
	assert gate.rz_angles_before[0] - gate.rz_angles_before[1] == pytest.approx(before[0] - before[1])
	assert gate.rz_angles_after[0] - gate.rz_angles_after[1] == pytest.approx(after[0] - after[1])

	q0, q1 = vqubit.LineQubit.range(2)
	b0, b1 = before
	a0, a1 = after
	circuit = vqubit.Circuit(
		vqubit.rz(b0)(q0),
		vqubit.rz(b1)(q1),
		vqubit.FSimGate(0.5, 0.2)(q0, q1),
		vqubit.rz(a0)(q0),
		vqubit.rz(a1)(q1),
	)
	assert up_to_phase(vqubit.unitary(circuit), vqubit.unitary(gate))


def test_phased_fsim_is_parameterized_and_resolves() -> None:
	parameter = vqubit.Parameter("t")
	gate = vqubit.PhasedFSimGate(parameter, 0.1, parameter, 0.3, 0.4)
	assert vqubit.is_parameterized(gate)
	with pytest.raises(TypeError):
		vqubit.unitary(gate)
	resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({parameter: 0.2}))
	assert resolved == vqubit.PhasedFSimGate(0.2, 0.1, 0.2, 0.3, 0.4)
	assert not vqubit.is_parameterized(resolved)


def test_uniform_superposition_gate_prepares_a_uniform_state() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gate = vqubit.UniformSuperpositionGate(3, 2)
	assert gate.m_value == 3
	assert gate.num_qubits == 2
	matrix = vqubit.unitary(gate)
	assert len(matrix) == 4
	for row in range(3):
		assert matrix[row][0] == pytest.approx(1 / math.sqrt(3))
	assert matrix[3][0] == pytest.approx(0)
	# The decomposition is non-empty and starts from |00>.
	operations = vqubit.decompose(gate.on(q0, q1))
	assert operations
	assert vqubit.is_unitary(matrix)


def test_uniform_superposition_gate_power_of_two_and_wider_register() -> None:
	# m = 4 spans all states of two qubits, so it is H⊗H.
	assert_matrix_close(
		vqubit.unitary(vqubit.UniformSuperpositionGate(4, 2)),
		vqubit.kron(vqubit.unitary(vqubit.H), vqubit.unitary(vqubit.H)),
	)
	# m = 5 on three qubits fills the first five amplitudes.
	matrix = vqubit.unitary(vqubit.UniformSuperpositionGate(5, 3))
	assert len(matrix) == 8
	for row in range(5):
		assert matrix[row][0] == pytest.approx(1 / math.sqrt(5))
	for row in range(5, 8):
		assert matrix[row][0] == pytest.approx(0)
	assert vqubit.circuit_diagram_info(vqubit.UniformSuperpositionGate(3, 2)).wire_symbols == ("USup", "USup")


def test_uniform_superposition_gate_validates_its_arguments() -> None:
	with pytest.raises(ValueError):
		vqubit.UniformSuperpositionGate(0, 1)
	with pytest.raises(ValueError):
		vqubit.UniformSuperpositionGate(-2, 2)
	with pytest.raises(ValueError):
		vqubit.UniformSuperpositionGate(9, 2)
	with pytest.raises(ValueError):
		vqubit.UniformSuperpositionGate(3.5, 2)  # ty: ignore[invalid-argument-type]
	# num_qubits 3 is wide enough for m = 5 (log2(5) rounds up to 3 bits).
	assert vqubit.UniformSuperpositionGate(5, 3).num_qubits == 3
	with pytest.raises(ValueError):
		vqubit.UniformSuperpositionGate(5, 2)
	assert vqubit.UniformSuperpositionGate(3, 2) != vqubit.UniformSuperpositionGate(3, 3)
	testing.assert_equivalent_repr(vqubit.UniformSuperpositionGate(3, 2))


def test_riswap_is_the_iswap_power_at_twice_the_angle() -> None:
	gate = vqubit.riswap(0.6)
	assert isinstance(gate, vqubit.ISwapPowGate)
	assert up_to_phase(vqubit.unitary(gate), vqubit.unitary(vqubit.ISwapPowGate(exponent=1.2 / math.pi)))
