"""Tests for the extra gates in :mod:`vqubit._gates_more`."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close, up_to_phase


def _unitary(gate: vqubit.Gate) -> list[list[complex]]:
	return vqubit.unitary(gate)


def test_controlled_pow_gates_at_exponent_one() -> None:
	assert_matrix_close(_unitary(vqubit.CXPowGate()), _unitary(vqubit.CNOT))
	assert_matrix_close(_unitary(vqubit.CYPowGate()), _unitary(vqubit.CY))
	assert_matrix_close(_unitary(vqubit.CZPowGate()), _unitary(vqubit.CZ))
	assert_matrix_close(_unitary(vqubit.CCXPowGate()), _unitary(vqubit.TOFFOLI))
	assert_matrix_close(_unitary(vqubit.CZPowGate(exponent=0)), vqubit.unitary(vqubit.IdentityGate(2)))


def test_cnot_pow_squares_to_cnot() -> None:
	half = _unitary(vqubit.CXPowGate(exponent=0.5))
	assert up_to_phase(vqubit.dot(half, half), _unitary(vqubit.CNOT))


def test_ccz_pow_is_diagonal() -> None:
	matrix = _unitary(vqubit.CCZPowGate())
	assert [matrix[i][i] for i in range(8)] == [1] * 7 + [-1]


def test_swap_pow_gates() -> None:
	assert_matrix_close(_unitary(vqubit.SwapPowGate()), _unitary(vqubit.SWAP))
	half = _unitary(vqubit.SwapPowGate(exponent=0.5))
	assert up_to_phase(vqubit.dot(half, half), _unitary(vqubit.SWAP))


def test_iswap_pow_gates() -> None:
	assert_matrix_close(_unitary(vqubit.ISwapPowGate()), _unitary(vqubit.ISWAP))
	assert_matrix_close(
		vqubit.dot(_unitary(vqubit.ISWAP_INV), _unitary(vqubit.ISWAP)), vqubit.unitary(vqubit.IdentityGate(2))
	)
	assert_matrix_close(
		vqubit.dot(_unitary(vqubit.SQRT_ISWAP_INV), _unitary(vqubit.SQRT_ISWAP_INV)), _unitary(vqubit.ISWAP_INV)
	)


def test_fixed_aliases_match() -> None:
	assert_matrix_close(_unitary(vqubit.CCNotPowGate()), _unitary(vqubit.TOFFOLI))
	assert vqubit.CNotPowGate is vqubit.CXPowGate
	assert_matrix_close(
		_unitary(vqubit.CCZ), [[1 if i == j and i < 7 else (-1 if i == j else 0) for j in range(8)] for i in range(8)]
	)


def test_diagonal_gates() -> None:
	# Angles, not entries: DiagonalGate([0, 0, 0, pi]) is diag(1, 1, 1, -1) = CZ.
	assert_matrix_close(_unitary(vqubit.DiagonalGate([0, 0, 0, math.pi])), _unitary(vqubit.CZ))
	assert vqubit.DiagonalGate([0, 0, 0, math.pi]).diag_angles_radians == (0, 0, 0, math.pi)
	assert_matrix_close(
		_unitary(vqubit.TwoQubitDiagonalGate([0.5, 1.0, 1.5, 2.0])),
		_unitary(vqubit.DiagonalGate([0.5, 1.0, 1.5, 2.0])),
	)
	assert vqubit.ThreeQubitDiagonalGate([0] * 8).num_qubits == 3
	with pytest.raises(ValueError):
		vqubit.DiagonalGate([1, 1, 1])


def test_qubit_permutation_gate() -> None:
	gate = vqubit.QubitPermutationGate([1, 0])
	assert gate.num_qubits == 2
	assert vqubit.is_unitary(_unitary(gate))
	# Swapping the two qubits is the SWAP gate.
	assert_matrix_close(_unitary(gate), _unitary(vqubit.SWAP))
	assert vqubit.QubitPermutationGate([1, 2, 0]).num_qubits == 3


def test_fourier_transform_gate() -> None:
	for n in (1, 2, 3):
		assert vqubit.is_unitary(_unitary(vqubit.QuantumFourierTransformGate(n)))
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.qft(q0, q1)
	assert operation.gate.num_qubits == 2


def test_lowercase_factories() -> None:
	assert_matrix_close(_unitary(vqubit.rx(0.3)), _unitary(vqubit.Rx(0.3)))
	assert_matrix_close(_unitary(vqubit.ry(0.3)), _unitary(vqubit.Ry(0.3)))
	assert_matrix_close(_unitary(vqubit.rz(0.3)), _unitary(vqubit.Rz(0.3)))
	assert isinstance(vqubit.ms(0.1), vqubit.MSGate)
	assert_matrix_close(_unitary(vqubit.cphase(math.pi)), _unitary(vqubit.CZ))
	assert vqubit.q(4) == vqubit.LineQubit(4)
	assert isinstance(vqubit.wait(10.0, 2), vqubit.WaitGate)


def test_reset_and_identity_each() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	resets = vqubit.reset_each(q0, q1)
	assert len(resets) == 2
	assert all(op.qubits == (q,) for op, q in zip(resets, (q0, q1), strict=True))
	identities = vqubit.identity_each(q0)
	assert len(identities) == 1
	assert isinstance(identities[0].gate, vqubit.IdentityGate)


def test_global_phase_operation() -> None:
	operation = vqubit.global_phase_operation(1j)
	assert abs(vqubit.unitary(operation)[0][0] - 1j) < 1e-9


def test_phased_iswap_pow_gate() -> None:
	# With phase exponent 0 it reduces to the iSWAP power gate.
	assert up_to_phase(
		_unitary(vqubit.PhasedISwapPowGate(phase_exponent=0.0, exponent=0.5)),
		_unitary(vqubit.ISwapPowGate(exponent=0.5)),
	)
	assert vqubit.is_unitary(_unitary(vqubit.PhasedISwapPowGate()))
	from vqubit import testing

	testing.assert_equivalent_repr(vqubit.PhasedISwapPowGate(phase_exponent=0.25, exponent=0.5))


def test_givens_factory() -> None:
	gate = vqubit.givens(0.7)
	assert isinstance(gate, vqubit.PhasedISwapPowGate)
	assert vqubit.is_unitary(_unitary(gate))
	assert up_to_phase(_unitary(gate), vqubit.unitary(vqubit.PhasedISwapPowGate(exponent=2 * 0.7 / math.pi)))


def test_qid_variants() -> None:
	assert vqubit.NamedQid("a").name == "a"
	assert vqubit.LineQid(3).x == 3
	assert vqubit.GridQid(1, 2).row == 1
	assert vqubit.GridQid(1, 2).col == 2
	assert vqubit.LineQid.range(3) == (vqubit.LineQid(0), vqubit.LineQid(1), vqubit.LineQid(2))
	assert isinstance(vqubit.NamedQid("a"), vqubit.Qid)


def test_riswap_and_phase_gradient() -> None:
	assert up_to_phase(_unitary(vqubit.riswap(0.7)), _unitary(vqubit.ISwapPowGate(exponent=2 * 0.7 / math.pi)))
	gradient = vqubit.PhaseGradientGate(num_qubits=2, exponent=1)
	assert vqubit.is_unitary(_unitary(gradient))
	assert _unitary(gradient)[0][0] == pytest.approx(1)
	assert _unitary(gradient)[2][2] == pytest.approx(-1)
	assert _unitary(gradient**2)[2][2] == pytest.approx(1)


def test_phased_fsim_gate() -> None:
	gate = vqubit.PhasedFSimGate(0.5, 0.1, 0.2, 0.3, 0.4)
	assert vqubit.is_unitary(_unitary(gate))
	rebuilt = vqubit.PhasedFSimGate.from_fsim_rz(0.5, 0.4, gate.rz_angles_before, gate.rz_angles_after)
	assert up_to_phase(_unitary(rebuilt), _unitary(gate))
	assert vqubit.PhasedFSimGate(0.5).theta == pytest.approx(0.5)


def test_uniform_superposition_gate() -> None:
	for m_value, num_qubits in [(1, 1), (2, 2), (3, 2), (5, 3), (7, 3), (8, 3)]:
		gate = vqubit.UniformSuperpositionGate(m_value, num_qubits)
		assert gate.num_qubits == num_qubits
		assert vqubit.is_unitary(vqubit.unitary(gate))
	with pytest.raises(ValueError):
		vqubit.UniformSuperpositionGate(0, 1)
	with pytest.raises(ValueError):
		vqubit.UniformSuperpositionGate(9, 2)


def test_extra_gate_reprs_round_trip() -> None:
	from vqubit import testing

	for gate in [
		vqubit.CXPowGate(exponent=0.5),
		vqubit.CY,
		vqubit.CCZ,
		vqubit.ISWAP_INV,
		vqubit.SwapPowGate(exponent=0.3),
		vqubit.DiagonalGate([1, 1, 1, -1]),
		vqubit.QubitPermutationGate([0, 2, 1, 3]),
		vqubit.QuantumFourierTransformGate(2),
	]:
		testing.assert_equivalent_repr(gate)
