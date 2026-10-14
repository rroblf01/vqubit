"""Tests for controlled gates and the Toffoli/Fredkin gates."""

from __future__ import annotations

import pytest

import vqubit
from tests.helpers import assert_matrix_close

I4 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]


def test_controlled_x_is_cnot() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	controlled = vqubit.X.controlled_by(q0)
	assert_matrix_close(vqubit.unitary(controlled.on(q0, q1)), vqubit.unitary(vqubit.CNOT))


def test_controlled_z_is_cz() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	controlled = vqubit.Z.controlled_by(q0)
	assert_matrix_close(vqubit.unitary(controlled.on(q0, q1)), vqubit.unitary(vqubit.CZ))


def test_control_values_zero() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	# Controls on a zero: X on q1 only when q0 is 0.
	gate = vqubit.X.controlled_by(q1, control_values=[0])
	matrix = vqubit.unitary(gate.on(q0, q1))
	expected = [
		[0j, 1 + 0j, 0j, 0j],
		[1 + 0j, 0j, 0j, 0j],
		[0j, 0j, 1 + 0j, 0j],
		[0j, 0j, 0j, 1 + 0j],
	]
	assert_matrix_close(matrix, expected)


def test_controlled_gate_num_qubits() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.X.controlled_by(q0).num_qubits == 2
	assert vqubit.X.controlled_by(q0, q1).num_qubits == 3


def test_controlled_gate_validates_control_values() -> None:
	with pytest.raises(ValueError):
		vqubit.ControlledGate(vqubit.X, num_controls=2, control_values=[True])
	with pytest.raises(ValueError):
		vqubit.ControlledGate(vqubit.X, num_controls=0)


def test_controlled_gate_inverse() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gate = vqubit.S.controlled_by(q0)
	inverted = vqubit.inverse(gate)
	assert_matrix_close(
		vqubit.unitary(inverted.on(q0, q1)), vqubit.unitary((vqubit.S**-1).controlled_by(q0).on(q0, q1))
	)


def test_controlled_gate_equality() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.X.controlled_by(q0) == vqubit.X.controlled_by(q0)
	assert vqubit.X.controlled_by(q0) != vqubit.Z.controlled_by(q0)


def test_ccx_unitary_is_a_basis_swap() -> None:
	expected = [row[:] for row in vqubit.unitary(vqubit.CCX)]
	# Only the block of the two set controls differs from the identity.
	for i in range(8):
		for j in range(8):
			if (i, j) in {(6, 6), (7, 7)}:
				assert expected[i][j] == 0j
			elif (i, j) == (6, 7) or (i, j) == (7, 6):
				assert expected[i][j] == 1 + 0j
			else:
				assert expected[i][j] == (1 + 0j if i == j else 0j)


def test_ccx_flips_the_target() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.X(q1), vqubit.CCX(q0, q1, q2))
	result = vqubit.Simulator().simulate(circuit)
	assert result.state_vector()[0b111] == pytest.approx(1 + 0j)


def test_cswap_unitary_is_a_basis_swap() -> None:
	matrix = vqubit.unitary(vqubit.CSWAP)
	for i in range(8):
		for j in range(8):
			if (i, j) in {(5, 5), (6, 6)}:
				assert matrix[i][j] == 0j
			elif (i, j) == (5, 6) or (i, j) == (6, 5):
				assert matrix[i][j] == 1 + 0j
			else:
				assert matrix[i][j] == (1 + 0j if i == j else 0j)


def test_cswap_swaps_when_control_is_set() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	# Prepare q0=1, q1=1: |110>. CSWAP swaps q1 and q2 -> |101>.
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.X(q1), vqubit.CSWAP(q0, q1, q2))
	result = vqubit.Simulator().simulate(circuit)
	assert result.state_vector()[0b101] == pytest.approx(1 + 0j)


def test_ccx_does_nothing_without_both_controls() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.CCX(q0, q1, q2))
	result = vqubit.Simulator().simulate(circuit)
	assert result.state_vector()[0b100] == pytest.approx(1 + 0j)
