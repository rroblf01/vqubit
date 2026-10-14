"""Tests for MatrixGate, pauli_expansion, phase_by and apply_unitary."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close

SQRT_HALF = math.sqrt(0.5)


def test_matrix_gate_wraps_a_matrix() -> None:
	gate = vqubit.MatrixGate(vqubit.unitary(vqubit.X))
	assert gate.num_qubits == 1
	assert_matrix_close(vqubit.unitary(gate), vqubit.unitary(vqubit.X))


def test_matrix_gate_equality() -> None:
	assert vqubit.MatrixGate(vqubit.unitary(vqubit.X)) == vqubit.MatrixGate(vqubit.unitary(vqubit.X))
	assert vqubit.MatrixGate(vqubit.unitary(vqubit.X)) != vqubit.MatrixGate(vqubit.unitary(vqubit.Z))


def test_matrix_gate_validates_shape() -> None:
	with pytest.raises(ValueError):
		vqubit.MatrixGate([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # 3x3 is not a power of two
	with pytest.raises(ValueError):
		vqubit.MatrixGate([[1, 0, 0], [0, 1, 0]])  # not square


def test_pauli_expansion_of_x() -> None:
	q0 = vqubit.LineQubit(0)
	expansion = vqubit.pauli_expansion(vqubit.X(q0))
	assert len(expansion) == 1
	assert expansion[vqubit.PauliString({q0: vqubit.Pauli.X})] == pytest.approx(1 + 0j)


def test_pauli_expansion_of_hadamard() -> None:
	q0 = vqubit.LineQubit(0)
	expansion = vqubit.pauli_expansion(vqubit.H(q0))
	assert expansion[vqubit.PauliString({q0: vqubit.Pauli.X})] == pytest.approx(SQRT_HALF)
	assert expansion[vqubit.PauliString({q0: vqubit.Pauli.Z})] == pytest.approx(SQRT_HALF)


def test_pauli_expansion_of_cnot() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	expansion = vqubit.pauli_expansion(vqubit.CNOT(q0, q1))
	# CNOT = (II + ZI + IX - ZX) / 2.
	assert expansion[vqubit.PauliString()] == pytest.approx(0.5)
	assert expansion[vqubit.PauliString({q0: vqubit.Pauli.Z})] == pytest.approx(0.5)
	assert expansion[vqubit.PauliString({q1: vqubit.Pauli.X})] == pytest.approx(0.5)
	assert expansion[vqubit.PauliString({q0: vqubit.Pauli.Z, q1: vqubit.Pauli.X})] == pytest.approx(-0.5)


def test_pauli_sum_str() -> None:
	q0 = vqubit.LineQubit(0)
	terms = {vqubit.PauliString({q0: vqubit.Pauli.X}): 1 + 0j}
	expansion = vqubit.PauliSum(terms)
	assert str(expansion) == "1*X(0)"
	assert vqubit.PauliSum()[vqubit.PauliString({q0: vqubit.Pauli.X})] == 0j
	assert len(vqubit.PauliSum(terms)) == 1


def test_phase_by_maps_x_to_y() -> None:
	# Z**0.5 X Z**-0.5 = Y.
	q0 = vqubit.LineQubit(0)
	phased = vqubit.phase_by(vqubit.X(q0), 0.5, 0)
	assert_matrix_close(vqubit.unitary(phased), vqubit.unitary(vqubit.Y))


def test_phase_by_validates_index() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.phase_by(vqubit.X(q0), 0.5, 3)


def test_apply_unitary() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.apply_unitary(vqubit.X(q0), [1 + 0j, 0j]) == [0j, 1 + 0j]
	hadamard = vqubit.apply_unitary(vqubit.H(q0), [1 + 0j, 0j])
	assert hadamard[0] == pytest.approx(SQRT_HALF)
	assert hadamard[1] == pytest.approx(SQRT_HALF)


def test_apply_unitary_validates_length() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.apply_unitary(vqubit.X(q0), [1 + 0j, 0j, 0j, 0j])
