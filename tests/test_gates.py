"""Tests for gates: unitaries, powers, inverses and equality."""

from __future__ import annotations

import cmath
import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close, identity_matrix
from vqubit._linalg import dagger, matmul

SQRT_HALF = math.sqrt(0.5)


def test_pauli_x_unitary() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.X), [[0j, 1 + 0j], [1 + 0j, 0j]])


def test_pauli_y_unitary() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.Y), [[0j, -1j], [1j, 0j]])


def test_pauli_z_unitary() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.Z), [[1 + 0j, 0j], [0j, -1 + 0j]])


def test_hadamard_unitary() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.H), [[SQRT_HALF + 0j, SQRT_HALF + 0j], [SQRT_HALF + 0j, -SQRT_HALF + 0j]])


def test_x_squared_is_identity() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.X**2), identity_matrix())


def test_h_squared_is_identity() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.H**2), identity_matrix())


def test_x_sqrt() -> None:
	expected = [[0.5 + 0.5j, 0.5 - 0.5j], [0.5 - 0.5j, 0.5 + 0.5j]]
	assert_matrix_close(vqubit.unitary(vqubit.X**0.5), expected)


def test_s_gate() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.S), [[1 + 0j, 0j], [0j, 1j]])


def test_t_gate() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.T), [[1 + 0j, 0j], [0j, cmath.exp(1j * math.pi / 4)]])


def test_rotations_are_unitary_and_invertible() -> None:
	for gate in [vqubit.Rx(0.3), vqubit.Ry(1.1), vqubit.Rz(-0.7)]:
		matrix = vqubit.unitary(gate)
		product = matmul(matrix, dagger(matrix))
		assert_matrix_close(product, identity_matrix())


def test_rx_matches_expectation() -> None:
	angle = 0.8
	expected = [
		[math.cos(angle / 2) + 0j, -1j * math.sin(angle / 2)],
		[-1j * math.sin(angle / 2), math.cos(angle / 2) + 0j],
	]
	assert_matrix_close(vqubit.unitary(vqubit.Rx(angle)), expected)


def test_inverse_of_power_gate() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.inverse(vqubit.X**0.5)), vqubit.unitary(vqubit.X**-0.5))


def test_inverse_of_rotation() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.inverse(vqubit.Rz(0.4))), vqubit.unitary(vqubit.Rz(-0.4)))


def test_two_qubit_unitaries() -> None:
	assert_matrix_close(
		vqubit.unitary(vqubit.CNOT),
		[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
	)
	assert_matrix_close(
		vqubit.unitary(vqubit.CZ),
		[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]],
	)
	assert_matrix_close(
		vqubit.unitary(vqubit.SWAP),
		[[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]],
	)


def test_gate_equality() -> None:
	assert vqubit.XPowGate() == vqubit.X
	assert vqubit.X != vqubit.Y
	assert vqubit.ZPowGate(exponent=0.5) == vqubit.S
	assert vqubit.CNOT == vqubit.CX


def test_gate_on() -> None:
	q0 = vqubit.LineQubit(0)
	op = vqubit.H.on(q0)
	assert op.gate == vqubit.H
	assert op.qubits == (q0,)


def test_gate_on_wrong_number_of_qubits() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.CNOT.on(q0)


def test_gate_on_each() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	operations = vqubit.H.on_each(q0, q1, q2)
	assert [op.qubits for op in operations] == [(q0,), (q1,), (q2,)]


def test_measurement_gate() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gate = vqubit.MeasurementGate("m", num_qubits=2)
	assert gate.num_qubits == 2
	assert gate.key == "m"
	assert gate.full_invert_mask == (False, False)
	op = gate.on(q0, q1)
	assert op.qubits == (q0, q1)


def test_measurement_invert_mask() -> None:
	gate = vqubit.MeasurementGate("m", num_qubits=3, invert_mask=(True,))
	assert gate.full_invert_mask == (True, False, False)
