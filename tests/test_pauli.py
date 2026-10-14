"""Tests for Pauli operators and Pauli strings."""

from __future__ import annotations

import pytest

import vqubit
from tests.helpers import assert_matrix_close

I2 = [[1 + 0j, 0j], [0j, 1 + 0j]]
X2 = [[0j, 1 + 0j], [1 + 0j, 0j]]
Y2 = [[0j, -1j], [1j, 0j]]
Z2 = [[1 + 0j, 0j], [0j, -1 + 0j]]


def test_pauli_matrices() -> None:
	assert_matrix_close(vqubit.Pauli.I.matrix(), I2)
	assert_matrix_close(vqubit.Pauli.X.matrix(), X2)
	assert_matrix_close(vqubit.Pauli.Y.matrix(), Y2)
	assert_matrix_close(vqubit.Pauli.Z.matrix(), Z2)


def test_pauli_on() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.Pauli.Z.on(q0) == vqubit.Z(q0)
	assert vqubit.Pauli.I.on(q0) == vqubit.I(q0)


def test_pauli_string_construction_drops_identity() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	ps = vqubit.PauliString({q0: vqubit.Pauli.X, q1: vqubit.Pauli.I})
	assert ps.qubits == (q0,)
	assert ps[q0] is vqubit.Pauli.X
	assert ps[q1] is vqubit.Pauli.I


def test_pauli_string_str() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	ps = vqubit.PauliString({q0: vqubit.Pauli.X, q1: vqubit.Pauli.Z})
	assert str(ps) == "X(0)*Z(1)"
	assert str(-1j * ps) == "-1j*X(0)*Z(1)"
	assert str(vqubit.PauliString()) == "I"


def test_pauli_string_multiplication() -> None:
	q0 = vqubit.LineQubit(0)
	x = vqubit.PauliString({q0: vqubit.Pauli.X})
	y = vqubit.PauliString({q0: vqubit.Pauli.Y})
	product = x * y
	assert product == vqubit.PauliString({q0: vqubit.Pauli.Z}, coefficient=1j)
	assert x * x == vqubit.PauliString()


def test_pauli_string_scalar_multiplication() -> None:
	q0 = vqubit.LineQubit(0)
	x = vqubit.PauliString({q0: vqubit.Pauli.X})
	assert (2 * x).coefficient == 2
	assert (x * 3j).coefficient == 3j
	assert (-x).coefficient == -1


def test_pauli_string_power() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	ps = vqubit.PauliString({q0: vqubit.Pauli.X, q1: vqubit.Pauli.Z})
	assert ps**2 == vqubit.PauliString()


def test_pauli_string_commutation() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	x0 = vqubit.PauliString({q0: vqubit.Pauli.X})
	z0 = vqubit.PauliString({q0: vqubit.Pauli.Z})
	z1 = vqubit.PauliString({q1: vqubit.Pauli.Z})
	assert not x0.commutes_with(z0)
	assert x0.commutes_with(z1)
	assert x0.commutes_with(x0)


def test_pauli_string_unitary() -> None:
	q0 = vqubit.LineQubit(0)
	ps = vqubit.PauliString({q0: vqubit.Pauli.Z})
	assert_matrix_close(vqubit.unitary(ps), Z2)


def test_pauli_string_validates_types() -> None:
	with pytest.raises(TypeError):
		vqubit.PauliString({vqubit.LineQubit(0): "X"})  # ty: ignore[invalid-argument-type]
	with pytest.raises(TypeError):
		vqubit.PauliString({0: vqubit.Pauli.X})  # ty: ignore[invalid-argument-type]


def test_pauli_string_value_equality() -> None:
	q0 = vqubit.LineQubit(0)
	a = vqubit.PauliString({q0: vqubit.Pauli.X})
	b = vqubit.PauliString({q0: vqubit.Pauli.X})
	assert a == b
	assert hash(a) == hash(b)
	assert a != vqubit.PauliString({q0: vqubit.Pauli.Y})
