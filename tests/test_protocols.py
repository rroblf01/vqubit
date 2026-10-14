"""Tests for the protocol functions."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close

SQRT_HALF = math.sqrt(0.5)


def test_unitary_of_gate() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.X), [[0j, 1 + 0j], [1 + 0j, 0j]])


def test_unitary_of_operation_uses_gate_matrix() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert_matrix_close(vqubit.unitary(vqubit.CNOT.on(q0, q1)), vqubit.unitary(vqubit.CNOT))


def test_unitary_of_bell_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	s = SQRT_HALF
	expected = [
		[s + 0j, 0j, s + 0j, 0j],
		[0j, s + 0j, 0j, s + 0j],
		[0j, s + 0j, 0j, -s + 0j],
		[s + 0j, 0j, -s + 0j, 0j],
	]
	assert_matrix_close(vqubit.unitary(circuit), expected)


def test_has_unitary() -> None:
	assert vqubit.has_unitary(vqubit.X)
	assert not vqubit.has_unitary(vqubit.measure(vqubit.LineQubit(0)))


def test_unitary_of_measurement_raises() -> None:
	with pytest.raises(TypeError):
		vqubit.unitary(vqubit.measure(vqubit.LineQubit(0)))


def test_unitary_of_rejects_unknown() -> None:
	with pytest.raises(TypeError):
		vqubit.unitary(42)


def test_num_qubits_and_shape() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.num_qubits(vqubit.CNOT) == 2
	assert vqubit.num_qubits(vqubit.CNOT(q0, q1)) == 2
	assert vqubit.num_qubits(vqubit.Circuit(vqubit.CNOT(q0, q1))) == 2
	assert vqubit.qid_shape(vqubit.CNOT) == (2, 2)


def test_inverse_of_operation() -> None:
	q0 = vqubit.LineQubit(0)
	op = vqubit.S.on(q0)
	inverted = vqubit.inverse(op)
	assert inverted.gate == vqubit.S**-1
	assert inverted.qubits == (q0,)


def test_inverse_of_circuit_reverses_order() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	inverted = vqubit.inverse(circuit)
	operations = list(inverted.all_operations())
	assert operations[0] == vqubit.CNOT(q0, q1)
	assert operations[1].gate == vqubit.inverse(vqubit.H)
	assert operations[1].qubits == (q0,)


def test_pow_of_gate() -> None:
	assert vqubit.pow(vqubit.X, 0.5) == vqubit.X**0.5


def test_decompose_of_circuit() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	assert vqubit.decompose(circuit) == [vqubit.H(q0)]
