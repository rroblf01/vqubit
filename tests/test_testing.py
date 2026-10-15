"""Tests for the vqubit.testing helpers."""

from __future__ import annotations

import pytest

import vqubit
from vqubit import testing


def test_assert_equivalent_repr() -> None:
	testing.assert_equivalent_repr(vqubit.LineQubit(3))
	testing.assert_equivalent_repr(vqubit.X)
	testing.assert_equivalent_repr(vqubit.Rx(0.3))
	testing.assert_equivalent_repr(vqubit.CNOT)
	testing.assert_equivalent_repr(vqubit.FSimGate(1.0, 0.5))


def test_assert_same_circuits() -> None:
	q0 = vqubit.LineQubit(0)
	testing.assert_same_circuits(vqubit.Circuit(vqubit.H(q0)), vqubit.Circuit(vqubit.H(q0)))
	with pytest.raises(AssertionError):
		testing.assert_same_circuits(vqubit.Circuit(vqubit.H(q0)), vqubit.Circuit(vqubit.X(q0)))


def test_assert_circuits_with_terminal_measurements_are_equivalent() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	a = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="m"))
	b = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	testing.assert_circuits_with_terminal_measurements_are_equivalent(a, b)
	with pytest.raises(AssertionError):
		testing.assert_circuits_with_terminal_measurements_are_equivalent(a, vqubit.Circuit(vqubit.X(q0)))


def test_assert_allclose_up_to_global_phase() -> None:
	testing.assert_allclose_up_to_global_phase(vqubit.unitary(vqubit.H), vqubit.unitary(vqubit.H))
	with pytest.raises(AssertionError):
		testing.assert_allclose_up_to_global_phase(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.Z))


def test_random_circuit_is_deterministic() -> None:
	qubits = vqubit.LineQubit.range(3)
	a = testing.random_circuit(qubits, 5, 0.5, random_state=0)
	b = testing.random_circuit(qubits, 5, 0.5, random_state=0)
	assert a == b
	assert len(list(a.all_operations())) > 0


def test_random_circuit_respects_gate_domain() -> None:
	qubits = vqubit.LineQubit.range(2)
	circuit = testing.random_circuit(qubits, 3, 1.0, random_state=1, gate_domain={1: [vqubit.H]})
	assert all(isinstance(op.gate, vqubit.HPowGate) for op in circuit.all_operations())
