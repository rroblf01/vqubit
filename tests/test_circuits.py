"""Tests for circuits and moment packing."""

from __future__ import annotations

import pytest

import vqubit


def test_circuit_packs_independent_gates_into_one_moment() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q1))
	assert len(circuit) == 1


def test_circuit_separates_dependent_gates() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	assert len(circuit) == 2
	assert circuit.moments[0].qubits == (q0,)
	assert circuit.moments[1].qubits == (q0, q1)


def test_circuit_all_qubits() -> None:
	q0, q2 = vqubit.LineQubit(0), vqubit.LineQubit(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q2))
	assert circuit.all_qubits() == (q0, q2)


def test_circuit_all_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	assert list(circuit.all_operations()) == [vqubit.H(q0), vqubit.CNOT(q0, q1)]


def test_circuit_iteration_and_indexing() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q0))
	assert len(circuit) == 2
	assert circuit[0] == vqubit.Moment([vqubit.H(q0)])
	assert circuit[1] == vqubit.Moment([vqubit.X(q0)])
	assert circuit[-1] == vqubit.Moment([vqubit.X(q0)])
	assert circuit[0:1] == vqubit.Circuit.from_moments([vqubit.Moment([vqubit.H(q0)])])


def test_circuit_equality() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	a = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	b = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	assert a == b
	assert a != vqubit.Circuit(vqubit.H(q0))


def test_circuit_append_list() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit()
	circuit.append([vqubit.H(q0), vqubit.X(q1)])
	assert len(circuit) == 1


def test_circuit_insert() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q0))
	circuit.insert(0, vqubit.Z(q0))
	# Z depends on the qubit, so it lands after the existing operations.
	assert circuit[0] == vqubit.Moment([vqubit.H(q0)])
	assert circuit[2] == vqubit.Moment([vqubit.Z(q0)])
	assert len(circuit) == 3


def test_circuit_new_strategy_forces_a_new_moment() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), strategy=vqubit.InsertStrategy.NEW)
	circuit.append(vqubit.X(q1), strategy=vqubit.InsertStrategy.NEW)
	assert len(circuit) == 2


def test_circuit_rejects_bad_content() -> None:
	with pytest.raises(TypeError):
		vqubit.Circuit("nope")  # type: ignore[arg-type]


def test_frozen_circuit_is_immutable() -> None:
	q0 = vqubit.LineQubit(0)
	frozen = vqubit.FrozenCircuit(vqubit.H(q0))
	with pytest.raises(TypeError):
		frozen.append(vqubit.X(q0))
	with pytest.raises(TypeError):
		frozen.insert(0, vqubit.X(q0))
