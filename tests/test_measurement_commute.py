"""Tests for the measurement and commutation protocols."""

from __future__ import annotations

import vqubit


def test_is_measurement() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.is_measurement(vqubit.measure(q0))
	assert vqubit.is_measurement(vqubit.MeasurementGate("m"))
	assert vqubit.is_measurement(vqubit.Circuit(vqubit.measure(q0)))
	assert not vqubit.is_measurement(vqubit.H)
	assert not vqubit.is_measurement(vqubit.H(q0))
	assert not vqubit.is_measurement(vqubit.Circuit(vqubit.H(q0)))


def test_measurement_key_names_of_operation() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.measurement_key_names(vqubit.measure(q0, key="alpha")) == {"alpha"}
	assert vqubit.measurement_key_names(vqubit.H(q0)) == set()


def test_measurement_key_names_of_moment_and_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	moment = vqubit.Moment([vqubit.measure(q0, key="a"), vqubit.measure(q1, key="b")])
	assert vqubit.measurement_key_names(moment) == {"a", "b"}
	circuit = vqubit.Circuit(vqubit.measure(q0, key="a"), vqubit.X(q1), vqubit.measure(q1, key="c"))
	assert vqubit.measurement_key_names(circuit) == {"a", "c"}


def test_measurement_keys_touched_is_an_alias() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="k"))
	assert vqubit.measurement_keys_touched(circuit) == vqubit.measurement_key_names(circuit)


def test_operations_on_disjoint_qubits_commute() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.commutes(vqubit.X(q0), vqubit.Z(q1))
	assert vqubit.definitely_commutes(vqubit.X(q0), vqubit.Z(q1))


def test_operations_on_the_same_qubit() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.commutes(vqubit.X(q0), vqubit.X(q0))
	assert not vqubit.commutes(vqubit.X(q0), vqubit.Z(q0))
	assert vqubit.commutes(vqubit.Z(q0), vqubit.S(q0))


def test_gates_commute() -> None:
	assert vqubit.commutes(vqubit.CNOT, vqubit.CNOT)
	assert not vqubit.commutes(vqubit.H, vqubit.X)


def test_findall_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"), vqubit.X(q1))
	found = list(circuit.findall_operations_with_gate_type(vqubit.MeasurementGate))
	assert len(found) == 1
	index, operation = found[0]
	assert index == 1
	gate = operation.gate
	assert isinstance(gate, vqubit.MeasurementGate)
	assert gate.key == "m"
	assert list(circuit.findall_operations(lambda op: op.qubits == (q1,))) == [(0, vqubit.X(q1))]
