"""Tests for the circuit transformers."""

from __future__ import annotations

import vqubit
from tests.helpers import assert_matrix_close


def test_drop_empty_moments() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit.from_moments([vqubit.Moment(), vqubit.Moment([vqubit.H(q0)]), vqubit.Moment()])
	dropped = vqubit.drop_empty_moments(circuit)
	assert len(dropped) == 1
	assert dropped[0] == vqubit.Moment([vqubit.H(q0)])


def test_drop_negligible_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.I(q0), vqubit.X(q1), vqubit.CNOT(q0, q1))
	dropped = vqubit.drop_negligible_operations(circuit)
	assert vqubit.I(q0) not in list(dropped.all_operations())
	assert len(list(dropped.all_operations())) == 3


def test_drop_negligible_operations_keeps_measurements() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="m"))
	assert len(list(vqubit.drop_negligible_operations(circuit).all_operations())) == 1


def test_expand_composite() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.CCX(q0, q1, q2))
	expanded = vqubit.expand_composite(circuit)
	assert all(not isinstance(op.gate, vqubit.CCXGate) for op in expanded.all_operations())
	assert_matrix_close(vqubit.unitary(expanded), vqubit.unitary(vqubit.CCX))


def test_stratified_circuit_orders_by_category() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.X(q0))
	one_qubit: list[vqubit.Operation] = []
	two_qubit: list[vqubit.Operation] = []
	for op in vqubit.stratified_circuit(
		circuit,
		lambda op: len(op.qubits) == 2,
		lambda op: len(op.qubits) == 1,
	).all_operations():
		(one_qubit if len(op.qubits) == 1 else two_qubit).append(op)
	assert two_qubit == [vqubit.CNOT(q0, q1)]
	assert one_qubit == [vqubit.H(q0), vqubit.X(q0)]


def test_stratified_circuit_accepts_gate_types() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.X(q0))
	stratified = vqubit.stratified_circuit(circuit, [vqubit.CNOTGate, vqubit.XPowGate])
	assert next(iter(stratified.all_operations())) == vqubit.CNOT(q0, q1)


def test_align_left_keeps_moments() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	assert vqubit.align_left(circuit) == circuit


def test_align_right_preserves_the_unitary() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q1), vqubit.CNOT(q0, q1))
	right = vqubit.align_right(circuit)
	assert set(right.all_operations()) == set(circuit.all_operations())
	assert_matrix_close(vqubit.unitary(right), vqubit.unitary(circuit))
