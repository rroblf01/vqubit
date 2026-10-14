"""Tests for compiling circuits to a target gate set."""

from __future__ import annotations

import pytest

import vqubit


def _up_to_phase(actual: list[list[complex]], expected: list[list[complex]]) -> bool:
	phase = None
	for j in range(len(expected[0])):
		if abs(expected[0][j]) > 1e-9:
			phase = actual[0][j] / expected[0][j]
			break
	if phase is None:
		return False
	return all(
		abs(actual[i][j] - phase * expected[i][j]) < 1e-9 for i in range(len(expected)) for j in range(len(expected))
	)


def test_cz_target_gateset_compiles() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.Rx(0.3)(q0), vqubit.H(q1), vqubit.CNOT(q0, q1), vqubit.CCX(q0, q1, q2))
	compiled = vqubit.optimize_for_target_gateset(circuit, vqubit.CZTargetGateset())
	operations = list(compiled.all_operations())
	assert all(isinstance(op.gate, (vqubit.CZGate, vqubit.PhasedXZGate)) for op in operations)
	assert any(isinstance(op.gate, vqubit.CZGate) for op in operations)
	assert _up_to_phase(vqubit.unitary(compiled), vqubit.unitary(circuit))


def test_cz_target_converts_cnot_to_cz() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	compiled = vqubit.optimize_for_target_gateset(vqubit.Circuit(vqubit.CNOT(q0, q1)), vqubit.CZTargetGateset())
	two_qubit = [op for op in compiled.all_operations() if len(op.qubits) == 2]
	assert len(two_qubit) == 1
	assert isinstance(two_qubit[0].gate, vqubit.CZGate)
	assert _up_to_phase(vqubit.unitary(compiled), vqubit.unitary(vqubit.CNOT))


def test_cnot_target_converts_cz_to_cnot() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	compiled = vqubit.optimize_for_target_gateset(vqubit.Circuit(vqubit.CZ(q0, q1)), vqubit.CNOTTargetGateset())
	two_qubit = [op for op in compiled.all_operations() if len(op.qubits) == 2]
	assert len(two_qubit) == 1
	assert isinstance(two_qubit[0].gate, vqubit.CNOTGate)
	assert _up_to_phase(vqubit.unitary(compiled), vqubit.unitary(vqubit.CZ))


def test_measurements_are_preserved() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	compiled = vqubit.optimize_for_target_gateset(circuit, vqubit.CZTargetGateset())
	assert vqubit.is_measurement(compiled)
	assert vqubit.measurement_key_names(compiled) == {"m"}


def test_unsupported_two_qubit_gate_raises() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	with pytest.raises(TypeError):
		vqubit.optimize_for_target_gateset(vqubit.Circuit(vqubit.ISWAP(q0, q1)), vqubit.CZTargetGateset())


def test_default_target_is_cz() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	compiled = vqubit.optimize_for_target_gateset(vqubit.Circuit(vqubit.CNOT(q0, q1)))
	assert any(isinstance(op.gate, vqubit.CZGate) for op in compiled.all_operations())
