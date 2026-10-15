"""Tests for the additional circuit transformers."""

from __future__ import annotations

import vqubit
from tests.helpers import up_to_phase


def test_merge_k_qubit_unitaries() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.H(q0), vqubit.CNOT(q0, q1))
	merged = vqubit.merge_k_qubit_unitaries(circuit, 1)
	operations = list(merged.all_operations())
	assert isinstance(operations[0].gate, vqubit.MatrixGate)
	assert len(operations) == 2
	assert up_to_phase(vqubit.unitary(merged), vqubit.unitary(circuit))


def test_merge_k_qubit_unitaries_keeps_measurements() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	merged = vqubit.merge_k_qubit_unitaries(circuit, 1)
	assert vqubit.is_measurement(merged)


def test_map_operations() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	mapped = vqubit.map_operations(circuit, lambda operation: vqubit.X(operation.qubits[0]))
	assert list(mapped.all_operations()) == [vqubit.X(q0)]


def test_map_operations_to_multiple() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	mapped = vqubit.map_operations(circuit, lambda operation: [vqubit.X(q0), vqubit.Y(q0)])
	assert len(list(mapped.all_operations())) == 2


def test_map_moments() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q0))
	mapped = vqubit.map_moments(circuit, lambda moment: moment)
	assert mapped == circuit


def test_toggle_tags() -> None:
	q0 = vqubit.LineQubit(0)
	tagged = vqubit.toggle_tags(vqubit.Circuit(vqubit.H(q0)), "layer1")
	operation = next(iter(tagged.all_operations()))
	assert isinstance(operation, vqubit.TaggedOperation)
	assert operation.tags == ("layer1",)
	assert "layer1" in str(tagged)


def test_synchronize_terminal_measurements() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="a"), vqubit.H(q1), vqubit.measure(q1, key="b"))
	synchronized = vqubit.synchronize_terminal_measurements(circuit)
	assert len(synchronized.moments[-1]) == 2
	assert vqubit.is_measurement(synchronized)


def test_synchronize_leaves_non_terminal_alone() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="m"), vqubit.X(q0))
	assert vqubit.synchronize_terminal_measurements(circuit) == circuit
