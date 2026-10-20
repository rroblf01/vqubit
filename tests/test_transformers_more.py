"""Tests for the additional circuit transformers."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import pytest

import vqubit
from tests.helpers import up_to_phase


def _tags(operation: vqubit.Operation) -> tuple[Any, ...]:
	assert isinstance(operation, vqubit.TaggedOperation)
	return operation.tags


def _merge_gates(op1: vqubit.Operation, op2: vqubit.Operation) -> vqubit.Operation:
	"""Compose two operations that act on the same qubits."""
	matrix = vqubit.dot(vqubit.unitary(op2.gate), vqubit.unitary(op1.gate))
	return vqubit.MatrixGate(matrix).on(*op1.qubits)


def _unrolled_unitary(circuit: vqubit.Circuit) -> Sequence[Sequence[complex]]:
	return vqubit.unitary(vqubit.unroll_circuit_op(circuit, deep=True))


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


def test_remove_tags_removes_a_subset_and_unwraps_the_operation() -> None:
	q0 = vqubit.LineQubit(0)
	tagged = vqubit.toggle_tags(vqubit.Circuit(vqubit.H(q0)), "a", "b")
	one_left = vqubit.remove_tags(tagged, ["a"])
	operation = next(iter(one_left.all_operations()))
	assert isinstance(operation, vqubit.TaggedOperation)
	assert operation.tags == ("b",)
	untagged = vqubit.remove_tags(tagged, ["a", "b"])
	operation = next(iter(untagged.all_operations()))
	assert not isinstance(operation, vqubit.TaggedOperation)
	assert operation == vqubit.H(q0)
	assert up_to_phase(vqubit.unitary(untagged), vqubit.unitary(tagged))


def test_remove_tags_leaves_unknown_tags_alone() -> None:
	q0 = vqubit.LineQubit(0)
	tagged = vqubit.toggle_tags(vqubit.Circuit(vqubit.H(q0)), "a", "b")
	unchanged = vqubit.remove_tags(tagged, ["missing"])
	operation = next(iter(unchanged.all_operations()))
	assert isinstance(operation, vqubit.TaggedOperation)
	assert operation.tags == ("a", "b")


def test_remove_tags_keeps_untagged_operations() -> None:
	q0 = vqubit.LineQubit(0)
	tagged = next(iter(vqubit.toggle_tags(vqubit.Circuit(vqubit.X(q0)), "keep").all_operations()))
	circuit = vqubit.Circuit.from_moments([vqubit.Moment([vqubit.H(q0)]), vqubit.Moment([tagged])])
	result = vqubit.remove_tags(circuit, ["keep"])
	operations = list(result.all_operations())
	assert len(operations) == 2
	assert not any(isinstance(operation, vqubit.TaggedOperation) for operation in operations)
	assert operations == [vqubit.H(q0), vqubit.X(q0)]


def test_index_tags_numbers_every_operation_in_order() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	indexed = vqubit.index_tags(circuit, start=5, prefix="t")
	assert [_tags(operation) for operation in indexed.all_operations()] == [("t5",), ("t6",)]
	assert up_to_phase(vqubit.unitary(indexed), vqubit.unitary(circuit))


def test_is_negligible_turn_uses_the_integer_tolerance() -> None:
	assert vqubit.is_negligible_turn(0.0)
	assert vqubit.is_negligible_turn(2.0)
	assert vqubit.is_negligible_turn(-3.000001, tolerance=1e-3)
	assert not vqubit.is_negligible_turn(2.5)
	assert not vqubit.is_negligible_turn(1.2, tolerance=0.1)


def test_unroll_circuit_op_expands_repetitions_and_nested_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	inner = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	wrapped = vqubit.CircuitOperation(inner, repetitions=2)
	unrolled = vqubit.unroll_circuit_op(vqubit.Circuit(wrapped))
	assert len(list(unrolled.all_operations())) == 4
	assert up_to_phase(vqubit.unitary(unrolled), vqubit.unitary(wrapped))
	outer = vqubit.CircuitOperation(vqubit.Circuit(wrapped), repetitions=2)
	shallow = vqubit.unroll_circuit_op(vqubit.Circuit(outer))
	assert len(list(shallow.all_operations())) == 2
	assert all(isinstance(operation, vqubit.CircuitOperation) for operation in shallow.all_operations())
	deep = vqubit.unroll_circuit_op(vqubit.Circuit(outer), deep=True)
	assert len(list(deep.all_operations())) == 8
	expected = vqubit.dot(vqubit.unitary(wrapped), vqubit.unitary(wrapped))
	assert up_to_phase(vqubit.unitary(deep), expected)


def test_unroll_circuit_op_greedy_variants_match_the_plain_unroll() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.CircuitOperation(vqubit.Circuit(vqubit.H(q0)), repetitions=2))
	expected = list(vqubit.unroll_circuit_op(circuit).all_operations())
	assert list(vqubit.unroll_circuit_op_greedy_earliest(circuit).all_operations()) == expected
	assert list(vqubit.unroll_circuit_op_greedy_frontier(circuit).all_operations()) == expected


def test_map_operations_and_unroll_accepts_operations_and_iterables() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	mapped = vqubit.map_operations_and_unroll(circuit, lambda operation: vqubit.X(operation.qubits[0]))
	assert list(mapped.all_operations()) == [vqubit.X(q0)]
	nested = vqubit.Circuit(vqubit.CircuitOperation(vqubit.Circuit(vqubit.H(q0)), repetitions=2))
	mapped = vqubit.map_operations_and_unroll(
		nested, lambda operation: [vqubit.X(operation.qubits[0]), vqubit.Y(operation.qubits[0])]
	)
	assert len(list(mapped.all_operations())) == 2
	assert len(list(vqubit.map_operations_and_unroll(nested, lambda operation: operation).all_operations())) == 2


def test_drop_terminal_measurements() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	terminal = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	assert list(vqubit.drop_terminal_measurements(terminal).all_operations()) == [vqubit.H(q0)]
	non_terminal = vqubit.Circuit(vqubit.measure(q0, key="m"), vqubit.X(q0))
	assert vqubit.is_measurement(vqubit.drop_terminal_measurements(non_terminal))
	partial = vqubit.Circuit(vqubit.measure(q0, q1, key="m"), vqubit.H(q1))
	assert len(list(vqubit.drop_terminal_measurements(partial).all_operations())) == 2
	repeated = vqubit.Circuit(vqubit.measure(q0, key="a"), vqubit.measure(q0, key="b"))
	assert list(vqubit.drop_terminal_measurements(repeated).all_operations()) == []


def test_drop_diagonal_before_measurement() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	measurement = vqubit.measure(q0, key="m")
	assert (
		len(list(vqubit.drop_diagonal_before_measurement(vqubit.Circuit(vqubit.Z(q0), measurement)).all_operations()))
		== 1
	)
	assert (
		len(list(vqubit.drop_diagonal_before_measurement(vqubit.Circuit(vqubit.X(q0), measurement)).all_operations()))
		== 2
	)
	assert (
		len(list(vqubit.drop_diagonal_before_measurement(vqubit.Circuit(vqubit.Z(q1), measurement)).all_operations()))
		== 2
	)
	separated = vqubit.Circuit(vqubit.Z(q0), vqubit.H(q0), measurement)
	assert len(list(vqubit.drop_diagonal_before_measurement(separated).all_operations())) == 3
	multi = vqubit.measure(q0, q1, key="m")
	assert (
		len(list(vqubit.drop_diagonal_before_measurement(vqubit.Circuit(vqubit.CZ(q0, q1), multi)).all_operations()))
		== 1
	)
	measurements = vqubit.Circuit(vqubit.measure(q0, key="a"), vqubit.measure(q0, key="b"))
	assert len(list(vqubit.drop_diagonal_before_measurement(measurements).all_operations())) == 2


def test_dephase_measurements_replaces_measure_reset_with_a_dephasing_channel() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"), vqubit.reset(q0))
	result = vqubit.dephase_measurements(circuit)
	operations = list(result.all_operations())
	assert len(operations) == 2
	assert isinstance(operations[0].gate, vqubit.HPowGate)
	assert isinstance(operations[1].gate, vqubit.PhaseDampingChannel)
	assert not vqubit.is_measurement(result)
	reset_only = vqubit.Circuit(vqubit.reset(q0))
	assert list(vqubit.dephase_measurements(reset_only).all_operations()) == list(reset_only.all_operations())


def test_dephase_measurements_keeps_a_plain_measurement() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	assert vqubit.is_measurement(vqubit.dephase_measurements(circuit))


def test_merge_operations_merges_adjacent_gates_and_keeps_the_unitary() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q0), vqubit.Y(q0))
	merged = vqubit.merge_operations(circuit, _merge_gates)
	operations = list(merged.all_operations())
	assert len(operations) == 1
	assert isinstance(operations[0].gate, vqubit.MatrixGate)
	assert up_to_phase(vqubit.unitary(merged), vqubit.unitary(circuit))
	assert len(list(vqubit.merge_operations(circuit, lambda left, right: None).all_operations())) == 3


def test_merge_operations_merges_two_qubit_gates_and_respects_tags() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.CNOT(q0, q1), vqubit.CNOT(q0, q1))
	merged = vqubit.merge_operations(circuit, _merge_gates)
	assert len(list(merged.all_operations())) == 1
	assert up_to_phase(vqubit.unitary(merged), vqubit.unitary(circuit))
	tagged = next(iter(vqubit.toggle_tags(vqubit.Circuit(vqubit.X(q0)), "keep").all_operations()))
	tagged_circuit = vqubit.Circuit.from_moments([vqubit.Moment([vqubit.H(q0)]), vqubit.Moment([tagged])])
	protected = vqubit.merge_operations(tagged_circuit, _merge_gates, tags_to_ignore=["keep"])
	assert len(list(protected.all_operations())) == 2


def test_merge_operations_protects_a_tagged_left_hand_operation() -> None:
	q0 = vqubit.LineQubit(0)
	tagged = next(iter(vqubit.toggle_tags(vqubit.Circuit(vqubit.H(q0)), "keep").all_operations()))
	circuit = vqubit.Circuit.from_moments([vqubit.Moment([tagged]), vqubit.Moment([vqubit.X(q0)])])
	merged = vqubit.merge_operations(circuit, _merge_gates, tags_to_ignore=["keep"])
	assert len(list(merged.all_operations())) == 2


def test_merge_operations_rejects_a_merged_operation_on_new_qubits() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	with pytest.raises(ValueError, match="subset"):
		vqubit.merge_operations(vqubit.Circuit(vqubit.X(q0), vqubit.X(q0)), lambda left, right: vqubit.CNOT(q0, q1))


def test_merge_operations_handles_operations_on_subsets_of_qubits() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.CNOT(q0, q1))
	merged = vqubit.merge_operations(
		circuit, lambda left, right: right if set(left.qubits).issubset(right.qubits) else None
	)
	assert list(merged.all_operations()) == [vqubit.CNOT(q0, q1)]


def test_merge_moments_uses_the_merge_function() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit.from_moments([vqubit.Moment([vqubit.H(q0)]), vqubit.Moment([vqubit.X(q0)])])
	merged = vqubit.merge_moments(circuit, lambda first, second: vqubit.Moment([vqubit.Y(q0)]))
	assert list(merged.all_operations()) == [vqubit.Y(q0)]
	kept = vqubit.merge_moments(circuit, lambda first, second: None)
	assert kept == circuit
	assert up_to_phase(vqubit.unitary(kept), vqubit.unitary(circuit))
	assert vqubit.merge_moments_batch(circuit, lambda first, second: vqubit.Moment([vqubit.Y(q0)])) == merged


def test_merge_moments_respects_tags_to_ignore() -> None:
	q0 = vqubit.LineQubit(0)
	tagged = next(iter(vqubit.toggle_tags(vqubit.Circuit(vqubit.X(q0)), "keep").all_operations()))
	circuit = vqubit.Circuit.from_moments([vqubit.Moment([vqubit.H(q0)]), vqubit.Moment([tagged])])
	merged = vqubit.merge_moments(circuit, lambda first, second: vqubit.Moment([vqubit.Y(q0)]), tags_to_ignore=["keep"])
	assert len(list(merged.all_operations())) == 2
	assert isinstance(list(merged.all_operations())[1], vqubit.TaggedOperation)
	unprotected = vqubit.merge_moments(circuit, lambda first, second: vqubit.Moment([vqubit.Y(q0)]))
	assert list(unprotected.all_operations()) == [vqubit.Y(q0)]


def test_merge_operations_to_circuit_op_groups_runs() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q0))
	merged = vqubit.merge_operations_to_circuit_op(circuit)
	operations = list(merged.all_operations())
	assert len(operations) == 1
	assert isinstance(operations[0], vqubit.CircuitOperation)
	assert up_to_phase(_unrolled_unitary(merged), vqubit.unitary(circuit))
	grouped = vqubit.Circuit.from_moments(
		[vqubit.Moment([vqubit.H(q0)]), vqubit.Moment([vqubit.H(q1)]), vqubit.Moment([vqubit.X(q0)])]
	)
	grouped_result = vqubit.merge_operations_to_circuit_op(grouped)
	assert len(list(grouped_result.all_operations())) == 3
	assert up_to_phase(_unrolled_unitary(grouped_result), vqubit.unitary(grouped))
	assert list(vqubit.merge_operations_to_circuit_op(vqubit.Circuit()).all_operations()) == []


def test_merge_k_qubit_unitaries_to_circuit_op_delegates() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q0), vqubit.CNOT(q0, q1))
	result = vqubit.merge_k_qubit_unitaries_to_circuit_op(circuit, 1)
	expected = vqubit.merge_k_qubit_unitaries(circuit, 1)
	assert list(result.all_operations()) == list(expected.all_operations())
	assert up_to_phase(vqubit.unitary(result), vqubit.unitary(circuit))


def test_merge_single_qubit_gates_to_phxz() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.T(q0), vqubit.CNOT(q0, q1), vqubit.H(q1))
	merged = vqubit.merge_single_qubit_gates_to_phxz(circuit)
	operations = list(merged.all_operations())
	assert len(operations) == 3
	assert isinstance(operations[0].gate, vqubit.PhasedXZGate)
	assert len(operations[1].qubits) == 2
	assert isinstance(operations[2].gate, vqubit.PhasedXZGate)
	assert up_to_phase(vqubit.unitary(merged), vqubit.unitary(circuit))


def test_merge_single_qubit_gates_to_phased_x_and_z() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.T(q0), vqubit.CNOT(q0, q1), vqubit.H(q1))
	merged = vqubit.merge_single_qubit_gates_to_phased_x_and_z(circuit)
	operations = list(merged.all_operations())
	assert len(operations) == 5
	assert isinstance(operations[0].gate, vqubit.PhasedXPowGate)
	assert isinstance(operations[1].gate, vqubit.ZPowGate)
	assert len(operations[2].qubits) == 2
	assert isinstance(operations[3].gate, vqubit.PhasedXPowGate)
	assert isinstance(operations[4].gate, vqubit.ZPowGate)
	assert up_to_phase(vqubit.unitary(merged), vqubit.unitary(circuit))


def test_merge_single_qubit_gates_flushes_on_measurement_and_non_unitary_gates() -> None:
	q0 = vqubit.LineQubit(0)
	measured = vqubit.merge_single_qubit_gates_to_phxz(vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m")))
	operations = list(measured.all_operations())
	assert len(operations) == 2
	assert isinstance(operations[0].gate, vqubit.PhasedXZGate)
	assert vqubit.is_measurement(operations[1])
	reset_circuit = vqubit.Circuit(vqubit.H(q0), vqubit.reset(q0), vqubit.X(q0))
	with_reset = vqubit.merge_single_qubit_gates_to_phxz(reset_circuit)
	assert len(list(with_reset.all_operations())) == 3
	assert isinstance(list(with_reset.all_operations())[1].gate, vqubit.ResetChannel)


def test_merge_single_qubit_gates_in_parallel_and_symbolized_variants() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit.from_moments(
		[vqubit.Moment([vqubit.H(q0), vqubit.Y(q1)]), vqubit.Moment([vqubit.X(q0), vqubit.Z(q1)])]
	)
	merged = vqubit.merge_single_qubit_gates_to_phxz(circuit)
	operations = list(merged.all_operations())
	assert len(operations) == 2
	assert all(isinstance(operation.gate, vqubit.PhasedXZGate) for operation in operations)
	assert up_to_phase(vqubit.unitary(merged), vqubit.unitary(circuit))
	assert vqubit.merge_single_qubit_gates_to_phxz_symbolized(circuit) == merged
	tagged = vqubit.symbolize_single_qubit_gates_by_indexed_tags(
		vqubit.Circuit(vqubit.H(q0), vqubit.X(q0)), unused=True
	)
	assert [_tags(operation) for operation in tagged.all_operations()] == [("sym0",), ("sym1",)]


def test_add_dynamical_decoupling_and_map_clean_are_no_ops() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	assert vqubit.add_dynamical_decoupling(circuit) == circuit
	assert vqubit.add_dynamical_decoupling(circuit, context=vqubit.TransformerContext(deep=True)) == circuit
	assert vqubit.map_clean_and_borrowable_qubits(circuit) == circuit
	assert vqubit.map_clean_and_borrowable_qubits(circuit, unused=True) == circuit


def test_eject_z_moves_rotations_through_commuting_gates() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.Z(q0), vqubit.CZ(q0, q1), vqubit.CZ(q1, q2))
	ejected = vqubit.eject_z(circuit)
	operations = list(ejected.all_operations())
	assert isinstance(operations[0].gate, vqubit.CZPowGate)
	assert operations[0].qubits == (q0, q1)
	assert isinstance(operations[1].gate, vqubit.ZPowGate)
	assert operations[1].qubits == (q0,)
	assert up_to_phase(vqubit.unitary(ejected), vqubit.unitary(circuit))
	assert vqubit.eject_phased_paulis(circuit) == ejected


def test_eject_z_leaves_non_commuting_and_leading_rz_gates_alone() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	non_commuting = vqubit.Circuit(vqubit.Z(q0), vqubit.X(q0))
	assert list(vqubit.eject_z(non_commuting).all_operations()) == list(non_commuting.all_operations())
	leading_rz = vqubit.Circuit(vqubit.rz(0.3)(q0), vqubit.CZ(q0, q1))
	assert list(vqubit.eject_z(leading_rz).all_operations()) == list(leading_rz.all_operations())
	assert up_to_phase(vqubit.unitary(vqubit.eject_z(leading_rz)), vqubit.unitary(leading_rz))


def test_defer_measurements_uses_a_fresh_ancilla() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	result = vqubit.defer_measurements(circuit)
	operations = list(result.all_operations())
	assert isinstance(operations[0].gate, vqubit.HPowGate)
	assert len(operations[1].qubits) == 2
	assert operations[1].qubits[0] == q0
	assert vqubit.is_measurement(operations[-1])
	measurement = operations[-1].gate
	assert isinstance(measurement, vqubit.MeasurementGate)
	assert str(measurement.key) == "m"
	assert set(result.all_qubits()) == {q0, vqubit.LineQubit(1)}


def test_defer_measurements_handles_multi_qubit_measurements() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.measure(q0, q1, key="both"))
	result = vqubit.defer_measurements(circuit)
	measurements = [operation for operation in result.all_operations() if operation.gate._is_measurement_()]
	# One deferred measurement, on the two ancillas and with the original key: a
	# single-qubit measurement per ancilla would repeat the key and lose the
	# classical bit ordering.
	assert len(measurements) == 1
	measurement = measurements[0].gate
	assert isinstance(measurement, vqubit.MeasurementGate)
	assert str(measurement.key) == "both"
	assert len(measurements[0].qubits) == 2
	assert len(result.all_qubits()) == 4
	no_measurement = vqubit.Circuit(vqubit.H(q0))
	assert vqubit.defer_measurements(no_measurement) == no_measurement


def test_transformer_logger_collects_and_formats_messages() -> None:
	logger = vqubit.TransformerLogger(tags_to_ignore=["a"])
	logger.log("plain")
	logger.log("value {}", 3)
	assert logger.messages() == ["plain", "value 3"]
	messages = logger.messages()
	messages.append("extra")
	assert logger.messages() == ["plain", "value 3"]
	assert logger.tags_to_ignore == ("a",)
	assert logger.context is None
	assert "2 messages" in repr(logger)


def test_transformer_context_carries_the_deep_flag_and_logger() -> None:
	logger = vqubit.TransformerLogger()
	context = vqubit.TransformerContext(deep=True, logger=logger, tags_to_ignore=["t"])
	assert context.deep is True
	assert context.logger is logger
	assert context.tags_to_ignore == ("t",)
	assert "deep=True" in repr(context)
	default = vqubit.TransformerContext()
	assert default.deep is False
	assert default.logger.context is default


def test_transformer_decorator_injects_a_context_and_forwards_kwargs() -> None:
	seen: dict[str, object] = {}

	@vqubit.transformer
	def record(circuit: vqubit.Circuit, *, context: vqubit.TransformerContext, times: int = 1, deep: bool = False):
		seen["deep"] = context.deep
		seen["is_context"] = isinstance(context, vqubit.TransformerContext)
		return vqubit.Circuit(list(circuit.all_operations()) * times)

	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	assert len(list(record(circuit, times=2).all_operations())) == 2
	assert seen == {"deep": False, "is_context": True}
	record(circuit, context=vqubit.TransformerContext(deep=True))
	assert seen["deep"] is True
	record(circuit, deep=True)
	assert seen["deep"] is True
	wrapped = getattr(record, "__wrapped__", None)
	assert wrapped is not None
	assert wrapped.__name__ == "record"


def test_create_transformer_with_kwargs_binds_and_allows_overrides() -> None:
	def repeat(circuit: vqubit.Circuit, *, context: object, n: int = 1) -> vqubit.Circuit:
		del context
		return vqubit.Circuit(list(circuit.all_operations()) * n)

	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	transformer = vqubit.create_transformer_with_kwargs(repeat, {"n": 2})
	assert len(list(transformer(circuit).all_operations())) == 2
	assert len(list(transformer(circuit, n=3).all_operations())) == 3
	assert len(list(transformer(circuit, context=vqubit.TransformerContext()).all_operations())) == 2
