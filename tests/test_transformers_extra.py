"""Tests for the extra circuit transformers."""

from __future__ import annotations

import vqubit
from tests.helpers import up_to_phase


def test_transformer_decorator() -> None:
	@vqubit.transformer
	def identity_transformer(circuit, *, context):
		assert isinstance(context, vqubit.TransformerContext)
		return circuit

	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	assert identity_transformer(circuit) == circuit
	assert identity_transformer(circuit, context=vqubit.TransformerContext(deep=True)) == circuit


def test_create_transformer_with_kwargs() -> None:
	def doubler(circuit, *, context, n=1):
		del context
		return vqubit.Circuit(list(circuit.all_operations()) * n)

	transformer = vqubit.create_transformer_with_kwargs(doubler, {"n": 3})
	q0 = vqubit.LineQubit(0)
	assert len(list(transformer(vqubit.Circuit(vqubit.H(q0))).all_operations())) == 3


def test_remove_and_index_tags() -> None:
	q0 = vqubit.LineQubit(0)
	tagged = vqubit.toggle_tags(vqubit.Circuit(vqubit.H(q0)), "a")
	untagged = vqubit.remove_tags(tagged, ["a"])
	assert next(iter(untagged.all_operations())) == vqubit.H(q0)
	indexed = vqubit.index_tags(vqubit.Circuit(vqubit.H(q0), vqubit.X(q0)), prefix="t")
	first = next(iter(indexed.all_operations()))
	assert isinstance(first, vqubit.TaggedOperation)


def test_is_negligible_turn() -> None:
	assert vqubit.is_negligible_turn(0.0)
	assert vqubit.is_negligible_turn(1.0)
	assert not vqubit.is_negligible_turn(0.5)


def test_unroll_circuit_op() -> None:
	q0 = vqubit.LineQubit(0)
	inner = vqubit.CircuitOperation(vqubit.Circuit(vqubit.H(q0)))
	unrolled = vqubit.unroll_circuit_op(vqubit.Circuit(inner))
	assert len(list(unrolled.all_operations())) == 1
	repeated = vqubit.CircuitOperation(vqubit.Circuit(vqubit.H(q0)), repetitions=3)
	assert len(list(vqubit.unroll_circuit_op(vqubit.Circuit(repeated)).all_operations())) == 3


def test_map_operations_and_unroll() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.CircuitOperation(vqubit.Circuit(vqubit.H(q0))))
	result = vqubit.map_operations_and_unroll(circuit, lambda op: op)
	assert len(list(result.all_operations())) == 1


def test_drop_terminal_measurements() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="m"))
	result = vqubit.drop_terminal_measurements(circuit)
	assert len(list(result.all_operations())) == 2


def test_drop_diagonal_before_measurement() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.Z(q0), vqubit.measure(q0, key="m"))
	assert len(list(vqubit.drop_diagonal_before_measurement(circuit).all_operations())) == 1


def test_merge_operations() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.X(q0))

	def merge_func(op1, op2):
		matrix = vqubit.dot(vqubit.unitary(op2.gate), vqubit.unitary(op1.gate))
		return vqubit.MatrixGate(matrix).on(op1.qubits[0])

	merged = vqubit.merge_operations(circuit, merge_func)
	assert len(list(merged.all_operations())) == 1
	assert up_to_phase(vqubit.unitary(merged), vqubit.unitary(circuit))


def test_merge_single_qubit_gates() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.T(q0))
	phxz = vqubit.merge_single_qubit_gates_to_phxz(circuit)
	assert len(list(phxz.all_operations())) == 1
	assert up_to_phase(vqubit.unitary(phxz), vqubit.unitary(circuit))
	phased_xz = vqubit.merge_single_qubit_gates_to_phased_x_and_z(circuit)
	assert up_to_phase(vqubit.unitary(phased_xz), vqubit.unitary(circuit))


def test_eject_z() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.Z(q0), vqubit.CZ(q0, q1))
	result = vqubit.eject_z(circuit)
	assert up_to_phase(vqubit.unitary(result), vqubit.unitary(circuit))


def test_defer_measurements() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	result = vqubit.defer_measurements(circuit)
	operations = list(result.all_operations())
	assert operations[-1].gate._is_measurement_()


def test_no_op_transformers() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	assert vqubit.map_clean_and_borrowable_qubits(circuit) == circuit
	assert vqubit.add_dynamical_decoupling(circuit) == circuit
