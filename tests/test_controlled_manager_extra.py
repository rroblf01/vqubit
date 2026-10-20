"""Tests for controlled operations and the qubit managers."""

from __future__ import annotations

import pytest

import vqubit
from tests.helpers import assert_matrix_close


def test_controlled_operation_reports_its_controls_target_and_values() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1))

	assert operation.controls == (q0,)
	assert operation.sub_operation == vqubit.X(q1)
	assert operation.qubits == (q0, q1)
	assert operation.control_values == (True,)
	assert operation.gate.num_qubits == 2


def test_controlled_operation_unitary_is_a_cnot() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1))

	assert_matrix_close(vqubit.unitary(operation), vqubit.unitary(vqubit.CNOT))


def test_controlled_operation_exposes_the_unitary_protocol() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1))

	assert_matrix_close(operation._unitary_(), vqubit.unitary(vqubit.CNOT))


def test_controlled_operation_fires_on_a_zero_control() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1), [0])
	expected = vqubit.X.controlled_by(q0, control_values=[0]).on(q0, q1)

	assert operation.control_values == (False,)
	assert_matrix_close(vqubit.unitary(operation), vqubit.unitary(expected))


def test_control_values_are_coerced_to_bools() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	operation = vqubit.ControlledOperation([q0, q1], vqubit.X(q2), [0, 1])

	assert operation.control_values == (False, True)
	assert operation.gate == vqubit.ControlledGate(vqubit.X, num_controls=2, control_values=(False, True))


def test_control_formula_with_one_satisfying_assignment() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	# ``q0 & !q1``: exactly one assignment, so it maps onto a single control value set.
	values = vqubit.SumOfProducts([q0, q1], [[(q0, True), (q1, False)]])
	operation = vqubit.ControlledOperation([q0, q1], vqubit.X(q2), values)
	expected = vqubit.X.controlled_by(q0, q1, control_values=[1, 0]).on(q0, q1, q2)

	assert operation.control_values is values
	assert_matrix_close(vqubit.unitary(operation), vqubit.unitary(expected))


def test_control_formula_without_a_satisfying_assignment_raises() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1), vqubit.SumOfProducts([q0], []))

	with pytest.raises(ValueError, match="no satisfying control assignment"):
		_ = operation.gate


def test_controlled_operation_gate_property() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1))

	assert operation.gate == vqubit.ControlledGate(vqubit.X, num_controls=1, control_values=(True,))
	assert_matrix_close(vqubit.unitary(operation.gate), vqubit.unitary(vqubit.CNOT))


def test_with_qubits_remaps_controls_first_then_target() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1))

	moved = operation.with_qubits(q2, q0)

	assert moved.controls == (q2,)
	assert moved.sub_operation == vqubit.X(q0)
	assert moved.qubits == (q2, q0)
	assert_matrix_close(vqubit.unitary(moved), vqubit.unitary(vqubit.CNOT))


def test_with_qubits_rejects_the_wrong_number_of_qubits() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1))

	with pytest.raises(ValueError):
		operation.with_qubits(q2)
	with pytest.raises(ValueError):
		operation.with_qubits(q2, q0, q1)


def test_controlled_operation_equality_and_hashing() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1))

	assert operation == vqubit.ControlledOperation([q0], vqubit.X(q1))
	assert operation != vqubit.ControlledOperation([q0], vqubit.X(q1), [0])
	assert operation != vqubit.ControlledOperation([q1], vqubit.X(q1))
	assert len({operation, vqubit.ControlledOperation([q0], vqubit.X(q1))}) == 1


def test_controlled_operation_str_and_repr() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1))

	assert str(operation) == "ControlledOperation((vqubit.LineQubit(0),), XPowGate(exponent=1))"
	assert repr(operation) == (
		"vqubit.ControlledOperation((vqubit.LineQubit(0),), "
		"vqubit.GateOperation(gate=vqubit.XPowGate(exponent=1.0, global_shift=0.0), "
		"qubits=(vqubit.LineQubit(1),)), control_values=(True,))"
	)


def test_nested_controlled_operations_give_a_toffoli() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	nested = vqubit.ControlledOperation([q0], vqubit.ControlledOperation([q1], vqubit.X(q2)))

	assert nested.qubits == (q0, q1, q2)
	assert_matrix_close(vqubit.unitary(nested), vqubit.unitary(vqubit.CCX))


def test_base_qubit_manager_operations_are_not_implemented() -> None:
	q0 = vqubit.LineQubit(0)
	manager = vqubit.QubitManager()

	with pytest.raises(NotImplementedError):
		manager.qalloc()
	with pytest.raises(NotImplementedError):
		manager.qfree(q0)
	with pytest.raises(NotImplementedError):
		manager.qborrow()
	with pytest.raises(NotImplementedError):
		manager.qrelease(q0)


def test_base_manager_allocate_and_free_delegate_to_the_hooks() -> None:
	class RecordingManager(vqubit.QubitManager):
		def __init__(self) -> None:
			self.allocated: list[vqubit.Qid] = []
			self.freed: list[vqubit.Qid] = []

		def qalloc(self) -> vqubit.Qid:
			qubit = vqubit.LineQubit(len(self.allocated))
			self.allocated.append(qubit)
			return qubit

		def qfree(self, qubit: vqubit.Qid) -> None:
			self.freed.append(qubit)

	manager = RecordingManager()

	qubits = manager.allocate(2)
	manager.free(qubits)

	assert qubits == [vqubit.LineQubit(0), vqubit.LineQubit(1)]
	assert manager.freed == qubits


def test_borrowable_qubit_strings_repr_order_and_equality() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	borrowed = vqubit.BorrowableQubit(q0)

	assert borrowed.qubit == q0
	assert str(borrowed) == "borrow(0)"
	assert repr(borrowed) == "vqubit.BorrowableQubit(vqubit.LineQubit(0))"
	assert borrowed == vqubit.BorrowableQubit(q0)
	assert borrowed != vqubit.BorrowableQubit(q1)
	assert borrowed < vqubit.BorrowableQubit(q1)
	assert hash(borrowed) == hash(vqubit.BorrowableQubit(q0))


def test_clean_qubit_strings_repr_order_and_equality() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	clean = vqubit.CleanQubit(q0)

	assert clean.qubit == q0
	assert isinstance(clean, vqubit.Qid)
	assert str(clean) == "clean(0)"
	assert repr(clean) == "vqubit.CleanQubit(vqubit.LineQubit(0))"
	assert clean == vqubit.CleanQubit(q0)
	assert clean != vqubit.CleanQubit(q1)
	assert clean < vqubit.CleanQubit(q1)


def test_a_borrowed_qubit_is_not_a_clean_qubit() -> None:
	q0 = vqubit.LineQubit(0)

	assert vqubit.BorrowableQubit(q0) != vqubit.CleanQubit(q0)


def test_simple_manager_allocates_after_the_provided_qubits() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	manager = vqubit.SimpleQubitManager(qubits=[q0, q1])

	assert manager.qubits == [q0, q1]
	assert manager.qalloc() == vqubit.LineQubit(2)
	assert manager.qubits == [q0, q1]


def test_simple_manager_qubits_property_is_a_copy() -> None:
	q0 = vqubit.LineQubit(0)
	manager = vqubit.SimpleQubitManager(qubits=[q0])

	manager.qubits.append(vqubit.LineQubit(9))

	assert manager.qubits == [q0]


def test_simple_manager_free_and_release_are_noops() -> None:
	manager = vqubit.SimpleQubitManager()

	first = manager.qalloc()
	manager.qfree(first)
	manager.qrelease(first)

	assert manager.qalloc() == vqubit.LineQubit(1)
	borrowed = manager.qborrow()
	assert isinstance(borrowed, vqubit.BorrowableQubit)
	assert borrowed.qubit == vqubit.LineQubit(2)


def test_simple_manager_counts_existing_qubits_against_max_size() -> None:
	q0 = vqubit.LineQubit(0)
	manager = vqubit.SimpleQubitManager(qubits=[q0], max_size=2)

	assert manager.qalloc() == vqubit.LineQubit(1)
	with pytest.raises(ValueError, match="max_size"):
		manager.qalloc()


def test_greedy_manager_recycles_provided_qubits_first_in_first_out() -> None:
	# Cirq's default (maximize_reuse=False) is FIFO: the qubit freed earliest
	# comes back first. ``maximize_reuse=True`` switches to LIFO.
	q0, q1 = vqubit.LineQubit.range(2)
	manager = vqubit.GreedyQubitManager(qubits=[q0, q1])

	assert manager.qalloc() == q0
	assert manager.qalloc() == q1
	assert manager.qalloc() == vqubit.LineQubit(2)

	lifo = vqubit.GreedyQubitManager(qubits=[q0, q1], maximize_reuse=True)
	assert lifo.qalloc() == q1
	assert lifo.qalloc() == q0


def test_greedy_manager_reuses_a_freed_qubit_only_once() -> None:
	manager = vqubit.GreedyQubitManager()
	manager.qfree(vqubit.LineQubit(5))  # Never allocated: ignored.

	qubit = manager.qalloc()
	manager.qfree(qubit)
	manager.qfree(qubit)  # A second free must not duplicate it in the pool.

	assert manager.qalloc() == qubit
	assert manager.qalloc() == vqubit.LineQubit(1)


def test_greedy_manager_borrow_and_release() -> None:
	manager = vqubit.GreedyQubitManager(max_size=1)

	borrowed = manager.qborrow()
	assert isinstance(borrowed, vqubit.BorrowableQubit)
	assert borrowed.qubit == vqubit.LineQubit(0)

	manager.qrelease(borrowed.qubit)
	assert manager.qalloc() == vqubit.LineQubit(0)


def test_greedy_manager_enforces_max_size() -> None:
	manager = vqubit.GreedyQubitManager(max_size=1)

	qubit = manager.qalloc()
	with pytest.raises(ValueError, match="max_size"):
		manager.qalloc()

	manager.qfree(qubit)
	assert manager.qalloc() == qubit
