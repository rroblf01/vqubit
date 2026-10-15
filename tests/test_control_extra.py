"""Tests for controlled operations, control values and Pauli interactions."""

from __future__ import annotations

import pytest

import vqubit
from tests.helpers import up_to_phase


def test_controlled_operation() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.ControlledOperation([q0], vqubit.X(q1))
	assert operation.qubits == (q0, q1)
	assert up_to_phase(vqubit.unitary(operation), vqubit.unitary(vqubit.CNOT))
	flipped = vqubit.ControlledOperation([q0], vqubit.X(q1), control_values=[False])
	assert vqubit.is_unitary(vqubit.unitary(flipped))
	assert not vqubit.allclose_up_to_global_phase(vqubit.unitary(flipped), vqubit.unitary(vqubit.CNOT))


def test_sum_of_products() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	formula = vqubit.SumOfProducts([q0, q1], [((q0, True), (q1, True)), ((q0, False), (q1, False))])
	assert (True, True) in formula
	assert (False, False) in formula
	assert (True, False) not in formula
	formula.validate((True, True))
	with pytest.raises(ValueError):
		formula.validate((True, False))


def test_product_of_sums() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	formula = vqubit.ProductOfSums([q0, q1], [((q0, True), (q1, True)), ((q0, False), (q1, True))])
	assert (True, True) in formula
	assert (False, True) in formula
	assert (True, False) not in formula


def test_abstract_control_values_is_base() -> None:
	assert issubclass(vqubit.SumOfProducts, vqubit.AbstractControlValues)
	assert issubclass(vqubit.ProductOfSums, vqubit.AbstractControlValues)


def test_pauli_interaction_gate() -> None:
	gate = vqubit.PauliInteractionGate(vqubit.Pauli.Z, False, vqubit.Pauli.X, False)
	assert gate.num_qubits == 2
	assert up_to_phase(vqubit.unitary(gate), vqubit.kron(vqubit.unitary(vqubit.Z), vqubit.unitary(vqubit.X)))
	inverted = vqubit.PauliInteractionGate(vqubit.Pauli.Z, True, vqubit.Pauli.X, False)
	assert vqubit.allclose_up_to_global_phase(vqubit.unitary(inverted), vqubit.unitary(gate))
