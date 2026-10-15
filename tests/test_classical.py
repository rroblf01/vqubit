"""Tests for classically controlled operations and classical values."""

from __future__ import annotations

import vqubit


def test_key_condition() -> None:
	condition = vqubit.KeyCondition("m", index=1)
	assert condition.key == "m"
	assert condition.index == 1
	negated = ~condition
	assert negated == ~condition
	assert ~negated == condition


def test_classically_controlled_operation() -> None:
	q0 = vqubit.LineQubit(0)
	operation = vqubit.ClassicallyControlledOperation([vqubit.KeyCondition("m")], vqubit.X(q0))
	assert operation.qubits == (q0,)
	assert operation.sub_operation == vqubit.X(q0)
	assert len(operation.conditions) == 1
	assert not vqubit.has_unitary(operation)


def test_if_operation() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.If(vqubit.KeyCondition("m"), vqubit.X(q0))
	assert isinstance(result, vqubit.ClassicallyControlledOperation)
	both = vqubit.If(vqubit.KeyCondition("m"), vqubit.X(q0), else_ops=[vqubit.Z(q0)])
	assert isinstance(both, tuple)
	assert len(both) == 2


def test_set_variable() -> None:
	operation = vqubit.SetVariable("v", 5)
	assert operation.key == "v"
	assert operation.value == 5
	assert operation.qubits == ()
