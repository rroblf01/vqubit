"""Tests for operations and moments."""

from __future__ import annotations

import pytest

import vqubit


def test_gate_operation_equality() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.X.on(q0) == vqubit.X.on(q0)
	assert vqubit.X.on(q0) != vqubit.Y.on(q0)
	assert vqubit.X.on(q0) != vqubit.X.on(vqubit.LineQubit(1))


def test_gate_operation_requires_qids() -> None:
	with pytest.raises(TypeError):
		vqubit.GateOperation(vqubit.X, (0,))  # ty: ignore[invalid-argument-type]


def test_moment_disjoint_qubits() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	moment = vqubit.Moment([vqubit.X(q0), vqubit.Y(q1)])
	assert moment.qubits == (q0, q1)
	assert len(moment) == 2
	assert list(moment) == [vqubit.X(q0), vqubit.Y(q1)]


def test_moment_rejects_overlap() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.Moment([vqubit.X(q0), vqubit.Y(q0)])


def test_moment_equality_ignores_order() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.Moment([vqubit.X(q0), vqubit.Y(q1)]) == vqubit.Moment([vqubit.Y(q1), vqubit.X(q0)])


def test_moment_with_operation() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	moment = vqubit.Moment([vqubit.X(q0)])
	extended = moment.with_operation(vqubit.Y(q1))
	assert extended.qubits == (q0, q1)
	assert moment.qubits == (q0,)


def test_moment_without_operations_touching() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	moment = vqubit.Moment([vqubit.X(q0), vqubit.Y(q1)])
	trimmed = moment.without_operations_touching([q0])
	assert trimmed.qubits == (q1,)


def test_empty_moment_is_falsy() -> None:
	assert not vqubit.Moment()
	assert bool(vqubit.Moment([vqubit.X(vqubit.LineQubit(0))]))
