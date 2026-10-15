"""Tests for op-tree helpers and gatesets."""

from __future__ import annotations

import pytest

import vqubit


def test_flatten_to_ops() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	tree = [circuit, vqubit.X(q1), [vqubit.Y(q0), None]]
	operations = vqubit.flatten_to_ops(tree)
	assert len(operations) == 4
	assert all(isinstance(op, vqubit.Operation) for op in operations)


def test_flatten_to_ops_or_moments() -> None:
	q0 = vqubit.LineQubit(0)
	items = vqubit.flatten_to_ops_or_moments([vqubit.Moment([vqubit.H(q0)]), vqubit.X(q0)])
	assert len(items) == 2
	assert isinstance(items[0], vqubit.Moment)
	assert isinstance(items[1], vqubit.Operation)


def test_freeze_op_tree() -> None:
	q0 = vqubit.LineQubit(0)
	frozen = vqubit.freeze_op_tree([vqubit.H(q0), [vqubit.X(q0)]])
	assert isinstance(frozen, tuple)
	assert len(frozen) == 2


def test_transform_op_tree() -> None:
	q0 = vqubit.LineQubit(0)
	transformed = vqubit.transform_op_tree(vqubit.H(q0), lambda op: (op.gate**-1).on(*op.qubits))
	assert transformed.gate == vqubit.H**-1


def test_flatten_rejects_non_op_tree() -> None:
	with pytest.raises(TypeError):
		vqubit.flatten_to_ops(123)


def test_gate_family() -> None:
	q0 = vqubit.LineQubit(0)
	family = vqubit.GateFamily(vqubit.X)
	assert vqubit.X(q0) in family
	assert vqubit.H(q0) not in family
	family.validate(vqubit.X(q0))
	with pytest.raises(ValueError):
		family.validate(vqubit.H(q0))
	predicate_family = vqubit.GateFamily(predicate=lambda gate: gate.num_qubits == 1)
	assert vqubit.H(q0) in predicate_family


def test_gateset() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gateset = vqubit.Gateset([vqubit.X, vqubit.CNOT])
	assert vqubit.X(q0) in gateset
	assert vqubit.CNOT(q0, q1) in gateset
	assert vqubit.H(q0) not in gateset
	bigger = gateset.with_added_gates(vqubit.H)
	assert vqubit.H(q0) in bigger
