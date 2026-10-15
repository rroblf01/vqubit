"""Tests for the extra gate families, parallel gates and tags."""

from __future__ import annotations

import vqubit
from tests.helpers import assert_matrix_close


def test_parallel_gate() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gate = vqubit.ParallelGate(vqubit.X)
	assert gate.num_qubits == 2
	assert_matrix_close(vqubit.unitary(gate), vqubit.kron(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.X)))
	operation = vqubit.parallel_gate_op(vqubit.H, q0, q1)
	assert operation.gate.num_qubits == 2


def test_any_unitary_gate_family() -> None:
	q0 = vqubit.LineQubit(0)
	family = vqubit.AnyUnitaryGateFamily()
	assert vqubit.X(q0) in family
	assert vqubit.bit_flip(0.1)(q0) not in family


def test_any_integer_power_gate_family() -> None:
	q0 = vqubit.LineQubit(0)
	family = vqubit.AnyIntegerPowerGateFamily(vqubit.X)
	assert (vqubit.X**2)(q0) in family
	assert (vqubit.X**-1)(q0) in family
	assert (vqubit.X**0.5)(q0) not in family
	assert vqubit.H(q0) not in family


def test_parallel_gate_family() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	family = vqubit.ParallelGateFamily(vqubit.X)
	assert vqubit.ParallelGate(vqubit.X)(q0, q1) in family
	assert vqubit.X(q0) in family
	assert vqubit.H(q0) not in family


def test_tags() -> None:
	assert vqubit.VirtualTag() == vqubit.VirtualTag()
	assert vqubit.VirtualTag() != vqubit.RoutingSwapTag()
	assert repr(vqubit.VirtualTag()) == "vqubit.VirtualTag()"
	assert isinstance(vqubit.InterchangeableQubitsGate, type)
