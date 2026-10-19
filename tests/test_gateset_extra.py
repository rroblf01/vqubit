"""Tests for the extra gate families, parallel gates and tags."""

from __future__ import annotations

from typing import Any

import pytest

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


def test_a_family_takes_a_gate_type_too() -> None:
	# ``ParallelGateFamily(ZPowGate)`` means *any* ``ZPowGate``, as in Cirq, while
	# an instance keeps meaning that one gate.
	q0, q1 = vqubit.LineQubit.range(2)
	family = vqubit.ParallelGateFamily(vqubit.ZPowGate)
	assert vqubit.Z in family
	assert vqubit.S in family
	assert vqubit.T in family
	assert vqubit.X not in family
	assert vqubit.ParallelGate(vqubit.Z)(q0, q1) in family
	assert vqubit.GateFamily(vqubit.ZPowGate).gate is vqubit.ZPowGate


def test_any_integer_power_family_takes_a_class_and_skips_symbols() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	family = vqubit.AnyIntegerPowerGateFamily(vqubit.CZPowGate)
	assert vqubit.CZ(q0, q1) in family
	assert (vqubit.CZ**2)(q0, q1) in family
	assert (vqubit.CZ**-1)(q0, q1) in family
	assert (vqubit.CZ**0.5)(q0, q1) not in family
	assert vqubit.CNOT(q0, q1) not in family
	t: Any = vqubit.symbol("t")
	assert (vqubit.CZ**t)(q0, q1) not in family  # a symbolic exponent never matches


def test_any_unitary_gate_family_can_restrict_the_width() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.X(q0) in vqubit.AnyUnitaryGateFamily(1)
	assert vqubit.X(q0) not in vqubit.AnyUnitaryGateFamily(2)
	assert vqubit.CZ(q0, q1) in vqubit.AnyUnitaryGateFamily(2)
	with pytest.raises(ValueError):
		vqubit.AnyUnitaryGateFamily(0)


def test_parallel_gate_family_caps_the_width() -> None:
	family = vqubit.ParallelGateFamily(vqubit.ZPowGate, max_parallel_allowed=1)
	assert vqubit.Z in family
	assert vqubit.ParallelGate(vqubit.Z) not in family


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
