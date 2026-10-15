"""Tests for the approximate-equality and measurement-key protocols."""

from __future__ import annotations

import vqubit


def test_approx_eq() -> None:
	assert vqubit.approx_eq(vqubit.X, vqubit.X)
	assert not vqubit.approx_eq(vqubit.X, vqubit.Y)
	assert vqubit.approx_eq(vqubit.CNOT, vqubit.CX)


def test_approx_eq_falls_back_to_equality() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.approx_eq(q0, vqubit.LineQubit(0))
	assert not vqubit.approx_eq(q0, vqubit.LineQubit(1))


def test_equal_up_to_global_phase() -> None:
	minus_identity = vqubit.MatrixGate([[-1 + 0j, 0j], [0j, -1 + 0j]])
	assert vqubit.equal_up_to_global_phase(vqubit.I, minus_identity)
	assert not vqubit.equal_up_to_global_phase(vqubit.I, vqubit.X)


def test_measurement_key_obj() -> None:
	q0 = vqubit.LineQubit(0)
	key = vqubit.measurement_key_obj(vqubit.measure(q0, key="alpha"))
	assert key is not None
	assert key.name == "alpha"
	assert str(key) == "alpha"
	assert vqubit.measurement_key_obj(vqubit.H) is None


def test_measurement_key_objs() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="a"), vqubit.X(q1), vqubit.measure(q1, key="b"))
	keys = vqubit.measurement_key_objs(circuit)
	assert {key.name for key in keys} == {"a", "b"}


def test_with_measurement_key_mapping() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="old"))
	renamed = vqubit.with_measurement_key_mapping(circuit, {"old": "new"})
	gate = next(iter(renamed.all_operations())).gate
	assert isinstance(gate, vqubit.MeasurementGate)
	assert gate.key == "new"


def test_control_keys_is_empty() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.control_keys(vqubit.measure(q0, key="m")) == frozenset()


def test_trace_distance_bound() -> None:
	assert vqubit.trace_distance_bound(vqubit.I) == 0.0
	assert vqubit.trace_distance_bound(vqubit.X) == 1.0
