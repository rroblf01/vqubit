"""Tests for devices and the local engine."""

from __future__ import annotations

import cmath
import math

import pytest

import vqubit


def test_line_device_connectivity() -> None:
	device = vqubit.LineDevice(3)
	assert device.qubits == frozenset(vqubit.LineQubit.range(3))
	assert set(device.neighbors(vqubit.LineQubit(1))) == {vqubit.LineQubit(0), vqubit.LineQubit(2)}
	assert device.neighbors(vqubit.LineQubit(0)) == [vqubit.LineQubit(1)]


def test_line_device_validates_operations() -> None:
	device = vqubit.LineDevice(2)
	device.validate_operation(vqubit.CNOT(vqubit.LineQubit(0), vqubit.LineQubit(1)))
	with pytest.raises(ValueError):
		device.validate_operation(vqubit.CNOT(vqubit.LineQubit(0), vqubit.LineQubit(5)))
	with pytest.raises(ValueError):
		device.validate_operation(vqubit.CZ(vqubit.LineQubit(0), vqubit.LineQubit(2)))


def test_grid_device() -> None:
	device = vqubit.GridDevice.rect(2, 3)
	assert len(device.qubits) == 6
	assert device.is_adjacent(vqubit.GridQubit(0, 0), vqubit.GridQubit(0, 1))
	assert not device.is_adjacent(vqubit.GridQubit(0, 0), vqubit.GridQubit(1, 1))


def test_sycamore_device_and_gate() -> None:
	device = vqubit.Sycamore()
	assert len(device.qubits) == 54
	assert device.metadata.name == "Sycamore"
	gate = vqubit.SycamoreGate()
	assert gate.num_qubits == 2
	matrix = vqubit.unitary(gate)
	assert matrix[1][2] == pytest.approx(-1j)
	assert matrix[2][1] == pytest.approx(-1j)
	assert matrix[3][3] == pytest.approx(cmath.exp(-1j * math.pi / 6))


def test_ion_device_is_all_to_all() -> None:
	device = vqubit.IonDevice(4)
	assert len(device.neighbors(vqubit.LineQubit(0))) == 3
	device.validate_operation(vqubit.CNOT(vqubit.LineQubit(0), vqubit.LineQubit(3)))


def test_engine_runs_locally() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	engine = vqubit.Engine(device=vqubit.LineDevice(2), seed=0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="m"))
	result = engine.run(circuit, repetitions=200)
	assert set(result.histogram(key="m")).issubset({0, 3})


def test_engine_validates_against_its_device() -> None:
	engine = vqubit.Engine(device=vqubit.LineDevice(2))
	with pytest.raises(ValueError):
		engine.run(vqubit.Circuit(vqubit.CNOT(vqubit.LineQubit(0), vqubit.LineQubit(4))))


def test_engine_sweep_and_sampler() -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.Rx(t)(q0), vqubit.measure(q0, key="m"))
	engine = vqubit.Engine(seed=0)
	results = engine.run_sweep(circuit, vqubit.Points([{"t": 0.0}, {"t": math.pi}]), repetitions=5)
	assert results[0].histogram(key="m") == {0: 5}
	assert results[1].histogram(key="m") == {1: 5}
	sampler = vqubit.get_engine(seed=0).get_sampler()
	assert sampler.run(vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))).histogram(key="m") == {1: 1}


def test_engine_qasm() -> None:
	q0 = vqubit.LineQubit(0)
	assert "x q[0];" in vqubit.Engine().qasm(vqubit.Circuit(vqubit.X(q0)))
