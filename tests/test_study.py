"""Tests for the study helpers (sampling, sweeps)."""

from __future__ import annotations

import math

import pytest

import vqubit


def test_points_from_dict() -> None:
	points = vqubit.Points({"t": 0.5})
	assert len(points) == 1
	assert points.param_resolvers[0]["t"] == 0.5


def test_points_from_list() -> None:
	points = vqubit.Points([{"t": 0.0}, {"t": 1.0}])
	assert [resolver["t"] for resolver in points.param_resolvers] == [0.0, 1.0]


def test_points_rejects_empty() -> None:
	with pytest.raises(ValueError):
		vqubit.Points([])


def test_linspace_values_stop_before_the_end() -> None:
	sweep = vqubit.Linspace({"t": 0.0}, {"t": 1.0}, 4)
	assert [resolver["t"] for resolver in sweep.param_resolvers] == [0.0, 0.25, 0.5, 0.75]
	assert len(sweep) == 4


def test_linspace_validates() -> None:
	with pytest.raises(ValueError):
		vqubit.Linspace({"t": 0.0}, {"t": 1.0}, 0)
	with pytest.raises(ValueError):
		vqubit.Linspace({"t": 0.0}, {"s": 1.0}, 3)


def test_sample_and_run() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	assert vqubit.sample(circuit, repetitions=5).histogram(key="m") == {1: 5}
	assert vqubit.run(circuit, repetitions=5).histogram(key="m") == {1: 5}


def test_run_sweep_applies_each_point() -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.Rx(t)(q0), vqubit.measure(q0, key="m"))
	results = vqubit.run_sweep(circuit, vqubit.Points([{"t": 0.0}, {"t": math.pi}]), repetitions=8)
	assert len(results) == 2
	assert results[0].histogram(key="m") == {0: 8}
	assert results[1].histogram(key="m") == {1: 8}


def test_run_sweep_over_linspace() -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.Rx(t)(q0), vqubit.measure(q0, key="m"))
	results = vqubit.run_sweep(circuit, vqubit.Linspace({"t": 0.0}, {"t": math.pi}, 3), repetitions=4)
	assert len(results) == 3


def test_sampler() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	sampler = vqubit.Sampler(seed=0)
	assert sampler.run(circuit, repetitions=3).histogram(key="m") == {1: 3}
	assert sampler.sample(circuit, repetitions=3).histogram(key="m") == {1: 3}


def test_run_sweep_rejects_unknown_sweep() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError):
		vqubit.run_sweep(vqubit.Circuit(vqubit.measure(q0)), "nope")  # ty: ignore[invalid-argument-type]
