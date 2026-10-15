"""Tests for the sweep algebra and run_batch."""

from __future__ import annotations

import math

import pytest

import vqubit


def test_single_sweep() -> None:
	sweep = vqubit.SingleSweep("t", [0.0, 1.0])
	assert [resolver["t"] for resolver in sweep.param_resolvers] == [0.0, 1.0]
	assert len(sweep) == 2


def test_product() -> None:
	sweep = vqubit.Product(vqubit.SingleSweep("a", [0, 1]), vqubit.SingleSweep("b", [2, 3]))
	assert len(sweep) == 4
	pairs = [(resolver["a"], resolver["b"]) for resolver in sweep.param_resolvers]
	assert pairs == [(0, 2), (0, 3), (1, 2), (1, 3)]


def test_zip() -> None:
	sweep = vqubit.Zip(vqubit.SingleSweep("a", [0, 1]), vqubit.SingleSweep("b", [2, 3]))
	assert [(resolver["a"], resolver["b"]) for resolver in sweep.param_resolvers] == [(0, 2), (1, 3)]


def test_zip_requires_equal_lengths() -> None:
	with pytest.raises(ValueError):
		_ = vqubit.Zip(vqubit.SingleSweep("a", [0, 1]), vqubit.SingleSweep("b", [2])).param_resolvers


def test_concat() -> None:
	sweep = vqubit.Concat(vqubit.SingleSweep("a", [0, 1]), vqubit.SingleSweep("a", [2]))
	assert [resolver["a"] for resolver in sweep.param_resolvers] == [0, 1, 2]


def test_to_sweep() -> None:
	assert len(vqubit.to_sweep({"a": [0, 1, 2]})) == 3
	assert len(vqubit.to_sweep(vqubit.ParamResolver({"a": 1}))) == 1
	assert len(vqubit.to_sweep(vqubit.Points([{"a": 0}, {"a": 1}]))) == 2
	with pytest.raises(TypeError):
		vqubit.to_sweep(5)


def test_run_batch_keeps_parameters() -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.Rx(t)(q0), vqubit.measure(q0, key="m"))
	results = vqubit.run_batch(circuit, vqubit.Points([{"t": 0.0}, {"t": math.pi}]), repetitions=8)
	assert len(results) == 2
	assert results[0].params["t"] == 0.0
	assert results[0].histogram(key="m") == {0: 8}
	assert results[1].histogram(key="m") == {1: 8}


def test_run_batch_over_product() -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.Rx(t)(q0), vqubit.measure(q0, key="m"))
	results = vqubit.run_batch(circuit, vqubit.Product(vqubit.SingleSweep("t", [0.0, math.pi])), repetitions=2)
	assert len(results) == 2
