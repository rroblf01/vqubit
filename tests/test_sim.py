"""Tests for the state-vector simulator."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_vector_close

SQRT_HALF = math.sqrt(0.5)


def test_simulate_x_on_zero() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.X(q0)))
	assert_vector_close(result.state_vector(), [0j, 1 + 0j])


def test_simulate_bell_state() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	result = vqubit.Simulator().simulate(circuit)
	assert_vector_close(result.state_vector(), [SQRT_HALF + 0j, 0j, 0j, SQRT_HALF + 0j])
	assert result.dirac_notation() == "0.707|00⟩ + 0.707|11⟩"


def test_simulate_with_initial_basis_state() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	result = vqubit.Simulator().simulate(vqubit.Circuit(), qubit_order=[q0, q1], initial_state=3)
	assert_vector_close(result.state_vector(), [0j, 0j, 0j, 1 + 0j])


def test_simulate_rejects_measurements() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0))
	with pytest.raises(ValueError):
		vqubit.Simulator().simulate(circuit)


def test_run_counts_have_the_right_total() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	result = vqubit.Simulator(seed=2).run(circuit, repetitions=500)
	histogram = result.histogram(key="m")
	assert sum(histogram.values()) == 500
	assert set(histogram).issubset({0, 1})


def test_run_is_deterministic_for_an_eigenstate() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	result = vqubit.Simulator().run(circuit, repetitions=10)
	assert result.histogram(key="m") == {1: 10}


def test_run_bell_state_only_produces_00_and_11() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="result"))
	result = vqubit.Simulator(seed=3).run(circuit, repetitions=200)
	assert set(result.histogram(key="result")).issubset({0, 3})


def test_run_respects_invert_mask() -> None:
	q0 = vqubit.LineQubit(0)
	gate = vqubit.MeasurementGate("m", num_qubits=1, invert_mask=(True,))
	circuit = vqubit.Circuit(gate.on(q0))
	result = vqubit.Simulator().run(circuit, repetitions=5)
	assert result.histogram(key="m") == {1: 5}


def test_qubit_order_is_respected() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	# X on q1 leaves |01> in the (q0, q1) order, which is basis index 1.
	result = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.X(q1)), qubit_order=[q0, q1])
	assert_vector_close(result.state_vector(), [0j, 1 + 0j, 0j, 0j])


def test_qubit_order_must_cover_the_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	with pytest.raises(ValueError):
		vqubit.Simulator().simulate(vqubit.Circuit(vqubit.CNOT(q0, q1)), qubit_order=[q0])


def test_repetitions_must_be_positive() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.Simulator().run(vqubit.Circuit(vqubit.measure(q0)), repetitions=0)


def test_state_vector_validates_length() -> None:
	with pytest.raises(ValueError):
		vqubit.StateVector([1 + 0j], [vqubit.LineQubit(0), vqubit.LineQubit(1)])


def test_result_default_histogram_with_single_key() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.Simulator().run(vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0)))
	assert result.histogram() == {1: 1}
	assert len(result) == 1
