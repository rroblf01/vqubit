"""Tests for the density-matrix simulator."""

from __future__ import annotations

import math

import pytest

import vqubit

SQRT_HALF = math.sqrt(0.5)


def test_pure_state_matches_state_vector() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	vector = vqubit.Simulator().simulate(circuit).state_vector()
	density = vqubit.DensityMatrixSimulator().simulate(circuit).density_matrix()
	for i in range(4):
		for j in range(4):
			assert abs(density[i][j] - vector[i] * vector[j].conjugate()) < 1e-9


def test_bell_state_density_matrix() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	density = vqubit.DensityMatrixSimulator().simulate(circuit).density_matrix()
	assert density[0][0] == pytest.approx(0.5)
	assert density[3][3] == pytest.approx(0.5)
	assert density[0][3] == pytest.approx(0.5)
	assert density[1][1] == pytest.approx(0)


def test_depolarizing_lowers_coherence() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	noisy = vqubit.Circuit(vqubit.H(q0), vqubit.depolarize(0.1)(q0), vqubit.CNOT(q0, q1))
	density = vqubit.DensityMatrixSimulator().simulate(noisy).density_matrix()
	assert density[0][3].real < 0.5
	assert sum(density[i][i].real for i in range(4)) == pytest.approx(1.0)


def test_amplitude_damping_relaxes_to_ground() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.amplitude_damp(1.0)(q0))
	density = vqubit.DensityMatrixSimulator().simulate(circuit).density_matrix()
	assert density[0][0] == pytest.approx(1.0)
	assert density[1][1] == pytest.approx(0.0)


def test_channel_on_one_of_two_qubits() -> None:
	q1 = vqubit.LineQubit(1)
	# Flip q1 fully; the density matrix should be |01> in (q0, q1) order.
	circuit = vqubit.Circuit(vqubit.bit_flip(1.0)(q1))
	density = vqubit.DensityMatrixSimulator().simulate(circuit).density_matrix()
	assert density[1][1] == pytest.approx(1.0)


def test_run_samples_measurements() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	result = vqubit.DensityMatrixSimulator(seed=1).run(circuit, repetitions=400)
	histogram = result.histogram(key="m")
	assert sum(histogram.values()) == 400
	assert set(histogram).issubset({0, 1})


def test_run_is_deterministic_for_a_basis_state() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	result = vqubit.DensityMatrixSimulator().run(circuit, repetitions=5)
	assert result.histogram(key="m") == {1: 5}


def test_simulate_rejects_measurements() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.DensityMatrixSimulator().simulate(vqubit.Circuit(vqubit.measure(q0)))


def test_trial_result_exposes_the_density_matrix() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.DensityMatrixSimulator().simulate(vqubit.Circuit(vqubit.H(q0)))
	assert result.final_density_matrix.qubits == (q0,)
	assert result.density_matrix()[0][0] == pytest.approx(0.5)


def test_repetitions_must_be_positive() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.DensityMatrixSimulator().run(vqubit.Circuit(vqubit.measure(q0)), repetitions=0)
