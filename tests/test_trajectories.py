"""Channels on the state-vector simulator, by trajectory sampling.

``cirq.Simulator`` runs a noisy circuit by picking one Kraus branch per
repetition, so the same circuit with the same seed gives the same *distribution*
as the density-matrix simulator. These tests check that distribution and the
per-trajectory invariants (the state is pure and normalised), including a channel
whose Kraus effects are not diagonal, which needs the coherences of the target
subsystem rather than its marginals.
"""

from __future__ import annotations

import math
from collections import Counter

import pytest

import vqubit


def _histogram(rows) -> dict[int, float]:
	counts: Counter[int] = Counter()
	for row in rows:
		value = 0
		for bit in row:
			value = (value << 1) | int(bit)
		counts[value] += 1
	return {key: count / len(rows) for key, count in counts.items()}


def _gap(left: dict[int, float], right: dict[int, float]) -> float:
	return max(abs(left.get(key, 0.0) - right.get(key, 0.0)) for key in set(left) | set(right))


def test_a_noisy_circuit_samples_like_the_density_matrix() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(
		vqubit.H(q0),
		vqubit.amplitude_damp(0.4)(q0),
		vqubit.CNOT(q0, q1),
		vqubit.depolarize(0.15)(q1),
		vqubit.measure(q0, q1, key="m"),
	)
	repetitions = 4000
	trajectories = _histogram(vqubit.Simulator(seed=3).run(circuit, repetitions=repetitions).measurements["m"])
	exact = _histogram(vqubit.DensityMatrixSimulator(seed=3).run(circuit, repetitions=repetitions).measurements["m"])
	assert _gap(trajectories, exact) < 0.05


def test_every_built_in_channel_runs_on_the_state_vector_simulator() -> None:
	q0 = vqubit.LineQubit(0)
	channels = [
		vqubit.depolarize(0.2),
		vqubit.bit_flip(0.2),
		vqubit.phase_flip(0.2),
		vqubit.bit_phase_flip(0.2),
		vqubit.amplitude_damp(0.2),
		vqubit.phase_damp(0.2),
		vqubit.generalized_amplitude_damp(0.2, 0.3),
		vqubit.reset,
	]
	for channel in channels:
		circuit = vqubit.Circuit(vqubit.H(q0), channel(q0), vqubit.measure(q0, key="m"))
		result = vqubit.Simulator(seed=1).run(circuit, repetitions=50)
		assert len(result.measurements["m"]) == 50
		state = vqubit.Simulator(seed=1).simulate(vqubit.Circuit(vqubit.H(q0), channel(q0))).state_vector()
		assert abs(sum(abs(z) ** 2 for z in state) - 1) < 1e-9, channel


def test_a_reset_forgets_the_state() -> None:
	q0 = vqubit.LineQubit(0)
	# ``X`` then ``reset`` then measure must always read zero.
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.reset(q0), vqubit.measure(q0, key="m"))
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=20)
	assert all(row == [0] for row in result.measurements["m"])


def test_a_channel_with_non_diagonal_effects() -> None:
	# With probability 1/2 nothing happens, 1/4 resets to |0> and 1/4 to |1>,
	# which fully randomises the measured qubit. ``K0†K0 + K1†K1 + K2†K2 = I``.
	channel = vqubit.KrausChannel(
		[
			[[math.sqrt(0.5), 0], [0, math.sqrt(0.5)]],
			[[0, math.sqrt(0.5)], [0, 0]],
			[[0, 0], [math.sqrt(0.5), 0]],
		]
	)
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), channel(q0), vqubit.measure(q0, q1, key="m"))
	repetitions = 4000
	trajectories = _histogram(vqubit.Simulator(seed=7).run(circuit, repetitions=repetitions).measurements["m"])
	exact = _histogram(vqubit.DensityMatrixSimulator(seed=7).run(circuit, repetitions=repetitions).measurements["m"])
	assert _gap(trajectories, exact) < 0.05
	# And the analytic answer: q0 is uniform and q1 is perfectly correlated with
	# the branch that does nothing, so the four outcomes come out nearly equal.
	assert set(trajectories) == {0, 1, 2, 3}
	assert all(abs(probability - 0.25) < 0.06 for probability in trajectories.values())


def test_a_channel_with_terminal_measurements_still_samples_per_trajectory() -> None:
	# The "simulate once, sample the distribution" shortcut is only valid for
	# unitary circuits; with a channel each repetition has its own state.
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.bit_flip(0.5)(q0), vqubit.measure(q0, key="m"))
	distribution = _histogram(vqubit.Simulator(seed=11).run(circuit, repetitions=4000).measurements["m"])
	assert abs(distribution.get(0, 0.0) - 0.5) < 0.06
	assert abs(distribution.get(1, 0.0) - 0.5) < 0.06


def test_step_by_step_simulation_applies_channels() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.reset(q0))
	steps = list(vqubit.Simulator(seed=0).simulate_moment_steps(circuit))
	assert [abs(z) for z in steps[-1].state_vector()] == pytest.approx([1.0, 0.0])
