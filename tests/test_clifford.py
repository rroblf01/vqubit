"""Tests for the Clifford (stabilizer) simulator."""

from __future__ import annotations

import pytest

import vqubit


def _aligned(vector: list[complex], reference: list[complex]) -> bool:
	# Compare up to a global phase: the tableau fixes the state only up to one.
	for i in range(len(reference)):
		if abs(reference[i]) > 1e-9:
			phase = vector[i] / reference[i]
			return all(abs(vector[j] - phase * reference[j]) < 1e-9 for j in range(len(reference)))
	return False


def test_bell_state_matches_the_state_vector_simulator() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	clifford = vqubit.CliffordSimulator().simulate(circuit).state_vector()
	exact = vqubit.Simulator().simulate(circuit).state_vector()
	assert _aligned(clifford, exact)


def test_bell_state_stabilizers() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	stabilizers = {str(s) for s in vqubit.CliffordSimulator().simulate(circuit).final_state.stabilizers()}
	assert stabilizers == {"X(0)*X(1)", "Z(0)*Z(1)"}


def test_ghz_state_stabilizers() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.CNOT(q1, q2))
	stabilizers = {str(s) for s in vqubit.CliffordSimulator().simulate(circuit).final_state.stabilizers()}
	assert stabilizers == {"X(0)*X(1)*X(2)", "Z(0)*Z(1)", "Z(1)*Z(2)"}


def test_cz_and_swap_are_clifford() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	clifford = vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.H(q0), vqubit.CZ(q0, q1))).state_vector()
	exact = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.H(q0), vqubit.CZ(q0, q1))).state_vector()
	assert _aligned(clifford, exact)

	swap_clifford = (
		vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.X(q0), vqubit.SWAP(q0, q1))).state_vector()
	)
	swap_exact = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.X(q0), vqubit.SWAP(q0, q1))).state_vector()
	assert _aligned(swap_clifford, swap_exact)


def test_deterministic_measurement() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	result = vqubit.CliffordSimulator().run(circuit, repetitions=10)
	assert result.histogram(key="m") == {1: 10}


def test_bell_measurement_only_produces_00_and_11() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="m"))
	result = vqubit.CliffordSimulator(seed=3).run(circuit, repetitions=300)
	assert set(result.histogram(key="m")).issubset({0, 3})


def test_simulate_rejects_measurements() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.measure(q0)))


def test_non_clifford_gate_raises() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError):
		vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.T(q0)))


def test_repetitions_must_be_positive() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.CliffordSimulator().run(vqubit.Circuit(vqubit.measure(q0)), repetitions=0)


def test_simulate_moment_steps_yields_the_state_after_each_moment() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.X(q1))
	steps = list(vqubit.CliffordSimulator().simulate_moment_steps(circuit))
	assert len(steps) == 3
	# Each step is the state of the circuit up to that moment.
	for index, step in enumerate(steps):
		prefix = vqubit.Circuit(list(circuit.moments)[: index + 1])
		expected = vqubit.CliffordSimulator().simulate(prefix, qubit_order=(q0, q1)).final_state
		assert all(
			abs(got - want) < 1e-9 for got, want in zip(step.state_vector(), expected.state_vector(), strict=True)
		)
		assert step.clifford_state.stabilizers() == expected.stabilizers()
	# The Bell moment is stabilized by XX and ZZ.
	assert sorted(str(stabilizer) for stabilizer in steps[1].clifford_state.stabilizers()) == [
		"X(0)*X(1)",
		"Z(0)*Z(1)",
	]
	with pytest.raises(ValueError, match="with measurements"):
		list(vqubit.CliffordSimulator().simulate_moment_steps(vqubit.Circuit(vqubit.measure(q0, key="m"))))
