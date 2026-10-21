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


def test_the_state_vector_agrees_with_the_state_vector_simulator() -> None:
	# The tableau's state used to be reconstructed with dense projectors, which
	# was O(2**3n): correct for a handful of qubits and a minute of work (or a
	# dead process) for twelve. A random circuit per size pins both the algebra
	# (including the `i` of a Y row) and the speed.
	import random

	import vqubit as _vqubit

	generator = random.Random(20261021)
	for num_qubits in (2, 3, 4, 5, 6, 8, 10):
		qubits = _vqubit.LineQubit.range(num_qubits)
		for _ in range(6):
			operations = []
			for _ in range(8):
				kind = generator.choice(["h", "s", "sdg", "x", "y", "z", "cnot", "cz", "swap"])
				if kind in ("cnot", "cz", "swap"):
					first, second = generator.sample(range(num_qubits), 2)
					gate = {"cnot": _vqubit.CNOT, "cz": _vqubit.CZ, "swap": _vqubit.SWAP}[kind]
					operations.append(gate(qubits[first], qubits[second]))
				else:
					gate = {
						"h": _vqubit.H,
						"s": _vqubit.S,
						"sdg": _vqubit.S**-1,
						"x": _vqubit.X,
						"y": _vqubit.Y,
						"z": _vqubit.Z,
					}[kind]
					operations.append(gate(qubits[generator.randrange(num_qubits)]))
			circuit = _vqubit.Circuit(operations)
			exact = list(_vqubit.Simulator().simulate(circuit).state_vector())
			tableau = list(_vqubit.CliffordSimulator().simulate(circuit).state_vector())
			assert _aligned(tableau, exact)


def test_the_state_vector_scales_to_twenty_qubits() -> None:
	import time

	qubits = vqubit.LineQubit.range(20)
	circuit = vqubit.Circuit([vqubit.H(qubit) for qubit in qubits])
	start = time.perf_counter()
	vector = vqubit.CliffordSimulator().simulate(circuit).state_vector()
	elapsed = time.perf_counter() - start
	assert len(vector) == 1 << 20
	assert abs(sum(abs(value) ** 2 for value in vector) - 1) < 1e-9
	# The dense-projector version needed O(2**3n) operations here; anything
	# under a couple of seconds can only be the linear algorithm.
	assert elapsed < 5.0


def test_a_huge_tableau_refuses_to_materialise_its_state_vector() -> None:
	qubits = vqubit.LineQubit.range(30)
	circuit = vqubit.Circuit([vqubit.H(qubit) for qubit in qubits])
	state = vqubit.CliffordSimulator().simulate(circuit).final_state
	with pytest.raises(ValueError, match="amplitudes"):
		state.state_vector()
	# The tableau itself is still available, which is the point of a stabilizer.
	assert len(state.stabilizers()) == 30


def test_the_tableau_takes_y_powers_by_hand() -> None:
	# ``apply_y`` has a branch per half-integer exponent; the simulator only ever
	# reaches the integer one, so the others are exercised here.
	for exponent in (0.5, 1.0, 1.5, 2.0):
		tableau = vqubit.CliffordTableau(1)
		tableau.apply_y(0, exponent)
		# A Y power is a Clifford, so the stabilisers it produces are well formed.
		stabilizers = tableau.stabilizer_pauli_strings()
		assert len(stabilizers) == 1
	with pytest.raises(ValueError, match="half integer"):
		vqubit.CliffordTableau(1).apply_y(0, 0.25)
