"""The Cliffords the native tableau has no opcode for.

Cirq treats ``X ** 0.5``, ``Y ** 0.5``, ``iSWAP``, ``XX``, ``YY`` and ``ZZ`` as
stabilizer operations, so the Clifford simulator does too, by rewriting each one
as a short word over the nine native instructions. The words are only correct up
to a global phase — which is all a tableau can see — so the tests compare the
state the tableau produces with the state the state-vector simulator produces,
and, when Cirq is installed, with Cirq's stabilizer simulator.
"""

from __future__ import annotations

import random

import pytest

import vqubit


def _phase_fit(candidate: list[complex], reference: list[complex]) -> float:
	"""The largest deviation between the two vectors after fitting a phase."""
	phase = 1 + 0j
	for value, goal in zip(candidate, reference, strict=True):
		if abs(goal) > 1e-9:
			phase = value / goal
			break
	return max(abs(value - phase * goal) for value, goal in zip(candidate, reference, strict=True))


def _clifford_gates() -> list[tuple[str, vqubit.Gate]]:
	q0, _q1 = vqubit.LineQubit.range(2)
	_ = q0
	return [
		("X**0.5", vqubit.X**0.5),
		("X**1.5", vqubit.X**1.5),
		("X**-0.5", vqubit.X**-0.5),
		("Y**0.5", vqubit.Y**0.5),
		("Y**1.5", vqubit.Y**1.5),
		("XX**0.5", vqubit.XX**0.5),
		("XX", vqubit.XX),
		("XX**1.5", vqubit.XX**1.5),
		("YY**0.5", vqubit.YY**0.5),
		("YY", vqubit.YY),
		("YY**1.5", vqubit.YY**1.5),
		("ZZ**0.5", vqubit.ZZ**0.5),
		("ZZ", vqubit.ZZ),
		("ZZ**1.5", vqubit.ZZ**1.5),
		("ISWAP", vqubit.ISWAP),
		("ISWAP**-1", vqubit.ISWAP**-1),
		("ISWAP**2", vqubit.ISWAP**2),
	]


@pytest.mark.parametrize(("name", "gate"), _clifford_gates(), ids=[name for name, _ in _clifford_gates()])
def test_a_decomposed_clifford_matches_the_state_vector_simulator(name: str, gate: vqubit.Gate) -> None:
	qubits = vqubit.LineQubit.range(gate.num_qubits)
	circuit = vqubit.Circuit(vqubit.H(qubits[0]), gate(*qubits))
	assert vqubit.has_stabilizer_effect(gate), name
	want = list(vqubit.Simulator().simulate(circuit).state_vector())
	got = list(vqubit.CliffordSimulator().simulate(circuit).state_vector())
	assert _phase_fit(got, want) < 1e-9, name


@pytest.mark.parametrize(("name", "gate"), _clifford_gates(), ids=[name for name, _ in _clifford_gates()])
def test_the_simulator_reports_them_as_supported(name: str, gate: vqubit.Gate) -> None:
	qubits = vqubit.LineQubit.range(gate.num_qubits)
	operation = gate(*qubits)
	assert vqubit.CliffordSimulator().is_supported_operation(operation), name
	# And the step machinery applies the same word.
	steps = list(vqubit.CliffordSimulator().simulate_moment_steps(vqubit.Circuit(operation)))
	matrix = operation.gate._unitary_()
	assert matrix is not None
	first_column = [row[0] for row in matrix]
	assert _phase_fit(list(steps[-1].state.state_vector()), first_column) < 1e-9


@pytest.mark.parametrize(
	"gate",
	[
		vqubit.X**0.25,
		vqubit.Y**0.25,
		vqubit.XX**0.25,
		vqubit.YY**0.25,
		vqubit.ZZ**0.25,
		vqubit.ISWAP**0.5,
		vqubit.T,
	],
	ids=["X**0.25", "Y**0.25", "XX**0.25", "YY**0.25", "ZZ**0.25", "ISWAP**0.5", "T"],
)
def test_the_rest_is_still_outside_the_stabilizer_formalism(gate: vqubit.Gate) -> None:
	qubits = vqubit.LineQubit.range(gate.num_qubits)
	assert not vqubit.has_stabilizer_effect(gate)
	assert not vqubit.CliffordSimulator().is_supported_operation(gate(*qubits))
	with pytest.raises(TypeError, match="doesn't support"):
		vqubit.CliffordSimulator().simulate(vqubit.Circuit(gate(*qubits)))


def test_random_circuits_over_the_new_gates_agree_with_the_state_vector() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	pool = [
		vqubit.X**0.5,
		vqubit.Y**0.5,
		vqubit.X**1.5,
		vqubit.Y**1.5,
		vqubit.XX**0.5,
		vqubit.YY**0.5,
		vqubit.ZZ**0.5,
		vqubit.XX,
		vqubit.YY,
		vqubit.ZZ,
		vqubit.ISWAP,
		vqubit.ISWAP**2,
		vqubit.ISWAP**-1,
		vqubit.H,
		vqubit.S,
		vqubit.CNOT,
	]
	generator = random.Random(20261021)
	for _ in range(40):
		operations = []
		for _ in range(8):
			gate = generator.choice(pool)
			operations.append(gate(*generator.sample([q0, q1], gate.num_qubits)))
		circuit = vqubit.Circuit(operations)
		want = list(vqubit.Simulator().simulate(circuit).state_vector())
		got = list(vqubit.CliffordSimulator().simulate(circuit).state_vector())
		assert _phase_fit(got, want) < 1e-9


def test_measurements_after_a_decomposed_clifford_are_correlated() -> None:
	# ``XX ** 0.5`` then a measurement of both qubits: the two bits must agree,
	# which they only do if the word really is the entangled gate.
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(
		vqubit.H(q0),
		vqubit.CNOT(q0, q1),
		(vqubit.YY**0.5)(q0, q1),
		vqubit.measure(q0, key="a"),
		vqubit.measure(q1, key="b"),
	)
	result = vqubit.CliffordSimulator(seed=3).run(circuit, repetitions=20)
	first = [row[0] for row in result.measurements["a"]]
	second = [row[0] for row in result.measurements["b"]]
	assert first == second


def test_cirq_agrees_on_the_stabilizer_effect_and_the_state() -> None:
	cirq = pytest.importorskip("cirq", exc_type=ImportError)
	numpy = pytest.importorskip("numpy", exc_type=ImportError)

	targets = {
		"X**0.5": cirq.X**0.5,
		"Y**0.5": cirq.Y**0.5,
		"XX**0.5": cirq.XX**0.5,
		"YY**0.5": cirq.YY**0.5,
		"ZZ**0.5": cirq.ZZ**0.5,
		"XX": cirq.XX,
		"YY": cirq.YY,
		"ZZ": cirq.ZZ,
		"ISWAP": cirq.ISWAP,
	}
	mine = {
		"X**0.5": vqubit.X**0.5,
		"Y**0.5": vqubit.Y**0.5,
		"XX**0.5": vqubit.XX**0.5,
		"YY**0.5": vqubit.YY**0.5,
		"ZZ**0.5": vqubit.ZZ**0.5,
		"XX": vqubit.XX,
		"YY": vqubit.YY,
		"ZZ": vqubit.ZZ,
		"ISWAP": vqubit.ISWAP,
	}
	cirq_qubits = cirq.LineQubit.range(2)
	my_qubits = vqubit.LineQubit.range(2)
	for name, cirq_gate in targets.items():
		assert cirq.has_stabilizer_effect(cirq_gate) == vqubit.has_stabilizer_effect(mine[name]), name
		count = cirq_gate.num_qubits()
		circuit_c = cirq.Circuit([cirq.H(cirq_qubits[0]), cirq_gate.on(*cirq_qubits[:count])])
		circuit_v = vqubit.Circuit([vqubit.H(my_qubits[0]), mine[name](*my_qubits[:count])])
		want = numpy.asarray(cirq.CliffordSimulator().simulate(circuit_c).final_state.state_vector())
		got = vqubit.CliffordSimulator().simulate(circuit_v).state_vector()
		assert _phase_fit(list(got), list(want)) < 1e-6, name
