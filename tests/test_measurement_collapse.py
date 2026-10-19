"""Regression tests: a measurement must collapse the state.

The state-vector simulator runs its terminal measurements in a single pass over
the final state, reusing it for every repetition. That pass used to sample each
measurement independently from its *own* marginal, which decorrelated separate
``measure`` calls: the two qubits of a Bell state measured one by one came out
uncorrelated, and measuring the same qubit twice could report two different bits.
The joint distribution over the union of the measured qubits is what collapses
correctly, so every assertion below is an invariant that must hold on every
shot, not a statistical threshold.
"""

from __future__ import annotations

import vqubit

REPETITIONS = 200


def _bits(result: vqubit.Result, key: str) -> list[int]:
	return [int(row[0]) for row in result.measurements[key]]


def test_separate_measurements_of_a_bell_pair_agree() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(
		vqubit.H(q0),
		vqubit.CNOT(q0, q1),
		vqubit.measure(q0, key="m"),
		vqubit.measure(q1, key="n"),
	)
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)
	assert _bits(result, "m") == _bits(result, "n")


def test_measurement_order_does_not_matter() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(
		vqubit.H(q0),
		vqubit.CNOT(q0, q1),
		vqubit.measure(q1, key="n"),
		vqubit.measure(q0, key="m"),
	)
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)
	assert _bits(result, "m") == _bits(result, "n")


def test_measuring_the_same_qubit_twice_agrees() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="a"), vqubit.measure(q0, key="b"))
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)
	assert _bits(result, "a") == _bits(result, "b")


def test_separate_measurements_of_a_ghz_state_agree() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(
		vqubit.H(q0),
		vqubit.CNOT(q0, q1),
		vqubit.CNOT(q1, q2),
		vqubit.measure(q0, key="a"),
		vqubit.measure(q1, key="b"),
		vqubit.measure(q2, key="c"),
	)
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)
	assert _bits(result, "a") == _bits(result, "b") == _bits(result, "c")


def test_marginals_are_still_random() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="a"), vqubit.measure(q0, key="b"))
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)
	column = _bits(result, "a")
	assert set(column) == {0, 1}


def test_joint_measurement_stays_correlated() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="both"))
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)
	assert {tuple(row) for row in result.measurements["both"]} <= {(0, 0), (1, 1)}


def test_interleaved_measurements_follow_the_collapse_too() -> None:
	# An operation after the measurements takes the circuit off the terminal fast
	# path; both paths must agree.
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(
		vqubit.H(q0),
		vqubit.CNOT(q0, q1),
		vqubit.measure(q0, key="m"),
		vqubit.measure(q1, key="n"),
		vqubit.I(q0),
	)
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)
	assert _bits(result, "m") == _bits(result, "n")


def test_pure_python_path_collapses_without_numpy(monkeypatch) -> None:
	import vqubit._sim as sim

	monkeypatch.setattr(sim, "HAS_NUMPY", False)
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(
		vqubit.H(q0),
		vqubit.CNOT(q0, q1),
		vqubit.measure(q0, key="m"),
		vqubit.measure(q1, key="n"),
	)
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)
	assert _bits(result, "m") == _bits(result, "n")


def test_a_circuit_without_measurements_gives_an_empty_result() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.Simulator().run(vqubit.Circuit(vqubit.H(q0)), repetitions=4)
	assert result.measurements == {}
	assert len(result) == 0
