"""Tests for the sweep, step-by-step and amplitude entry points of the simulators.

These are the calls Cirq users reach for when they port code: ``run_sweep`` and
``simulate_sweep`` over a :class:`~vqubit.Sweep`, ``simulate_moment_steps`` to
watch the state evolve, ``compute_amplitudes`` for selected basis states and
``run_batch`` for a batch of circuits. Each simulator answers all of them through
the same mixin, so the tests run over the three of them.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

import vqubit
from tests.helpers import assert_matrix_close

SIMULATORS = [vqubit.Simulator, vqubit.DensityMatrixSimulator, vqubit.CliffordSimulator]


@pytest.mark.parametrize("simulator_class", SIMULATORS)
def test_run_sweep_walks_the_sweep_and_tags_each_result(simulator_class: type) -> None:
	q0 = vqubit.LineQubit(0)
	t: Any = vqubit.Parameter("t")
	circuit = vqubit.Circuit((vqubit.X**t)(q0), vqubit.measure(q0, key="m"))
	# ``X ** 0`` is the identity and ``X ** 1`` is ``X``, so the point of the
	# sweep decides the deterministic outcome.
	points = [{"t": 0.0}, {"t": 1.0}]
	results = simulator_class(seed=0).run_sweep(circuit, vqubit.Points(points), repetitions=3)
	assert [dict(result.params)["t"] for result in results] == [0.0, 1.0]
	assert results[0].histogram(key="m") == {0: 3}
	assert results[1].histogram(key="m") == {1: 3}


def test_run_sweep_iter_yields_lazily() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	iterator = vqubit.Simulator(seed=0).run_sweep_iter(circuit, vqubit.Points([{"t": 1.0}]), repetitions=2)
	assert hasattr(iterator, "__next__")
	assert len(list(iterator)) == 1


@pytest.mark.parametrize("simulator_class", [vqubit.Simulator, vqubit.DensityMatrixSimulator])
def test_simulate_sweep_returns_the_state_of_every_point(simulator_class: type) -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.rx(t)(q0))
	results = simulator_class().simulate_sweep(circuit, vqubit.Points([{"t": 0.0}, {"t": math.pi}]))
	assert len(results) == 2
	first = results[0].state_vector() if simulator_class is vqubit.Simulator else results[0].density_matrix()
	assert len(first) == 2
	assert [dict(result.params)["t"] for result in results] == [0.0, math.pi]


def test_simulate_sweep_iter_is_lazy() -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.rx(t)(q0))
	iterator = vqubit.Simulator().simulate_sweep_iter(circuit, vqubit.Points([{"t": 0.5}]))
	assert hasattr(iterator, "__next__")
	assert len(list(iterator)) == 1


def test_simulate_moment_steps_yields_one_state_per_moment() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.S(q1))
	steps = list(vqubit.Simulator().simulate_moment_steps(circuit))
	assert len(steps) == 3
	# The first step is the state after the Hadamard, so |00> and |11> have an
	# amplitude of one half and the others none.
	assert_matrix_close(vqubit.unitary(vqubit.S(q1).gate), vqubit.unitary(vqubit.S))
	first = steps[0].state_vector()
	assert abs(first[0]) > 0.5 and abs(first[2]) > 0.5
	assert steps[0].qubits == (q0, q1)
	# Later snapshots are independent copies, so the earlier ones do not move.
	assert abs(steps[0].state_vector()[3]) < 1e-12
	assert abs(steps[-1].state_vector()[3]) > 0.5


def test_density_matrix_moment_steps() -> None:
	q0 = vqubit.LineQubit(0)
	steps = list(vqubit.DensityMatrixSimulator().simulate_moment_steps(vqubit.Circuit(vqubit.X(q0))))
	assert len(steps) == 1
	matrix = steps[0].density_matrix()
	assert abs(matrix[0][0]) < 1e-12
	assert abs(matrix[1][1] - 1) < 1e-12


@pytest.mark.parametrize("simulator_class", [vqubit.Simulator, vqubit.DensityMatrixSimulator])
def test_simulate_moment_steps_refuses_measurements(simulator_class: type) -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	with pytest.raises(ValueError):
		list(simulator_class().simulate_moment_steps(circuit))


def test_compute_amplitudes_picks_the_asked_for_basis_states() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	amplitudes = vqubit.Simulator().compute_amplitudes(circuit, [0, 1, 3])
	assert abs(amplitudes[0] - 2**-0.5) < 1e-9
	assert abs(amplitudes[1]) < 1e-12
	assert abs(amplitudes[2] - 2**-0.5) < 1e-9


def test_compute_amplitudes_resolves_parameters() -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.rx(t)(q0))
	amplitudes = vqubit.Simulator().compute_amplitudes(circuit, [0, 1], vqubit.ParamResolver({"t": math.pi}))
	assert abs(amplitudes[0]) < 1e-12


def test_run_batch_takes_a_circuit_list_with_its_own_repetitions() -> None:
	q0 = vqubit.LineQubit(0)
	one = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	zero = vqubit.Circuit(vqubit.I(q0), vqubit.measure(q0, key="m"))
	batches = vqubit.Simulator(seed=0).run_batch([zero, one], repetitions=[2, 5])
	assert len(batches) == 2 and all(len(batch) == 1 for batch in batches)
	assert batches[0][0].histogram(key="m") == {0: 2}
	assert batches[1][0].histogram(key="m") == {1: 5}


def test_run_batch_rejects_mismatched_lengths() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="m"))
	with pytest.raises(ValueError):
		vqubit.Simulator().run_batch([circuit, circuit], repetitions=[1])


def test_run_batch_runs_an_empty_params_list() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	batches = vqubit.Simulator(seed=0).run_batch([circuit], params_list=[vqubit.Points([{"t": 1.0}])])
	assert batches[0][0].histogram(key="m") == {1: 1}


def test_param_resolver_is_the_second_positional_argument() -> None:
	# ``simulator.simulate(circuit, resolver)`` is the Cirq idiom.
	q0 = vqubit.LineQubit(0)
	t: Any = vqubit.Parameter("t")
	circuit = vqubit.Circuit((vqubit.X**t)(q0))
	assert_matrix_close(
		vqubit.unitary(vqubit.X),
		[[0, 1], [1, 0]],
	)
	result = vqubit.Simulator().simulate(circuit, vqubit.ParamResolver({"t": 1}))
	assert abs(result.state_vector()[1] - 1) < 1e-12
	counts = vqubit.Simulator(seed=0).run(circuit, vqubit.ParamResolver({"t": 2}), 3)
	assert len(counts.measurements) == 0
	assert counts.params["t"] == 2


def test_a_circuit_with_symbols_needs_a_resolver() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.rx(vqubit.Parameter("t"))(q0))
	with pytest.raises(ValueError):
		vqubit.Simulator().simulate(circuit)
	with pytest.raises(ValueError):
		vqubit.DensityMatrixSimulator().simulate(circuit)
	with pytest.raises(ValueError):
		vqubit.CliffordSimulator().simulate(circuit)


def test_clifford_is_supported_operation() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	simulator = vqubit.CliffordSimulator()
	assert simulator.is_supported_operation(vqubit.H(q0))
	assert simulator.is_supported_operation(vqubit.CZ(q0, q1))
	assert simulator.is_supported_operation(vqubit.measure(q0, key="m"))
	assert not simulator.is_supported_operation(vqubit.T(q0))
	assert not simulator.is_supported_operation(vqubit.rx(0.3)(q0))


def test_the_clifford_simulator_refuses_what_it_cannot_run() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError):
		vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.T(q0)))
	with pytest.raises(TypeError):
		vqubit.CliffordSimulator().run(vqubit.Circuit(vqubit.T(q0), vqubit.measure(q0, key="m")))


def test_has_stabilizer_effect_answers_instead_of_raising() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.has_stabilizer_effect(vqubit.H(q0))
	assert vqubit.has_stabilizer_effect(vqubit.I(q0))
	assert vqubit.has_stabilizer_effect(vqubit.measure(q0, key="m"))
	assert not vqubit.has_stabilizer_effect(vqubit.T(q0))
	assert not vqubit.has_stabilizer_effect(vqubit.depolarize(0.1)(q0))


def test_results_carry_the_cirq_attributes() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.X(q0)))
	assert result.qubit_map == {q0: 0}
	assert result.measurements == {}
	assert len(result.final_state_vector) == 2
	density = vqubit.DensityMatrixSimulator().simulate(vqubit.Circuit(vqubit.X(q0)))
	assert density.qubit_map == {q0: 0}
	assert density.measurements == {}
