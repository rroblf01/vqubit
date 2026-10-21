"""Tests for the step-by-step simulation API."""

from __future__ import annotations

import math

import pytest

import vqubit
import vqubit_native as native
from tests.helpers import assert_matrix_close, assert_vector_close


def test_simulators_implement_protocols() -> None:
	assert isinstance(vqubit.Simulator(), vqubit.SimulatesSamples)
	assert isinstance(vqubit.Simulator(), vqubit.SimulatesFinalState)
	assert isinstance(vqubit.DensityMatrixSimulator(), vqubit.SimulatesSamples)
	assert issubclass(vqubit.Simulator, vqubit.SimulatorBase)


def test_state_vector_step_result() -> None:
	qubits = vqubit.LineQubit.range(2)
	result = vqubit.StateVectorStepResult(qubits=qubits, measurements={}, state_vector=[1, 0, 0, 0])
	assert result.qubits == tuple(qubits)
	assert result.dirac_notation() == "1.000|00⟩"
	assert result.sample([qubits[0]], repetitions=3, seed=0) == [0, 0, 0]
	assert_matrix_close(result.density_matrix(), [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])


def test_state_vector_simulation_state() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	state = vqubit.StateVectorSimulationState(qubits=(q0, q1))
	state.apply_operation(vqubit.H(q0))
	assert abs(state.state_vector()[0]) == pytest.approx(1 / math.sqrt(2))
	state.apply_operation(vqubit.CNOT(q0, q1))
	assert_vector_close(state.state_vector(), [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)])


def test_density_matrix_simulation_state() -> None:
	q0 = vqubit.LineQubit(0)
	state = vqubit.DensityMatrixSimulationState(qubits=(q0,))
	state.apply_operation(vqubit.X(q0))
	assert_matrix_close(state.density_matrix(), [[0, 0], [0, 1]])
	state.apply_operation(vqubit.bit_flip(0.0)(q0))
	assert_matrix_close(state.density_matrix(), [[0, 0], [0, 1]])


def test_clifford_tableau_simulation_state() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	tableau = native.NativeTableau(2)
	state = vqubit.CliffordTableauSimulationState(tableau=tableau, qubits=(q0, q1))
	state.apply_operation(vqubit.H(q0))
	state.apply_operation(vqubit.CNOT(q0, q1))
	assert_vector_close(state.state_vector(), [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)])


def test_stabilizer_state_ch_form() -> None:
	state = vqubit.StabilizerStateChForm([1, 0], 1)
	assert state.num_qubits == 1
	assert state.state_vector() == [1 + 0j, 0j]
	assert_matrix_close(state.density_matrix(), [[1, 0], [0, 0]])


def test_state_vector_trial_result() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.StateVectorTrialResult(None, {}, qubits=(q0,), state_vector=[1, 0])
	assert result.params is None
	assert result.dirac_notation() == "1.000|0⟩"


def test_clifford_trial_result() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	final_state = vqubit.CliffordSimulator().simulate(circuit).final_state
	result = vqubit.CliffordTrialResult(None, {}, qubits=(q0, q1), final_state=final_state)
	assert result.final_state is final_state


def test_stabilizer_sampler() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	result = vqubit.StabilizerSampler(seed=0).run(circuit, repetitions=5)
	assert result.measurements["m"] == [[1]] * 5


def test_measure_and_sample_functions() -> None:
	bell = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	outcomes, collapsed = vqubit.measure_state_vector(bell, [0], seed=0)
	assert all(outcome in (0, 1) for outcome in outcomes)
	assert abs(sum(abs(z) ** 2 for z in collapsed) - 1) < 1e-9
	outcomes = vqubit.sample_state_vector(bell, [0, 1], repetitions=10, seed=0)
	assert all(o in (0, 3) for o in outcomes)
	rho = vqubit.density_matrix_from_state_vector(bell)
	rho_outcomes, collapsed_rho = vqubit.measure_density_matrix(rho, [0], seed=0)
	assert all(outcome in (0, 1) for outcome in rho_outcomes)
	assert len(collapsed_rho) == 4
	assert len(vqubit.sample_density_matrix(rho, [0], repetitions=5, seed=0)) == 5


def test_final_state_helpers() -> None:
	q0 = vqubit.LineQubit(0)
	simulation = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.X(q0)))
	assert vqubit.final_state_vector(simulation) == [0j, 1 + 0j]
	density = vqubit.DensityMatrixSimulator().simulate(vqubit.Circuit(vqubit.X(q0)))
	assert_matrix_close(vqubit.final_density_matrix(density), [[0, 0], [0, 1]])


def test_sample_sweep() -> None:
	q0 = vqubit.LineQubit(0)
	parameter = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.rx(parameter)(q0), vqubit.measure(q0, key="m"))
	results = vqubit.sample_sweep(
		circuit, vqubit.Points([{parameter: 0.0}, {parameter: 3.141592653589793}]), repetitions=4, seed=0
	)
	assert len(results) == 2
	assert results[0].measurements["m"] == [[0]] * 4
