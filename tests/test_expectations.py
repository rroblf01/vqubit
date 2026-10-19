"""Tests for expectation values: the Pauli maths and the simulator entry points.

Expectation values are what variational algorithms ask for, and they come in two
flavours, as in Cirq: exact ones off a simulated final state, and estimates from
sampling a circuit's measurements.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

import vqubit
from tests.helpers import assert_matrix_close

SQRT_HALF = math.sqrt(0.5)


def bell_state() -> Any:
	q0, q1 = vqubit.LineQubit.range(2)
	return vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))


def test_pauli_expectation_from_state_vector() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	state = vqubit.Simulator().simulate(bell_state()).state_vector()
	qubit_map = {q0: 0, q1: 1}
	# The Bell state is stabilized by ZZ and XX, and has no Z (or X) on one qubit.
	assert (
		abs(vqubit.PauliString({q0: vqubit.Z, q1: vqubit.Z}).expectation_from_state_vector(state, qubit_map) - 1)
		< 1e-12
	)
	assert (
		abs(vqubit.PauliString({q0: vqubit.X, q1: vqubit.X}).expectation_from_state_vector(state, qubit_map) - 1)
		< 1e-12
	)
	assert abs(vqubit.PauliString({q0: vqubit.Z}).expectation_from_state_vector(state, qubit_map)) < 1e-12
	assert (
		abs(
			vqubit.PauliString({q0: vqubit.Z, q1: vqubit.Z}, coefficient=2.5).expectation_from_state_vector(
				state, qubit_map
			)
			- 2.5
		)
		< 1e-12
	)


def test_pauli_expectation_reads_the_state_in_any_qubit_order() -> None:
	# ``|+>`` on q0 and ``|0>`` on q1: the qubit map says which index of the state
	# each qubit holds, so swapping it swaps which qubit is read.
	q0, q1 = vqubit.LineQubit.range(2)
	result = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.H(q0)), qubit_order=[q1, q0])
	state = result.state_vector()
	qubit_map = result.qubit_map
	assert abs(vqubit.PauliString({q0: vqubit.X}).expectation_from_state_vector(state, qubit_map) - 1) < 1e-12
	assert abs(vqubit.PauliString({q1: vqubit.X}).expectation_from_state_vector(state, qubit_map)) < 1e-12
	swapped = {read: qubit_map[held] for read, held in zip((q0, q1), (q1, q0), strict=True)}
	assert abs(vqubit.PauliString({q0: vqubit.X}).expectation_from_state_vector(state, swapped)) < 1e-12
	assert abs(vqubit.PauliString({q1: vqubit.X}).expectation_from_state_vector(state, swapped) - 1) < 1e-12


def test_pauli_expectation_from_density_matrix_matches_the_state_vector() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.rz(0.4)(q1))
	state = vqubit.Simulator().simulate(circuit).state_vector()
	density = vqubit.DensityMatrixSimulator().simulate(circuit).density_matrix()
	observable = vqubit.PauliString({q0: vqubit.Y, q1: vqubit.Z})
	assert (
		abs(
			observable.expectation_from_state_vector(state, {q0: 0, q1: 1})
			- observable.expectation_from_density_matrix(density, {q0: 0, q1: 1})
		)
		< 1e-9
	)


def test_pauli_sum_expectation_is_the_sum_of_its_terms() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	state = vqubit.Simulator().simulate(circuit).state_vector()
	observable = 0.5 * vqubit.PauliString({q0: vqubit.Z}) + 0.3 * vqubit.PauliString({q0: vqubit.Z, q1: vqubit.Z}) - 1.2
	assert abs(observable.expectation_from_state_vector(state, {q0: 0, q1: 1}) - (0.3 - 1.2)) < 1e-12
	density = vqubit.DensityMatrixSimulator().simulate(circuit).density_matrix()
	assert abs(observable.expectation_from_density_matrix(density, {q0: 0, q1: 1}) - (0.3 - 1.2)) < 1e-12


def test_a_pauli_string_takes_a_pauli_gate_or_letter() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.PauliString({q0: vqubit.Z}) == vqubit.PauliString({q0: vqubit.Pauli.Z})
	assert vqubit.PauliString({q0: "x"}) == vqubit.PauliString({q0: vqubit.X})
	assert vqubit.PauliString({q0: vqubit.X**3}) == vqubit.PauliString({q0: vqubit.X})
	with pytest.raises(TypeError):
		vqubit.PauliString({q0: vqubit.H})


def test_expectation_refuses_a_non_hermitian_observable() -> None:
	q0 = vqubit.LineQubit(0)
	state = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.X(q0))).state_vector()
	with pytest.raises(NotImplementedError):
		vqubit.PauliString({q0: vqubit.X}, coefficient=1j).expectation_from_state_vector(state, {q0: 0})


def test_expectation_refuses_an_unnormalized_state_or_a_bad_map() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	observable = vqubit.PauliString({q0: vqubit.Z, q1: vqubit.Z})
	with pytest.raises(ValueError):
		observable.expectation_from_state_vector([2.0, 0j, 0j, 0j], {q0: 0, q1: 1})
	with pytest.raises(ValueError):
		observable.expectation_from_state_vector([1.0, 0j, 0j, 0j], {q0: 0})
	with pytest.raises(ValueError):
		vqubit.PauliString({q0: vqubit.Z}).expectation_from_state_vector([1.0, 0j, 0j, 0j], {q0: 7})
	# The checks are optional, which is what the simulator methods use; a state of
	# norm 2 then gives four times the expectation value.
	assert (
		abs(observable.expectation_from_state_vector([2.0, 0j, 0j, 0j], {q0: 0, q1: 1}, check_preconditions=False) - 4)
		< 1e-12
	)


@pytest.mark.parametrize("simulator_class", [vqubit.Simulator, vqubit.DensityMatrixSimulator])
def test_simulate_expectation_values(simulator_class: type) -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	observables = [vqubit.PauliString({q0: vqubit.Z, q1: vqubit.Z}), vqubit.PauliString({q0: vqubit.Z})]
	values = simulator_class().simulate_expectation_values(bell_state(), observables)
	assert abs(values[0] - 1) < 1e-9
	assert abs(values[1]) < 1e-9
	# A single observable needs no list, as in Cirq.
	single = simulator_class().simulate_expectation_values(bell_state(), observables[0])
	assert abs(single[0] - 1) < 1e-9


def test_simulate_expectation_values_sweep() -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	circuit = vqubit.Circuit(vqubit.rx(t)(q0))
	observable = vqubit.PauliString({q0: vqubit.Z})
	values = vqubit.Simulator().simulate_expectation_values_sweep(
		circuit, observable, vqubit.Points([{"t": 0.0}, {"t": math.pi}])
	)
	assert abs(values[0][0] - 1) < 1e-9
	assert abs(values[1][0] + 1) < 1e-9
	iterator = vqubit.Simulator().simulate_expectation_values_sweep_iter(
		circuit, observable, vqubit.Points([{"t": 0.0}])
	)
	assert [list(row) for row in iterator] == [[values[0][0]]]


def test_simulate_expectation_values_refuses_terminal_measurements() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	observable = vqubit.PauliString({q0: vqubit.Z})
	with pytest.raises(ValueError):
		vqubit.Simulator().simulate_expectation_values(circuit, observable)
	# ``permit_terminal_measurements`` lets the collapse happen, as Cirq's does,
	# so the value is the one of the collapsed state (+1 or -1), not of |+>.
	values = vqubit.Simulator(seed=0).simulate_expectation_values(
		circuit, observable, permit_terminal_measurements=True
	)
	assert abs(abs(values[0]) - 1) < 1e-9


def test_sample_expectation_values_approaches_the_exact_value() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.rz(0.7)(q1))
	observables = [vqubit.PauliString({q0: vqubit.X, q1: vqubit.X}), vqubit.PauliString({q0: vqubit.Z})]
	exact = vqubit.Simulator().simulate_expectation_values(circuit, observables)
	sampled = vqubit.Simulator(seed=5).sample_expectation_values(circuit, observables, num_samples=4000)[0]
	for estimate, value in zip(sampled, exact, strict=True):
		assert abs(estimate - value) < 0.1
	# The Clifford simulator samples too.
	clifford_circuit = bell_state()
	estimate = vqubit.CliffordSimulator(seed=5).sample_expectation_values(
		clifford_circuit, vqubit.PauliString({q0: vqubit.Z, q1: vqubit.Z}), num_samples=2000
	)
	assert abs(estimate[0][0] - 1) < 1e-9


def test_sample_expectation_values_validates_its_arguments() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	observable = vqubit.PauliString({q0: vqubit.Z})
	with pytest.raises(ValueError):
		vqubit.Simulator().sample_expectation_values(circuit, observable, num_samples=0)
	with pytest.raises(ValueError):
		vqubit.Simulator().sample_expectation_values(circuit, [], num_samples=10)
	measured = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	with pytest.raises(ValueError):
		vqubit.Simulator().sample_expectation_values(measured, observable, num_samples=10)


def test_bloch_vector_of_and_density_matrix_of() -> None:
	q0 = vqubit.LineQubit(0)
	zero = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.I(q0)))
	plus = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.H(q0)))
	assert_matrix_close([zero.bloch_vector_of(q0)], [[0.0, 0.0, 1.0]])
	assert_matrix_close([plus.bloch_vector_of(q0)], [[1.0, 0.0, 0.0]])
	assert_matrix_close(zero.density_matrix_of([q0]), [[1.0, 0.0], [0.0, 0.0]])
	assert zero.qubit_map == {q0: 0}


def test_the_trial_result_reports_what_it_simulated() -> None:
	q0 = vqubit.LineQubit(0)
	t = vqubit.Parameter("t")
	result = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.rx(t)(q0)), vqubit.ParamResolver({"t": 0.5}))
	assert result.params["t"] == 0.5
	assert result.measurements == {}
	assert len(result.final_state_vector) == 2
