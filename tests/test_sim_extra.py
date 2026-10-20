"""Extra behavioural tests for the simulator modules.

These cover the corners of :mod:`vqubit._sim` and :mod:`vqubit._sim_step`
that the other simulator tests leave out: the state and step-result classes,
the mutable simulation states, the pure-Python fallbacks and the error paths.
Every number below is chosen so that it can be checked by hand.
"""

from __future__ import annotations

import asyncio
import math
import random
from typing import Any

import pytest

import vqubit
import vqubit_native as native
from tests.helpers import assert_matrix_close, assert_vector_close

SQRT_HALF = math.sqrt(0.5)


class _NoUnitaryGate(vqubit.Gate):
	"""A gate vqubit knows nothing about: no unitary and no Kraus operators."""

	def _num_qubits_(self) -> int:
		return 1


class _AlwaysOne(random.Random):
	"""A PRNG whose draw is at the very top of ``[0, 1)``."""

	def random(self) -> float:
		return 1.0


def _bell_circuit() -> Any:
	q0, q1 = vqubit.LineQubit.range(2)
	return vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))


def test_state_vector_carries_its_qubits_and_formats_amplitudes() -> None:
	q0 = vqubit.LineQubit(0)
	state = vqubit.StateVector([1 + 0j, 0j], [q0])
	assert state.qubits == (q0,)
	assert state.state_vector() == [1 + 0j, 0j]
	assert str(state) == "1.000|0⟩"
	assert repr(state) == "vqubit.StateVector(1.000|0⟩)"
	# A non-zero imaginary part is printed too, rather than dropped.
	assert vqubit.StateVector([0j, 1j], [q0]).dirac_notation() == "0.000+1.000j|1⟩"


def test_simulate_accepts_a_state_vector_as_initial_state() -> None:
	q0 = vqubit.LineQubit(0)
	initial = vqubit.StateVector([1 + 0j, 0j], [q0])
	result = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.X(q0)), initial_state=initial)
	assert_vector_close(list(result.state_vector()), [0j, 1 + 0j])


def test_initial_state_must_fit_the_number_of_amplitudes() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError, match="out of range"):
		vqubit.Simulator().simulate(vqubit.Circuit(), qubit_order=[q0], initial_state=2)
	with pytest.raises(ValueError, match="length 3"):
		vqubit.Simulator().simulate(vqubit.Circuit(), qubit_order=[q0], initial_state=[1, 0, 0])
	# A plain list is a valid initial state as well.
	result = vqubit.Simulator().simulate(vqubit.Circuit(), qubit_order=[q0], initial_state=[0, 1])
	assert_vector_close(list(result.state_vector()), [0j, 1 + 0j])


def test_simulation_result_describes_its_final_state() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	result = vqubit.Simulator().simulate(_bell_circuit())
	assert result.qubits == (q0, q1)
	assert result.final_state.qubits == (q0, q1)
	assert_vector_close(list(result.final_state_vector), [SQRT_HALF, 0, 0, SQRT_HALF])
	assert str(result) == "0.707|00⟩ + 0.707|11⟩"
	assert repr(result) == "vqubit.SimulationResult(vqubit.StateVector(0.707|00⟩ + 0.707|11⟩))"
	assert result.qubit_map == {q0: 0, q1: 1}


def test_get_state_containing_qubit_returns_the_whole_state() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	result = vqubit.Simulator().simulate(_bell_circuit())
	# vqubit keeps one state for the whole circuit, so asking for either qubit
	# gives that state back (with both qubits in it).
	state = result.get_state_containing_qubit(q1)
	assert state.qubits == (q0, q1)
	assert_vector_close(state.state_vector(), [SQRT_HALF, 0, 0, SQRT_HALF])


def test_result_reports_params_histogram_and_repr() -> None:
	plain = vqubit.Result({"m": [[1]]}, params={"t": 2})
	assert plain.params["t"] == 2
	assert plain.histogram("m") == {1: 1}
	result = vqubit.Result({"a": [[0], [1]], "b": [[1], [0]]})
	assert len(result) == 2
	# Without a key, a result with two measurements cannot fold into one number.
	with pytest.raises(ValueError, match="needs a key"):
		result.histogram()
	assert repr(result) == "vqubit.Result(keys=['a', 'b'])"
	assert result.histogram("a") == {0: 1, 1: 1}


def test_simulate_falls_back_when_the_native_view_is_unavailable(monkeypatch: pytest.MonkeyPatch) -> None:
	import vqubit._sim as sim

	monkeypatch.setattr(sim, "native_complex_view", lambda *args, **kwargs: None)
	result = sim.Simulator().simulate(_bell_circuit())
	assert isinstance(result.state_vector(), list)
	assert_vector_close(result.state_vector(), [SQRT_HALF, 0, 0, SQRT_HALF])


def test_density_simulator_falls_back_without_the_native_view(monkeypatch: pytest.MonkeyPatch) -> None:
	import vqubit._sim as sim

	monkeypatch.setattr(sim, "native_complex_view", lambda *args, **kwargs: None)
	result = sim.DensityMatrixSimulator().simulate(_bell_circuit())
	assert_matrix_close(
		result.density_matrix(),
		[[0.5, 0, 0, 0.5], [0, 0, 0, 0], [0, 0, 0, 0], [0.5, 0, 0, 0.5]],
	)


def test_compute_amplitudes_takes_numpy_bitstrings() -> None:
	np = pytest.importorskip("numpy", exc_type=ImportError)
	amplitudes = vqubit.Simulator().compute_amplitudes(_bell_circuit(), np.array([0, 3]))
	assert_vector_close(amplitudes, [SQRT_HALF, SQRT_HALF])


def test_moment_steps_can_start_from_another_basis_state() -> None:
	q0 = vqubit.LineQubit(0)
	steps = list(vqubit.Simulator().simulate_moment_steps(vqubit.Circuit(vqubit.X(q0)), initial_state=1))
	assert len(steps) == 1
	# X on |1> is |0>.
	assert_vector_close(list(steps[0].state_vector()), [1 + 0j, 0j])


def test_run_samples_from_a_given_initial_state() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="m"))
	result = vqubit.Simulator(seed=0).run(circuit, repetitions=3, initial_state=1)
	assert result.measurements["m"] == [[1]] * 3


def test_density_simulator_accepts_initial_states() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.DensityMatrixSimulator().simulate(vqubit.Circuit(), qubit_order=[q0], initial_state=[0, 1])
	assert_matrix_close(result.density_matrix(), [[0, 0], [0, 1]])
	steps = list(vqubit.DensityMatrixSimulator().simulate_moment_steps(vqubit.Circuit(vqubit.X(q0)), initial_state=1))
	assert_matrix_close(steps[-1].density_matrix(), [[1, 0], [0, 0]])


def test_simulate_rejects_a_non_unitary_operation() -> None:
	q0 = vqubit.LineQubit(0)
	# A reset is a channel: it has Kraus operators but no unitary, so the
	# state-vector simulator must refuse it instead of guessing.
	with pytest.raises(TypeError, match="has no unitary"):
		vqubit.Simulator().simulate(vqubit.Circuit(vqubit.reset(q0)))


def test_run_rejects_a_gate_without_a_unitary() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError, match="has no unitary"):
		vqubit.Simulator().run(vqubit.Circuit(_NoUnitaryGate().on(q0)))


def test_density_simulator_rejects_a_gate_without_unitary_or_kraus() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(_NoUnitaryGate().on(q0))
	with pytest.raises(TypeError, match="neither a unitary nor Kraus"):
		vqubit.DensityMatrixSimulator().simulate(circuit)
	with pytest.raises(TypeError, match="neither a unitary nor Kraus"):
		vqubit.DensityMatrixSimulator().run(circuit)


def test_async_entry_points_wrap_the_synchronous_ones() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	simulator = vqubit.Simulator(seed=0)
	run = asyncio.run(simulator.run_async(circuit, repetitions=2))
	assert run.measurements["m"] == [[1], [1]]
	sweep = asyncio.run(simulator.run_sweep_async(circuit, vqubit.Points([{"t": 1.0}]), 2))
	assert sweep[0].measurements["m"] == [[1], [1]]
	batch = asyncio.run(simulator.run_batch_async([circuit], repetitions=2))
	assert batch[0][0].measurements["m"] == [[1], [1]]


def test_density_matrix_objects_and_trial_results() -> None:
	q0 = vqubit.LineQubit(0)
	matrix = vqubit.DensityMatrix([[1 + 0j, 0j], [0j, 0j]], [q0])
	assert matrix.qubits == (q0,)
	assert_matrix_close(matrix.density_matrix(), [[1, 0], [0, 0]])
	assert matrix.dirac_notation() == "1.000|0⟩"
	assert repr(matrix) == "vqubit.DensityMatrix(1 qubits)"
	result = vqubit.DensityMatrixSimulator().simulate(vqubit.Circuit(vqubit.X(q0)))
	assert result.qubits == (q0,)
	assert repr(result) == "vqubit.DensityMatrixTrialResult(vqubit.DensityMatrix(1 qubits))"


def test_an_identity_observable_is_sampled_as_its_coefficient() -> None:
	q0 = vqubit.LineQubit(0)
	# The empty Pauli string is the identity, so sampling it is exact: no run
	# is needed and the answer is the coefficient itself.
	identity = 2.0 * vqubit.PauliString({})
	values = vqubit.Simulator(seed=0).sample_expectation_values(vqubit.Circuit(vqubit.H(q0)), identity, num_samples=5)
	assert values == [[2.0]]


def test_sampling_helpers_fall_back_to_the_last_outcome() -> None:
	import vqubit._sim as sim
	import vqubit._sim_step as step

	rng = _AlwaysOne()
	assert step._draw([0.5, 0.5], rng) == 1
	assert sim._draw_outcome([0.5, 0.5], rng) == 1


def test_step_result_reports_its_qubit_map_and_repr() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	result = vqubit.StateVectorStepResult(qubits=(q0, q1), measurements={}, state_vector=[1, 0, 0, 0])
	assert result.qubit_map == {q0: 0, q1: 1}
	assert repr(result) == "vqubit.StepResult(2 qubits)"


def test_a_density_matrix_step_has_no_state_vector() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.DensityMatrixStepResult(qubits=(q0,), measurements={}, density_matrix=[[1, 0], [0, 0]])
	assert_matrix_close(result.density_matrix(), [[1, 0], [0, 0]])
	with pytest.raises(ValueError, match="no state vector"):
		result.state_vector()


def test_set_state_vector_replaces_the_step_state() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.StateVectorStepResult(qubits=(q0,), measurements={}, state_vector=[1, 0])
	result.set_state_vector([0, 1])
	assert_vector_close(result.state_vector(), [0, 1])
	assert result.sample([q0], repetitions=3, seed=0) == [1, 1, 1]


def test_sample_measurement_ops_collapses_each_shot() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	step = vqubit.StateVectorStepResult(qubits=(q0, q1), measurements={}, state_vector=[SQRT_HALF, 0, 0, SQRT_HALF])
	# Measuring q0 and then q1 on a Bell pair collapses the first outcome onto
	# the second, so the two columns must agree shot by shot.
	records = step.sample_measurement_ops(
		[vqubit.measure(q0, key="a"), vqubit.measure(q1, key="b")], repetitions=20, seed=0
	)
	assert len(records["a"]) == len(records["b"]) == 20
	assert [row[0] for row in records["a"]] == [row[0] for row in records["b"]]
	assert {tuple(row) for row in records["a"]} <= {(0,), (1,)}


def test_state_vector_step_bloch_vector_and_reduced_density_matrix() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	# |+> on q0 and |0> on q1, in the (q0, q1) order.
	plus = vqubit.StateVectorStepResult(qubits=(q0, q1), measurements={}, state_vector=[SQRT_HALF, 0, SQRT_HALF, 0])
	assert plus.bloch_vector_of(q0) == pytest.approx([1.0, 0.0, 0.0])
	assert plus.bloch_vector_of(q1) == pytest.approx([0.0, 0.0, 1.0])
	assert_matrix_close(plus.density_matrix_of([q0]), [[0.5, 0.5], [0.5, 0.5]])
	assert_matrix_close(plus.density_matrix_of([q1]), [[1, 0], [0, 0]])


def test_state_vector_trial_result_uses_the_mixin_string_forms() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.StateVectorTrialResult(None, {}, qubits=(q0,), state_vector=[1, 0])
	assert result.dirac_notation(1) == "1.0|0⟩"
	assert str(result) == "1.000|0⟩"


def test_clifford_simulator_step_result_exposes_the_clifford_state() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	final = vqubit.CliffordSimulator().simulate(_bell_circuit()).final_state
	step = vqubit.CliffordSimulatorStepResult(qubits=(q0, q1), measurements={}, state=final)
	assert step.state is final
	assert step.clifford_state is final
	assert_vector_close(step.state_vector(), [SQRT_HALF, 0, 0, SQRT_HALF])


def test_simulation_state_base_class_is_abstract_but_mergeable() -> None:
	q0 = vqubit.LineQubit(0)
	state = vqubit.SimulationState((q0,))
	assert state.qubits == (q0,)
	assert state.create_merged_state() is state
	assert repr(state) == "vqubit.SimulationState(1 qubits)"
	with pytest.raises(NotImplementedError):
		state.apply_operation(vqubit.X(q0))
	with pytest.raises(NotImplementedError):
		state.measure([q0], "m")
	with pytest.raises(NotImplementedError):
		state.reset([q0])


def test_state_vector_simulation_state_measures_and_rejects_channels() -> None:
	q0 = vqubit.LineQubit(0)
	state = vqubit.StateVectorSimulationState(qubits=(q0,))
	state.apply_operation(vqubit.X(q0))
	state.apply_operation(vqubit.measure(q0, key="m"))
	assert state.log_of_measurement_results == {"m": 1}
	assert_vector_close(state.state_vector(), [0j, 1 + 0j])
	assert state.measure([q0], "again") == 1
	assert state.log_of_measurement_results == {"m": 1, "again": 1}
	# A channel has no unitary, so the state-vector state cannot apply it.
	with pytest.raises(TypeError, match="has no unitary"):
		state.apply_operation(vqubit.reset(q0))


def test_density_matrix_simulation_state_measures_and_rejects_unknown_gates() -> None:
	q0 = vqubit.LineQubit(0)
	state = vqubit.DensityMatrixSimulationState(qubits=(q0,))
	state.apply_operation(vqubit.X(q0))
	state.apply_operation(vqubit.measure(q0, key="m"))
	assert state.log_of_measurement_results == {"m": 1}
	assert_matrix_close(state.density_matrix(), [[0, 0], [0, 1]])
	assert state.measure([q0], "again") == 1
	assert state.log_of_measurement_results == {"m": 1, "again": 1}
	# A gate with neither a unitary nor Kraus operators cannot be applied either.
	with pytest.raises(TypeError, match="neither a unitary nor Kraus"):
		state.apply_operation(_NoUnitaryGate().on(q0))


def test_clifford_tableau_state_exposes_the_tableau_and_rejects_other_gates() -> None:
	q0 = vqubit.LineQubit(0)
	tableau = native.NativeTableau(1)
	state = vqubit.CliffordTableauSimulationState(tableau=tableau, qubits=(q0,))
	assert state.tableau is tableau
	state.apply_operation(vqubit.H(q0))
	assert_vector_close(state.state_vector(), [SQRT_HALF, SQRT_HALF])
	with pytest.raises(TypeError, match="not a Clifford gate"):
		state.apply_operation(vqubit.T(q0))


def test_stabilizer_state_ch_form_compares_by_amplitudes_and_reprs() -> None:
	one = vqubit.StabilizerStateChForm([1, 0], 1)
	same = vqubit.StabilizerStateChForm([1, 0])
	other = vqubit.StabilizerStateChForm([0, 1], 1)
	assert one == same
	assert one != other
	assert repr(one) == "vqubit.StabilizerStateChForm(1 qubits)"


def test_stabilizer_ch_form_state_evolves_measures_and_rejects_channels() -> None:
	q0 = vqubit.LineQubit(0)
	state = vqubit.StabilizerChFormSimulationState(qubits=(q0,))
	assert_vector_close(state.state_vector(), [1 + 0j, 0j])
	state.apply_operation(vqubit.X(q0))
	assert_vector_close(state.state_vector(), [0j, 1 + 0j])
	state.apply_operation(vqubit.measure(q0, key="m"))
	assert state.log_of_measurement_results == {"m": 1}
	assert_vector_close(state.state_vector(), [0j, 1 + 0j])
	assert state.measure([q0], "again") == 1
	assert state.log_of_measurement_results == {"m": 1, "again": 1}
	with pytest.raises(TypeError, match="has no unitary"):
		state.apply_operation(vqubit.reset(q0))


def test_simulation_product_state_exposes_its_factors() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	first = vqubit.StateVectorSimulationState(qubits=(q0,))
	second = vqubit.StateVectorSimulationState(qubits=(q1,))
	product = vqubit.SimulationProductState([first, second])
	assert product.states == (first, second)
	assert repr(product) == "vqubit.SimulationProductState(2 states)"


def test_stabilizer_sampler_sample_matches_run() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
	result = vqubit.StabilizerSampler(seed=0).sample(circuit, repetitions=3)
	assert result.measurements["m"] == [[1]] * 3


def test_measurement_helpers_can_write_into_an_output_buffer() -> None:
	bell = [SQRT_HALF, 0, 0, SQRT_HALF]
	out: list[complex] = []
	collapsed, outcome = vqubit.measure_state_vector(bell, [0], out=out, seed=0)
	assert out == collapsed
	assert outcome in (0, 1)
	rho = vqubit.density_matrix_from_state_vector(bell)
	rho_out: list[list[complex]] = []
	collapsed_rho, rho_outcome = vqubit.measure_density_matrix(rho, [0], out=rho_out, seed=0)
	assert rho_out == collapsed_rho
	assert rho_outcome in (0, 1)
