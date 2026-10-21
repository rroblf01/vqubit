"""The Cirq-compatible shapes, pinned without needing Cirq installed.

Every one of these is a call that works in Cirq and used to fail here (or to
return something else): the point is that the *shape* stays right, so the tests
run with the wheel and pytest only.
"""

from __future__ import annotations

import math

import pytest

import vqubit


def test_the_final_helpers_take_a_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0))
	vector = vqubit.final_state_vector(circuit)
	assert [round(abs(z), 6) for z in vector] == [0.707107, 0.707107]
	# The same call with a qubit order and an initial state.
	assert len(vqubit.final_state_vector(circuit, qubit_order=[q0, q1], initial_state=0)) == 4
	# And with terminal measurements, which a state vector cannot carry otherwise.
	measured = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	assert len(vqubit.final_state_vector(measured, ignore_terminal_measurements=True)) == 2
	with pytest.raises(ValueError, match="measurements"):
		vqubit.final_state_vector(measured)
	matrix = vqubit.final_density_matrix(vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1)))
	assert len(matrix) == 4 and len(matrix[0]) == 4
	# A simulation result is still accepted, which is what the helpers took.
	result = vqubit.Simulator().simulate(circuit)
	assert vqubit.final_state_vector(result) == vqubit.final_state_vector(circuit)
	density_result = vqubit.DensityMatrixSimulator().simulate(circuit)
	assert vqubit.final_density_matrix(density_result) == vqubit.final_density_matrix(circuit)


def test_the_density_result_exposes_an_array() -> None:
	q0 = vqubit.LineQubit(0)
	result = vqubit.DensityMatrixSimulator().simulate(vqubit.Circuit(vqubit.H(q0)))
	matrix = result.final_density_matrix
	# A nested list of lists, not a wrapper: `matrix[0][0]` has to work.
	assert matrix[0][0] == pytest.approx(0.5)
	assert [list(row) == list(result.density_matrix()[index]) for index, row in enumerate(matrix)] == [True, True]
	assert result.qubits == (q0,)
	assert result.dirac_notation() == "0.707|0⟩ + 0.707|1⟩"


def test_the_measurement_helpers_return_the_outcomes_first() -> None:
	bell = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	outcomes, collapsed = vqubit.measure_state_vector(bell, [0], seed=0)
	assert isinstance(outcomes, list) and len(outcomes) == 1 and outcomes[0] in (0, 1)
	assert len(collapsed) == 4
	rho = vqubit.density_matrix_from_state_vector(bell)
	rho_outcomes, collapsed_rho = vqubit.measure_density_matrix(rho, [0], seed=0)
	assert isinstance(rho_outcomes, list) and len(rho_outcomes) == 1
	assert len(collapsed_rho) == 4


def test_to_valid_state_vector_takes_the_number_of_qubits() -> None:
	assert vqubit.to_valid_state_vector(2, 2) == [0j, 0j, 1 + 0j, 0j]
	assert vqubit.to_valid_state_vector(2, qid_shape=(2, 2)) == [0j, 0j, 1 + 0j, 0j]
	# The old positional shape keeps working.
	assert vqubit.to_valid_state_vector(2, (2, 2)) == [0j, 0j, 1 + 0j, 0j]
	assert len(vqubit.to_valid_density_matrix(0, 1)) == 2


def test_operations_carry_tags() -> None:
	q0 = vqubit.LineQubit(0)
	operation = vqubit.H(q0)
	assert operation.tags == ()
	tagged = operation.with_tags("a", "b")
	assert tagged.tags == ("a", "b")
	assert isinstance(tagged, vqubit.TaggedOperation)
	assert tagged.untagged == operation
	assert tagged.with_tags("c").tags == ("a", "b", "c")
	assert operation.with_tags() is operation


def test_a_circuit_has_tags_and_adds_and_multiplies() -> None:
	q0, _q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0))
	assert circuit.tags == ()
	assert circuit.with_tags("t").tags == ("t",)
	assert vqubit.Circuit(vqubit.X(q0), tags=("t",)).tags == ("t",)
	# + concatenates circuits, appends operations, and * repeats.
	assert len(circuit + vqubit.Circuit(vqubit.X(q0))) == 2
	assert len(circuit + vqubit.X(q0)) == 2
	# Cirq has no `operation + circuit`; here the operation comes first.
	assert len(vqubit.X(q0) + circuit) == 2
	assert len(circuit * 3) == 3
	assert len(3 * circuit) == 3
	assert len(sum([circuit, circuit], vqubit.Circuit())) == 2


def test_the_diagram_takes_a_qubit_order() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0))
	assert circuit.to_text_diagram(qubit_order=[q1, q0]) == "1: ───────\n0: ───H───"
	assert circuit.to_text_diagram() == "0: ───H───"
	# The other Cirq arguments are accepted (the renderer is ASCII and does not
	# use them today).
	assert circuit.to_text_diagram(use_unicode_characters=False, precision=6, transpose=False)


def test_a_circuit_operation_repeats_and_resolves_parameters() -> None:
	q0 = vqubit.LineQubit(0)
	operation = vqubit.CircuitOperation(vqubit.FrozenCircuit(vqubit.X(q0)))
	assert operation.repeat(3).repetitions == 3
	assert operation.repeat(2).circuit == operation.circuit
	parameterised = vqubit.CircuitOperation(
		vqubit.FrozenCircuit(vqubit.XPowGate(exponent=vqubit.Parameter("t"))(q0)),
		param_resolver=vqubit.ParamResolver({"t": 0.5}),
	)
	assert parameterised.circuit == vqubit.Circuit((vqubit.X**0.5)(q0))


def test_measure_takes_an_invert_mask() -> None:
	q0 = vqubit.LineQubit(0)
	gate = vqubit.measure(q0, key="m", invert_mask=(True,)).gate
	assert isinstance(gate, vqubit.MeasurementGate)
	assert gate.invert_mask == (True,)


def test_stratified_circuit_takes_the_categories_by_name() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.H(q0))
	named = vqubit.stratified_circuit(circuit, categories=[vqubit.HPowGate, vqubit.XPowGate])
	positional = vqubit.stratified_circuit(circuit, vqubit.HPowGate, vqubit.XPowGate)
	assert named == positional


def test_cz_target_gateset_accepts_the_partial_cz_flag() -> None:
	assert vqubit.CZTargetGateset().allow_partial_czs is False
	assert vqubit.CZTargetGateset(allow_partial_czs=True).allow_partial_czs is True
	# The other Cirq arguments are accepted and kept.
	gateset = vqubit.CZTargetGateset(additional_gates=(), preserve_moment_structure=False, reorder_operations=True)
	assert gateset.preserve_moment_structure is False
	assert gateset.reorder_operations is True


def test_json_round_trips_every_new_gate_type() -> None:
	q0 = vqubit.LineQubit(0)
	matrix = [[1 + 0j, 0j], [0j, 1j]]
	gates = [
		vqubit.bit_flip(0.1),
		vqubit.phase_flip(0.2),
		vqubit.BitPhaseFlipChannel(0.3),
		vqubit.depolarize(0.1, 2),
		vqubit.amplitude_damp(0.2),
		vqubit.phase_damp(0.3),
		vqubit.generalized_amplitude_damp(0.2, 0.4),
		vqubit.ResetChannel(),
		vqubit.KrausChannel([matrix]),
		vqubit.MixedUnitaryChannel([(0.5, vqubit.unitary(vqubit.I)), (0.5, vqubit.unitary(vqubit.X))]),
		vqubit.WaitGate(0.3),
		vqubit.MatrixGate(vqubit.unitary(vqubit.H)),
		vqubit.FSimGate(theta=0.3, phi=0.6),
	]
	for gate in gates:
		restored = vqubit.read_json(vqubit.to_json(gate))
		assert restored == gate, type(gate).__name__
	# A tagged operation survives too, with the tag.
	operation = vqubit.H(q0).with_tags("a")
	restored = vqubit.read_json(vqubit.to_json(operation))
	assert isinstance(restored, vqubit.TaggedOperation)
	assert restored.tags == ("a",)
	assert isinstance(operation, vqubit.TaggedOperation)
	assert restored.untagged == operation.untagged
	# And the JSON names are Cirq's, not new ones.
	assert '"cirq_type": "BitFlipChannel"' in vqubit.to_json(vqubit.bit_flip(0.1))
	assert '"p": 0.1' in vqubit.to_json(vqubit.bit_flip(0.1))


def test_a_channel_whose_effects_are_not_diagonal() -> None:
	# ``K0`` has an off-diagonal term, so ``K0†K0`` is not diagonal and the
	# trajectory sampler has to look at the coherences of its qubit. The second
	# operator completes the channel: ``K0†K0 + K1†K1 = I``.
	root = math.sqrt(0.5)
	first = [[root, 0.2], [0.0, root]]
	# ``I - K0†K0``, factorised by hand (it is PSD: determinant 0.21).
	second = [[math.sqrt(0.5), -0.1414213562373095], [0.0, math.sqrt(0.46 - 0.02)]]
	channel = vqubit.KrausChannel([first, second])
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), channel(q0), vqubit.measure(q0, key="m"))
	state = vqubit.Simulator(seed=3).run(circuit, repetitions=2000)
	assert all(row[0] in (0, 1) for row in state.measurements["m"])
	# The state after the channel is normalised in every trajectory.
	for seed in range(20):
		vector = vqubit.Simulator(seed=seed).simulate(vqubit.Circuit(vqubit.H(q0), channel(q0))).state_vector()
		assert abs(sum(abs(z) ** 2 for z in vector) - 1) < 1e-9


def test_the_circuit_operators_refuse_what_they_cannot_do() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	with pytest.raises(TypeError):
		_ = circuit + 1
	# `c * 1.5` and `c * True` are refused, which Python turns into a TypeError.
	with pytest.raises(TypeError):
		_ = circuit * 1.5  # type: ignore[unsupported-operator]
	with pytest.raises(TypeError):
		_ = circuit * True  # type: ignore[unsupported-operator]
	# A QubitOrder object, not just a list, decides the diagram's order.
	assert circuit.to_text_diagram(qubit_order=vqubit.QubitOrder.DEFAULT) == "0: ───H───"


def test_a_two_qubit_channel_uses_the_coherences_of_both_qubits() -> None:
	# The target subsystem is not the last qubits, so the bit scattering has to
	# interleave the "rest" and the target indices correctly.
	mixed = vqubit.KrausChannel([vqubit.unitary(vqubit.CNOT)])
	q0, _q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.H(q0), mixed(q0, q2), vqubit.measure(q0, q2, key="m"))
	result = vqubit.Simulator(seed=5).run(circuit, repetitions=500)
	# The Bell pair ends up perfectly correlated, whatever the branch.
	assert all(row[0] == row[1] for row in result.measurements["m"])
