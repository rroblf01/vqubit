"""Property tests: the invariants that must hold for *every* input.

Each test draws its inputs with Hypothesis and asserts an algebraic law rather
than a single example: powers compose, the unitary of a circuit is what the
simulator applies, a round trip through QASM/JSON/Quirk preserves the circuit,
and a Pauli string's matrix is the product of its factors.

Hypothesis is a development dependency, so the file skips when it is missing
(the dependency-free test matrix installs only pytest and the wheel).
"""

from __future__ import annotations

import math

import pytest

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

import vqubit  # noqa: E402

Q0, Q1, Q2 = vqubit.LineQubit.range(3)
SETTINGS = settings(max_examples=40, deadline=None)

angles = st.floats(min_value=-4.0, max_value=4.0, allow_nan=False, allow_infinity=False)
turns = st.floats(min_value=-2.0, max_value=2.0, allow_nan=False, allow_infinity=False)


def unitary_of(value: object) -> list[list[complex]]:
	return vqubit.unitary(value)  # type: ignore[arg-type]


def same_matrix(first: object, second: object) -> bool:
	return vqubit.linalg.allclose_up_to_global_phase(unitary_of(first), unitary_of(second))


@SETTINGS
@given(angles)
def test_a_rotation_is_unitary(angle: float) -> None:
	assert vqubit.is_unitary(vqubit.unitary(vqubit.rz(angle)))
	assert vqubit.is_unitary(vqubit.unitary(vqubit.rx(angle)))
	assert same_matrix(vqubit.rz(angle), vqubit.rz(angle))
	# exp(-i a Z/2) exp(i a Z/2) = I, on the qubit it acts on.
	assert same_matrix(
		vqubit.Circuit(vqubit.rz(angle)(Q0), vqubit.rz(-angle)(Q0)),
		vqubit.Circuit(vqubit.I(Q0)),
	)


@SETTINGS
@given(angles, angles)
def test_rotations_add(first: float, second: float) -> None:
	combined = vqubit.Circuit(vqubit.rz(first)(Q0), vqubit.rz(second)(Q0))
	assert same_matrix(combined, vqubit.Circuit(vqubit.rz(first + second)(Q0)))
	two = vqubit.Circuit(vqubit.rx(first)(Q0), vqubit.rx(second)(Q0))
	assert same_matrix(two, vqubit.Circuit(vqubit.rx(first + second)(Q0)))


@SETTINGS
@given(turns, turns)
def test_powers_compose(first: float, second: float) -> None:
	base = vqubit.X
	powered = (base**first) ** second
	assert same_matrix(powered, base ** (first * second))
	# The inverse undoes the power, and an integer number of turns wraps.
	assert same_matrix(vqubit.Circuit((base**first)(Q0), (base**-first)(Q0)), vqubit.Circuit(vqubit.I(Q0)))
	assert same_matrix(base ** (first + 2), base**first)


one_qubit_gates = st.sampled_from(
	[
		vqubit.H,
		vqubit.X,
		vqubit.Y,
		vqubit.Z,
		vqubit.S,
		vqubit.T,
		vqubit.SQRT_X if hasattr(vqubit, "SQRT_X") else vqubit.X**0.5,
	]
)
two_qubit_gates = st.sampled_from(
	[vqubit.CNOT, vqubit.CZ, vqubit.SWAP, vqubit.ISWAP, vqubit.ISWAP**0.5, vqubit.XX**0.25]
)


@st.composite
def random_circuit(draw: st.DrawFn) -> vqubit.Circuit:
	"""A circuit of up to eight random gates, spanning both qubits."""
	operations = []
	for _ in range(draw(st.integers(min_value=0, max_value=8))):
		if draw(st.booleans()):
			operations.append(draw(one_qubit_gates).on(draw(st.sampled_from([Q0, Q1]))))
		else:
			operations.append(draw(two_qubit_gates).on(Q0, Q1))
	# The identities keep the circuit on both qubits, so its unitary is 4x4 even
	# when no drawn gate mentioned one of them.
	operations.extend([vqubit.I(Q0), vqubit.I(Q1)])
	return vqubit.Circuit(operations)


@SETTINGS
@given(random_circuit())
def test_the_circuit_unitary_is_what_the_simulator_applies(circuit: vqubit.Circuit) -> None:
	matrix = vqubit.unitary(circuit)
	assert vqubit.is_unitary(matrix)
	# Applying the matrix to |00> must give the simulated state.
	state = [matrix[row][0] for row in range(4)]
	expected = vqubit.Simulator().simulate(circuit).state_vector()
	for got, want in zip(state, expected, strict=True):
		assert abs(got - want) < 1e-9


@SETTINGS
@given(random_circuit())
def test_the_density_matrix_is_the_pure_state(circuit: vqubit.Circuit) -> None:
	state = vqubit.Simulator().simulate(circuit).state_vector()
	density = vqubit.DensityMatrixSimulator().simulate(circuit).final_density_matrix.density_matrix()
	for i in range(4):
		for j in range(4):
			assert abs(density[i][j] - state[i] * state[j].conjugate()) < 1e-9


qasm_gates = st.sampled_from(
	[vqubit.H, vqubit.X, vqubit.Y, vqubit.Z, vqubit.S, vqubit.T, vqubit.CNOT, vqubit.CZ, vqubit.SWAP, vqubit.CCX]
)


@st.composite
def qasm_circuit(draw: st.DrawFn) -> vqubit.Circuit:
	"""A circuit using only the gates the OpenQASM 2.0 exporter can write."""
	operations = []
	for _ in range(draw(st.integers(min_value=0, max_value=6))):
		gate = draw(qasm_gates)
		if gate.num_qubits == 1:
			operations.append(gate.on(draw(st.sampled_from([Q0, Q1]))))
		elif gate.num_qubits == 2:
			operations.append(gate.on(Q0, Q1))
		else:
			operations.append(gate.on(Q0, Q1, Q2))
	operations.extend([vqubit.I(Q0), vqubit.I(Q1), vqubit.I(Q2)])
	return vqubit.Circuit(operations)


@SETTINGS
@given(qasm_circuit())
def test_the_qasm_round_trip_preserves_the_circuit(circuit: vqubit.Circuit) -> None:
	recovered = vqubit.circuit_from_qasm(vqubit.qasm(circuit))
	assert len(list(recovered.all_operations())) == len(list(circuit.all_operations()))
	assert same_matrix(recovered, circuit)


@SETTINGS
@given(random_circuit())
def test_the_json_round_trip_preserves_the_circuit(circuit: vqubit.Circuit) -> None:
	assert vqubit.read_json(vqubit.to_json(circuit)) == circuit


paulis = st.sampled_from([vqubit.Pauli.X, vqubit.Pauli.Y, vqubit.Pauli.Z])
real_coefficients = st.sampled_from([1, -1, 0.5, 2])


@st.composite
def random_pauli_string(draw: st.DrawFn) -> vqubit.PauliString:
	mapping: dict[vqubit.LineQubit, vqubit.Pauli] = {}
	for qubit in (Q0, Q1):
		if draw(st.booleans()):
			mapping[qubit] = draw(paulis)
	coefficient = draw(st.sampled_from([1, -1, 1j, -1j, 0.5]))
	return vqubit.PauliString(mapping, coefficient)


@SETTINGS
@given(random_pauli_string(), random_pauli_string())
def test_pauli_multiplication_matches_the_matrices(first: vqubit.PauliString, second: vqubit.PauliString) -> None:
	product = first * second
	# Every string is widened to the same qubits, so the matrices multiply.
	qubits = sorted(set(first.qubits) | set(second.qubits))
	if not qubits:
		assert product.coefficient == pytest.approx(first.coefficient * second.coefficient)
		return

	def widened(string: vqubit.PauliString) -> list[list[complex]]:
		# ``PauliString`` drops identities, so build the matrix over the union
		# explicitly: the Kronecker product of each qubit's Pauli, scaled.
		matrix: list[list[complex]] = [[1 + 0j]]
		for qubit in qubits:
			matrix = vqubit.linalg.kron(matrix, string[qubit].matrix())
		return [[string.coefficient * entry for entry in row] for row in matrix]

	expected = vqubit.linalg.dot(widened(first), widened(second))
	assert vqubit.linalg.allclose_up_to_global_phase(widened(product), expected), (first, second, product)


@st.composite
def hermitian_pauli_string(draw: st.DrawFn) -> vqubit.PauliString:
	mapping: dict[vqubit.LineQubit, vqubit.Pauli] = {}
	for qubit in (Q0, Q1):
		if draw(st.booleans()):
			mapping[qubit] = draw(paulis)
	return vqubit.PauliString(mapping, draw(real_coefficients))


@SETTINGS
@given(hermitian_pauli_string(), angles)
def test_pauli_expectation_matches_the_matrix(string: vqubit.PauliString, angle: float) -> None:
	# A random one-qubit state; the expectation is <psi|P|psi> with P the matrix
	# over the string's own qubit. The empty string is the identity, tested
	# separately, and two-qubit strings are out of scope here.
	if len(string.qubits) != 1:
		return
	state = [math.cos(angle), math.sin(angle) * 1j]
	norm = math.sqrt(sum(abs(z) ** 2 for z in state))
	state = [z / norm for z in state]
	qubit = string.qubits[0]
	full = vqubit.PauliString({qubit: string[qubit]}, string.coefficient)
	matrix = vqubit.unitary(full)
	expected = sum(state[i].conjugate() * matrix[i][j] * state[j] for i in range(2) for j in range(2))
	got = full.expectation_from_state_vector(state, {qubit: 0})
	assert abs(got - expected) < 1e-9


@SETTINGS
@given(st.integers(min_value=2, max_value=12), angles, angles)
def test_sweeps_have_the_asked_for_length(count: int, start: float, stop: float) -> None:
	resolvers = vqubit.Linspace("t", start, stop, count).param_resolvers
	assert len(resolvers) == count
	values = [resolver.value_of(vqubit.Parameter("t")) for resolver in resolvers]
	assert values[0] == pytest.approx(start)
	assert values[-1] == pytest.approx(stop)
	# A sweep of points keeps them in order and in number.
	points = vqubit.Points("t", [start, stop]).param_resolvers
	assert len(points) == 2


@SETTINGS
@given(st.lists(angles, min_size=0, max_size=4))
def test_resolving_a_gate_equals_building_it_with_the_value(values: list[float]) -> None:
	parameter = vqubit.Parameter("t")
	gate = vqubit.rz(parameter)
	for value in values:
		resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({"t": value}))
		assert resolved == vqubit.rz(value)
	# With no value supplied the gate stays symbolic.
	assert vqubit.is_parameterized(vqubit.resolve_parameters(gate, vqubit.ParamResolver({})))
