"""Extra protocol tests: unitary/decomposition, keys, commutators and Kraus.

These complement ``test_protocols.py`` and ``test_advanced_protocols.py`` by
covering the protocol paths that accept *every* kind of value (gates,
operations, moments, circuits, Pauli strings and plain objects) and the
fallback implementations used when numpy is absent.
"""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close, assert_vector_close
from vqubit import _compat

SQRT_HALF = math.sqrt(0.5)


def test_unitary_of_pauli_string_is_its_matrix() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	string = vqubit.PauliString({q0: vqubit.Pauli.X, q1: vqubit.Pauli.Z})
	# Qubit order is the sorted order of the string's qubits: X on q0, Z on q1.
	assert_matrix_close(
		vqubit.unitary(string),
		[
			[0j, 0j, 1 + 0j, 0j],
			[0j, 0j, 0j, -1 + 0j],
			[1 + 0j, 0j, 0j, 0j],
			[0j, -1 + 0j, 0j, 0j],
		],
	)


def test_unitary_of_empty_pauli_string_is_the_coefficient() -> None:
	# A Pauli string with no qubits is a 1x1 matrix, and the coefficient counts.
	assert_matrix_close(vqubit.unitary(vqubit.PauliString()), [[1 + 0j]])
	assert_matrix_close(vqubit.unitary(vqubit.PauliString(coefficient=-2j)), [[-2j]])


def test_unitary_of_circuit_operation_repeats_the_circuit() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	# H is its own inverse, so two repetitions are the identity.
	assert_matrix_close(vqubit.unitary(vqubit.CircuitOperation(circuit, repetitions=2)), vqubit.unitary(vqubit.I))
	assert_matrix_close(vqubit.unitary(vqubit.CircuitOperation(circuit, repetitions=3)), vqubit.unitary(vqubit.H))


def test_unitary_falls_back_to_a_gates_decomposition() -> None:
	# UniformSuperpositionGate has no _unitary_, only _decompose_; the protocol
	# must build the matrix from the decomposition.
	gate = vqubit.UniformSuperpositionGate(3, 2)
	matrix = vqubit.unitary(gate)
	assert len(matrix) == 4
	# The first column is the uniform superposition of |00>, |01>, |10>.
	for i in range(3):
		assert matrix[i][0] == pytest.approx(1 / math.sqrt(3))
	assert matrix[3][0] == pytest.approx(0)


def test_unitary_of_a_channel_is_not_available() -> None:
	with pytest.raises(TypeError):
		vqubit.unitary(vqubit.bit_flip(0.1))
	assert not vqubit.has_unitary(vqubit.bit_flip(0.1))
	assert not vqubit.has_unitary(vqubit.measure(vqubit.LineQubit(0)))


def test_unitary_of_a_circuit_with_measurements_raises() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError):
		vqubit.unitary(vqubit.Circuit(vqubit.measure(q0)))
	with pytest.raises(TypeError):
		vqubit.unitary(vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0)))


def test_unitary_of_a_circuit_with_a_channel_raises() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError):
		vqubit.unitary(vqubit.Circuit(vqubit.bit_flip(0.1)(q0)))


def test_unitary_without_numpy_matches_the_numpy_result(monkeypatch: pytest.MonkeyPatch) -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	expected = vqubit.unitary(circuit)
	monkeypatch.setattr(_compat, "HAS_NUMPY", False)
	fallback = vqubit.unitary(circuit)
	assert_matrix_close(fallback, expected)


def test_unitary_fallback_embeds_a_gate_on_the_right_qubits(monkeypatch: pytest.MonkeyPatch) -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.CNOT(q1, q2))
	expected = vqubit.unitary(circuit)
	monkeypatch.setattr(_compat, "HAS_NUMPY", False)
	fallback = vqubit.unitary(circuit)
	assert_matrix_close(fallback, expected)
	# X on q0 flips the high bit; the CNOT then flips q2 when q1 is set. So the
	# input |001> (= index 1) leaves q1 clear and lands on |101> (= index 5).
	assert fallback[5][1] == pytest.approx(1)
	assert fallback[1][5] == pytest.approx(1)


def test_a_gate_whose_decompose_has_the_wrong_signature_has_no_unitary() -> None:
	class BadDecomposeGate(vqubit.Gate):
		def _num_qubits_(self) -> int:
			return 1

		def _decompose_(self) -> list[vqubit.GateOperation]:  # ty: ignore[invalid-method-override]
			return []

	assert not vqubit.has_unitary(BadDecomposeGate())
	with pytest.raises(TypeError):
		vqubit.unitary(BadDecomposeGate())


def test_num_qubits_and_qid_shape_of_operations_and_empty_circuits() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.num_qubits(vqubit.Circuit()) == 0
	assert vqubit.qid_shape(vqubit.Circuit()) == ()
	assert vqubit.num_qubits(vqubit.measure(q0)) == 1
	assert vqubit.num_qubits(vqubit.bit_flip(0.1)) == 1
	assert vqubit.qid_shape(vqubit.bit_flip(0.1)) == (2,)
	with pytest.raises(TypeError):
		vqubit.num_qubits(42)
	with pytest.raises(TypeError):
		vqubit.qid_shape(vqubit.PauliString({q0: vqubit.Pauli.X}))


def test_inverse_of_gates_and_operations() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.inverse(vqubit.H) == vqubit.H
	assert vqubit.inverse(vqubit.S) == vqubit.S**-1
	inverted = vqubit.inverse(vqubit.T(q0))
	assert inverted.qubits == (q0,)
	assert_matrix_close(vqubit.unitary(inverted), vqubit.unitary(vqubit.T**-1))
	# The Toffoli is its own inverse.
	assert_matrix_close(vqubit.unitary(vqubit.inverse(vqubit.TOFFOLI)), vqubit.unitary(vqubit.TOFFOLI))


def test_inverse_is_undefined_for_measurements_and_channels() -> None:
	q0 = vqubit.LineQubit(0)
	for value in (vqubit.measure(q0), vqubit.MeasurementGate("m"), vqubit.bit_flip(0.1), 42):
		with pytest.raises(TypeError):
			vqubit.inverse(value)
	with pytest.raises(TypeError):
		vqubit.inverse(vqubit.Circuit(vqubit.measure(q0)))


def test_pow_of_operations_and_unsupported_values() -> None:
	q0 = vqubit.LineQubit(0)
	half = vqubit.pow(vqubit.X(q0), 0.5)
	assert half.qubits == (q0,)
	assert_matrix_close(vqubit.unitary(half), vqubit.unitary(vqubit.X**0.5))
	assert_matrix_close(vqubit.unitary(vqubit.pow(vqubit.X, 0)), vqubit.unitary(vqubit.I))
	with pytest.raises(TypeError):
		vqubit.pow(vqubit.Circuit(), 2)
	with pytest.raises(TypeError):
		vqubit.pow(42, 2)
	# A gate that is not a power gate cannot be raised to a fractional power.
	with pytest.raises(NotImplementedError):
		vqubit.pow(vqubit.QubitPermutationGate([1, 0]), 0.5)


def test_decompose_of_a_gate_returns_the_gate_itself() -> None:
	assert vqubit.decompose(vqubit.TOFFOLI) == [vqubit.TOFFOLI]
	assert vqubit.decompose_once(vqubit.TOFFOLI) == [vqubit.TOFFOLI]


def test_decompose_expands_operations_recursively() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	operations = vqubit.decompose(vqubit.CSWAP(q0, q1, q2))
	assert len(operations) > 3
	# Everything left is primitive, and no decomposition escapes the qubits.
	assert {type(op.gate) for op in operations} <= {vqubit.HPowGate, vqubit.CXPowGate, vqubit.ZPowGate}
	for op in operations:
		assert set(op.qubits) <= {q0, q1, q2}


def test_decompose_once_stops_after_a_single_level() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	once = vqubit.decompose_once(vqubit.CSWAP(q0, q1, q2))
	assert len(once) == 3
	assert [type(op.gate) for op in once] == [vqubit.CXPowGate, vqubit.CCXPowGate, vqubit.CXPowGate]
	# A primitive operation is returned unchanged.
	assert vqubit.decompose_once(vqubit.H(q0)) == [vqubit.H(q0)]


def test_decompose_of_circuits_and_invalid_input() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	assert vqubit.decompose(circuit) == [vqubit.H(q0)]
	assert vqubit.decompose_once(circuit) == [vqubit.H(q0)]
	with pytest.raises(TypeError):
		vqubit.decompose(42)
	with pytest.raises(TypeError):
		vqubit.decompose_once(42)


def test_is_measurement_accepts_every_value_kind() -> None:
	q0 = vqubit.LineQubit(0)
	measurement = vqubit.measure(q0, key="m")
	assert vqubit.is_measurement(measurement)
	assert vqubit.is_measurement(measurement.gate)
	assert vqubit.is_measurement(vqubit.Moment([measurement]))
	assert vqubit.is_measurement(vqubit.Circuit(measurement))
	assert not vqubit.is_measurement(vqubit.H(q0))
	assert not vqubit.is_measurement(vqubit.Moment())
	assert not vqubit.is_measurement(42)


def test_measurement_key_names_of_circuits_moments_and_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(
		vqubit.measure(q0, key="a"),
		vqubit.X(q1),
		vqubit.measure(q1, key="a"),
	)
	assert vqubit.measurement_key_names(circuit) == {"a"}
	assert vqubit.measurement_keys_touched(circuit) == {"a"}
	assert vqubit.measurement_key_names(vqubit.Moment([vqubit.measure(q0, key="m")])) == {"m"}
	assert vqubit.measurement_key_names(vqubit.measure(q0, key="n")) == {"n"}
	assert vqubit.measurement_key_names(vqubit.MeasurementGate("g")) == {"g"}
	# Non-measurements and plain objects name no keys at all.
	assert vqubit.measurement_key_names(vqubit.H) == set()
	assert vqubit.measurement_key_names(42) == set()
	assert vqubit.measurement_key_names(vqubit.Circuit()) == set()


def test_measurement_key_objects_of_values() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	operation = vqubit.measure(q0, key="a")
	key = vqubit.measurement_key_obj(operation)
	assert key is not None
	assert key == vqubit.MeasurementKey("a")
	assert vqubit.measurement_key_obj(operation.gate) == vqubit.MeasurementKey("a")
	assert vqubit.measurement_key_obj(vqubit.H) is None
	assert vqubit.measurement_key_obj(42) is None

	circuit = vqubit.Circuit(operation, vqubit.measure(q1, key="a"), vqubit.X(q0))
	assert vqubit.measurement_key_objs(circuit) == frozenset({vqubit.MeasurementKey("a")})
	assert vqubit.measurement_key_objs(vqubit.Moment([operation])) == frozenset({vqubit.MeasurementKey("a")})
	assert vqubit.measurement_key_objs(vqubit.H) == frozenset()
	assert vqubit.measurement_key_objs(42) == frozenset()
	assert vqubit.control_keys(operation) == frozenset()
	assert vqubit.control_keys(42) == frozenset()


def test_with_measurement_key_mapping_renames_keys_on_all_value_kinds() -> None:
	q0 = vqubit.LineQubit(0)
	operation = vqubit.measure(q0, key="old")
	renamed = vqubit.with_measurement_key_mapping(
		operation, {vqubit.MeasurementKey("old"): vqubit.MeasurementKey("new")}
	)
	assert isinstance(renamed, vqubit.GateOperation)
	assert renamed.gate == vqubit.MeasurementGate("new")
	assert renamed.qubits == (q0,)

	# The invert mask survives a rename, and an unmapped key is left alone.
	gated = vqubit.with_measurement_key_mapping(vqubit.MeasurementGate("k", invert_mask=[True]), {"k": "k2"})
	assert isinstance(gated, vqubit.MeasurementGate)
	assert gated.key == "k2"
	assert gated.invert_mask == (True,)
	assert vqubit.with_measurement_key_mapping(gated, {"other": "x"}) == gated

	circuit = vqubit.Circuit(operation)
	keys = vqubit.measurement_key_names(vqubit.with_measurement_key_mapping(circuit, {"old": "new"}))
	assert keys == {"new"}

	# Values that do not carry a measurement key pass through unchanged.
	assert vqubit.with_measurement_key_mapping(vqubit.H, {"old": "new"}) == vqubit.H
	assert vqubit.with_measurement_key_mapping(42, {"old": "new"}) == 42


def test_commutes_for_gates_operations_and_disjoint_qubits() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.commutes(vqubit.X, vqubit.X)
	assert not vqubit.commutes(vqubit.X, vqubit.Z)
	assert vqubit.commutes(vqubit.X, vqubit.X**0.5)
	# Operations on disjoint qubits always commute...
	assert vqubit.commutes(vqubit.X(q0), vqubit.Z(q1))
	# ...while X and Z on the same qubit do not.
	assert not vqubit.commutes(vqubit.X(q0), vqubit.Z(q0))
	assert vqubit.commutes(vqubit.X(q0), vqubit.pow(vqubit.X(q0), 0.5))
	# Undecidable cases are reported as not commuting.
	assert not vqubit.commutes(vqubit.X, vqubit.CNOT)
	assert not vqubit.commutes(42, 42)
	assert vqubit.definitely_commutes(vqubit.X, vqubit.X)
	assert not vqubit.definitely_commutes(vqubit.X, vqubit.Z)


def test_is_parameterized_and_resolve_parameters_across_value_kinds() -> None:
	parameter = vqubit.Parameter("t")
	q0 = vqubit.LineQubit(0)
	resolver = vqubit.ParamResolver({parameter: 0.25})

	assert vqubit.is_parameterized(parameter)
	assert not vqubit.is_parameterized(3)
	assert vqubit.is_parameterized(vqubit.rx(parameter))
	assert vqubit.is_parameterized(vqubit.rx(parameter)(q0))
	assert vqubit.is_parameterized(vqubit.Moment([vqubit.rx(parameter)(q0)]))
	assert vqubit.is_parameterized(vqubit.Circuit(vqubit.rx(parameter)(q0)))

	assert vqubit.resolve_parameters(parameter, resolver) == 0.25
	assert vqubit.resolve_parameters(3, resolver) == 3
	resolved_gate = vqubit.resolve_parameters(vqubit.rx(parameter), resolver)
	assert isinstance(resolved_gate, vqubit.RxGate)
	assert resolved_gate.radians == pytest.approx(0.25)
	resolved_op = vqubit.resolve_parameters(vqubit.rx(parameter)(q0), resolver)
	assert resolved_op.qubits == (q0,)
	assert not vqubit.is_parameterized(resolved_op)
	resolved_moment = vqubit.resolve_parameters(vqubit.Moment([vqubit.rx(parameter)(q0)]), resolver)
	assert not vqubit.is_parameterized(resolved_moment)
	resolved_circuit = vqubit.resolve_parameters(vqubit.Circuit(vqubit.rx(parameter)(q0)), resolver)
	assert not vqubit.is_parameterized(resolved_circuit)
	assert len(resolved_circuit) == 1


def test_pauli_expansion_uses_operations_qubits_and_honours_atol() -> None:
	named = vqubit.NamedQubit("alpha")
	expansion = vqubit.pauli_expansion(vqubit.Z(named))
	assert len(expansion) == 1
	assert expansion[vqubit.PauliString({named: vqubit.Pauli.Z})] == pytest.approx(1 + 0j)

	# X on two parallel qubits is the single Pauli string X⊗X.
	q0, q1 = vqubit.LineQubit.range(2)
	parallel = vqubit.pauli_expansion(vqubit.ParallelGate(vqubit.X))
	assert len(parallel) == 1
	assert parallel[vqubit.PauliString({q0: vqubit.Pauli.X, q1: vqubit.Pauli.X})] == pytest.approx(1 + 0j)

	# With a tolerance above H's 1/sqrt(2) coefficients, nothing is kept.
	assert len(vqubit.pauli_expansion(vqubit.H(q0), atol=0.9)) == 0


def test_phase_by_phases_one_qubit_of_a_two_qubit_gate() -> None:
	# Conjugating X⊗X by Z**0.5 on the second qubit turns the second X into Y.
	phased = vqubit.phase_by(vqubit.ParallelGate(vqubit.X), 0.5, 1)
	assert isinstance(phased, vqubit.MatrixGate)
	assert vqubit.num_qubits(phased) == 2
	assert_matrix_close(vqubit.unitary(phased), vqubit.kron(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.Y)))

	# With no phase the gate is unchanged.
	assert_matrix_close(vqubit.unitary(vqubit.phase_by(vqubit.X, 0.5, 0)), vqubit.unitary(vqubit.Y))
	with pytest.raises(ValueError):
		vqubit.phase_by(vqubit.CNOT, 0.5, -1)
	with pytest.raises(ValueError):
		vqubit.phase_by(vqubit.CNOT, 0.5, 2)


def test_phase_by_requires_a_numeric_phase_exponent() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError):
		vqubit.phase_by(vqubit.X(q0), vqubit.Parameter("t"), 0)  # ty: ignore[invalid-argument-type]


def test_apply_unitary_to_circuits_and_two_qubit_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.Z(q0))
	# X|0> = |1>, then Z|1> = -|1>.
	assert_vector_close(vqubit.apply_unitary(circuit, [1 + 0j, 0j]), [0j, -1 + 0j])
	# CNOT maps |10> (= index 2) to |11> (= index 3).
	assert_vector_close(
		vqubit.apply_unitary(vqubit.CNOT(q0, q1), [0j, 0j, 1 + 0j, 0j]),
		[0j, 0j, 0j, 1 + 0j],
	)
	assert vqubit.apply_unitary(vqubit.X(q0), [1 + 0j, 0j]) == [0j, 1 + 0j]


def test_kraus_of_gates_operations_and_invalid_values() -> None:
	q0 = vqubit.LineQubit(0)
	unitary_kraus = vqubit.kraus(vqubit.X(q0))
	assert len(unitary_kraus) == 1
	assert_matrix_close(unitary_kraus[0], vqubit.unitary(vqubit.X))
	channel_kraus = vqubit.kraus(vqubit.bit_flip(0.1)(q0))
	assert len(channel_kraus) == 2
	assert vqubit.has_kraus(vqubit.bit_flip(0.1))
	assert vqubit.has_kraus(vqubit.X)
	assert not vqubit.has_kraus(vqubit.measure(q0))
	assert not vqubit.has_kraus(42)
	with pytest.raises(TypeError):
		vqubit.kraus(vqubit.measure(q0))
	with pytest.raises(TypeError):
		vqubit.kraus(42)


def test_approx_eq_compares_matrices_with_a_tolerance() -> None:
	near_x = vqubit.MatrixGate(
		[
			[1e-7 + 0j, 1 + 0j],
			[1 + 0j, 1e-7 + 0j],
		]
	)
	assert vqubit.approx_eq(near_x, vqubit.X, atol=1e-5)
	assert not vqubit.approx_eq(near_x, vqubit.X, atol=1e-9)
	# Non-unitary values fall back to value equality.
	assert vqubit.approx_eq(vqubit.bit_flip(0.1), vqubit.bit_flip(0.1))
	assert not vqubit.approx_eq(vqubit.bit_flip(0.1), vqubit.bit_flip(0.2))


def test_equal_up_to_global_phase_handles_phases_and_dimensions() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.equal_up_to_global_phase(vqubit.S, vqubit.Z**0.5)
	assert not vqubit.equal_up_to_global_phase(vqubit.S, vqubit.Z**0.25)
	minus_i = vqubit.MatrixGate([[-1j + 0j, 0j], [0j, -1j + 0j]])
	assert vqubit.equal_up_to_global_phase(minus_i, vqubit.I)
	# Different dimensions can never be equal, with or without a phase.
	assert not vqubit.equal_up_to_global_phase(vqubit.X(q0), vqubit.CNOT(q0, q1))


def test_trace_distance_bound_uses_the_hook_and_the_unitary() -> None:
	assert vqubit.trace_distance_bound(vqubit.I) == 0.0
	assert vqubit.trace_distance_bound(vqubit.X) == 1.0
	# Rx(theta) has trace cos(theta/2), so the bound is |sin(theta/2)|.
	assert vqubit.trace_distance_bound(vqubit.Rx(math.pi / 2)) == pytest.approx(SQRT_HALF)
	# A value with no unitary is maximally far from the identity.
	assert vqubit.trace_distance_bound(vqubit.bit_flip(0.1)) == 1.0

	class FixedBound:
		def _trace_distance_bound_(self) -> float:
			return 0.25

	assert vqubit.trace_distance_bound(FixedBound()) == 0.25
