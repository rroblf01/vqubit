"""Tests for interop import, neutral atoms, extra protocols and testing helpers."""

from __future__ import annotations

import pytest

import vqubit
from vqubit import testing


def test_quirk_roundtrip() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	url = vqubit.quirk_url(circuit)
	recovered = vqubit.quirk_url_to_circuit(url)
	assert testing._without_measurements(recovered).all_qubits() is not None
	assert len(list(recovered.all_operations())) == 2


def test_quirk_json_import() -> None:
	circuit = vqubit.quirk_json_to_circuit('{"cols":[["H",1],["•","X"]],"qubits":["a","b"]}')
	operations = list(circuit.all_operations())
	assert len(operations) == 2
	assert operations[1].gate == vqubit.CNOT


def test_neutral_atoms() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.is_native_neutral_atom_gate(vqubit.CZ)
	assert vqubit.is_native_neutral_atom_gate(vqubit.CZPowGate(exponent=0.5))
	assert not vqubit.is_native_neutral_atom_gate(vqubit.X)
	assert vqubit.is_native_neutral_atom_op(vqubit.CZ(q0, q1))
	assert not vqubit.is_native_neutral_atom_op(vqubit.X(q0))


def test_circuit_diagram_info() -> None:
	q0 = vqubit.LineQubit(0)
	info = vqubit.circuit_diagram_info(vqubit.H(q0))
	assert info.wire_symbols == ("H",)
	args = vqubit.CircuitDiagramInfoArgs()
	assert args.use_unicode_characters is True


def test_mixture_protocol() -> None:
	mixed = vqubit.MixedUnitaryChannel([(0.5, vqubit.unitary(vqubit.I)), (0.5, vqubit.unitary(vqubit.X))])
	assert vqubit.has_mixture(mixed)
	m = vqubit.mixture(mixed)
	assert m is not None
	vqubit.validate_mixture(m)
	with pytest.raises(ValueError):
		vqubit.validate_mixture([(0.4, vqubit.unitary(vqubit.I))])


def test_parameter_protocols() -> None:
	parameter = vqubit.Parameter("t")
	gate = vqubit.rx(parameter)
	assert vqubit.parameter_names(gate) == {"t"}
	assert vqubit.parameter_symbols(gate) == {parameter}
	resolved = vqubit.resolve_parameters_once(gate, vqubit.ParamResolver({"t": 0.5}))
	assert not vqubit.is_parameterized(resolved)


def test_apply_args_and_act_on() -> None:
	q0 = vqubit.LineQubit(0)
	state = vqubit.StateVectorSimulationState(qubits=(q0,))
	vqubit.act_on(vqubit.X(q0), state)
	assert state.state_vector() == [0j, 1 + 0j]
	args = vqubit.ApplyUnitaryArgs([1 + 0j, 0j], axes=[0])
	applied = vqubit.apply_unitaries([vqubit.X], args)
	assert applied == [0j, 1 + 0j]


def test_key_paths_and_json_helpers() -> None:
	assert str(vqubit.with_key_path("m", ["a", "b"])) == "a.b.m"
	assert str(vqubit.with_key_path_prefix("m", ["p"])) == "p.m"
	assert vqubit.measurement_key_name(vqubit.measure(vqubit.LineQubit(0), key="k")) == "k"
	assert vqubit.json_namespace(vqubit.X) == ""
	assert vqubit.obj_to_dict_helper(vqubit.OpIdentifier(vqubit.XPowGate, 1), ["num_qubits"]) == {"num_qubits": 1}


def test_supports_protocols() -> None:
	assert isinstance(vqubit.X, vqubit.SupportsUnitary)
	assert isinstance(vqubit.X, vqubit.SupportsExplicitNumQubits)
	assert isinstance(vqubit.X, vqubit.SupportsParameterization)


def test_trace_distance_from_angle_list() -> None:
	assert vqubit.trace_distance_from_angle_list([0, 0]) == pytest.approx(0)
	assert vqubit.trace_distance_from_angle_list([0, 3.141592653589793]) == pytest.approx(1)


def test_testing_helpers() -> None:
	testing.assert_implements_consistent_protocols(vqubit.X)
	testing.assert_has_consistent_apply_unitary(vqubit.CNOT)
	testing.assert_decompose_is_consistent_with_unitary(vqubit.TOFFOLI)
	testing.assert_pauli_expansion_is_consistent_with_unitary(vqubit.X)
	testing.assert_qasm_is_consistent_with_unitary(vqubit.Circuit(vqubit.H(vqubit.LineQubit(0))))
	testing.assert_has_consistent_qid_shape(vqubit.X)
	tester = testing.equals_tester()
	tester(vqubit.X, vqubit.X)
