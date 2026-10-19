"""Tests for JSON serialization."""

from __future__ import annotations

from typing import Any

import pytest

import vqubit


def _roundtrip(value: object) -> Any:
	return vqubit.read_json(vqubit.to_json(value))


def test_qubit_roundtrip() -> None:
	for qubit in [vqubit.LineQubit(3), vqubit.GridQubit(1, 2), vqubit.NamedQubit("alice")]:
		assert _roundtrip(qubit) == qubit


def test_gate_roundtrip() -> None:
	for gate in [vqubit.X, vqubit.S, vqubit.Rx(0.7), vqubit.CNOT, vqubit.ISWAP, vqubit.CCX]:
		assert _roundtrip(gate) == gate


def test_controlled_gate_roundtrip() -> None:
	q0 = vqubit.LineQubit(0)
	gate = vqubit.X.controlled_by(q0, control_values=[0])
	assert _roundtrip(gate) == gate


def test_measurement_gate_preserves_invert_mask() -> None:
	gate = vqubit.MeasurementGate("m", num_qubits=2, invert_mask=(True,))
	decoded = _roundtrip(gate)
	assert decoded == gate
	assert decoded.invert_mask == (True,)


def test_operation_moment_circuit_roundtrip() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.Rz(0.3)(q2), vqubit.measure(q0, q1, key="m"))
	assert _roundtrip(circuit) == circuit
	moment = vqubit.Moment([vqubit.X(q0)])
	assert _roundtrip(moment) == moment
	operation = vqubit.CNOT(q0, q1)
	assert _roundtrip(operation) == operation


def test_frozen_circuit_roundtrip() -> None:
	q0 = vqubit.LineQubit(0)
	frozen = vqubit.FrozenCircuit(vqubit.H(q0))
	decoded = _roundtrip(frozen)
	assert isinstance(decoded, vqubit.FrozenCircuit)
	assert decoded == frozen


def test_pauli_string_roundtrip() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	ps = vqubit.PauliString({q0: vqubit.Pauli.X, q1: vqubit.Pauli.Z}, coefficient=1j)
	assert _roundtrip(ps) == ps


def test_parameter_and_expression_roundtrip() -> None:
	t = vqubit.Parameter("t")
	assert _roundtrip(t) == t
	gate = vqubit.Rx((t + 1) * 2)
	assert _roundtrip(gate) == gate


def test_gzip_roundtrip() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	assert vqubit.read_json_gzip(vqubit.to_json_gzip(circuit)) == circuit


def test_indentation() -> None:
	text = vqubit.to_json(vqubit.LineQubit(1), indent=2)
	assert "\n" in text


def test_cirq_type_helpers() -> None:
	assert vqubit.json_cirq_type(vqubit.X) == "XPowGate"
	assert vqubit.json_cirq_type(vqubit.LineQubit) == "LineQubit"
	assert vqubit.cirq_type_from_json("LineQubit") is vqubit.LineQubit
	assert vqubit.cirq_type_from_json("nope") is None


def test_the_plain_gates_are_written_under_their_cirq_names() -> None:
	# Cirq writes ``{"cirq_type": "CXPowGate", "exponent": 1.0, ...}`` for
	# ``cirq.CNOT``, so a file should be readable by either library.
	assert vqubit.json_cirq_type(vqubit.CNOT) == "CXPowGate"
	assert vqubit.json_cirq_type(vqubit.CZ) == "CZPowGate"
	assert vqubit.json_cirq_type(vqubit.SWAP) == "SwapPowGate"
	assert vqubit.json_cirq_type(vqubit.CCX) == "CCXPowGate"
	assert vqubit.json_cirq_type(vqubit.ISWAP) == "ISwapPowGate"
	assert '"exponent": 1.0' in vqubit.to_json(vqubit.CNOT)


def test_the_historical_gate_names_still_load() -> None:
	# Before the plain gates became power gates, ``CNOT`` was written as
	# ``CNOTGate``; those files still have to read back.
	assert vqubit.read_json('{"cirq_type": "CNOTGate"}') == vqubit.CNOT
	assert vqubit.read_json('{"cirq_type": "CZGate"}') == vqubit.CZ
	assert vqubit.read_json('{"cirq_type": "SWAPGate"}') == vqubit.SWAP
	assert vqubit.read_json('{"cirq_type": "CCXGate"}') == vqubit.CCX
	assert vqubit.read_json('{"cirq_type": "ISWAPGate"}') == vqubit.ISWAP
	assert vqubit.read_json('{"cirq_type": "CSWAPGate"}') == vqubit.CSWAP


def test_unknown_type_raises() -> None:
	with pytest.raises(ValueError):
		vqubit.read_json('{"cirq_type": "Nope"}')


def test_unsupported_value_raises() -> None:
	with pytest.raises(TypeError):
		vqubit.to_json(object())
