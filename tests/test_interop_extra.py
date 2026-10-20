"""Extra tests for the Quirk and JSON interop helpers.

``tests/test_interop.py`` covers the shape of a Quirk URL and
``tests/test_json.py`` the common JSON values. This file documents the whole
Quirk gate set in both directions (including the errors for gates Quirk has no
spelling for) and the JSON representation of the less common gates,
measurements and malformed input.
"""

from __future__ import annotations

import base64
import json
import math
from typing import Any

import pytest

import vqubit

Q0, Q1, Q2 = vqubit.LineQubit.range(3)


def _quirk_payload(url: str) -> dict[str, Any]:
	payload = url.split("circuit=", 1)[1]
	payload += "=" * (-len(payload) % 4)
	decoded = json.loads(base64.urlsafe_b64decode(payload))
	assert isinstance(decoded, dict)
	return decoded


def _columns(url: str) -> list[list[Any]]:
	return _quirk_payload(url)["cols"]


def _roundtrip(value: Any) -> Any:
	return vqubit.read_json(vqubit.to_json(value))


@pytest.mark.parametrize(
	("operation", "column"),
	[
		(vqubit.H(Q0), ["H"]),
		(vqubit.X(Q0), ["X"]),
		((vqubit.X**2)(Q0), ["X^2"]),
		((vqubit.Y**-1)(Q0), ["Y^-1"]),
		(vqubit.S(Q0), ["Z^0.5"]),
		(vqubit.T(Q0), ["Z^0.25"]),
		((vqubit.Z**-0.5)(Q0), ["Z^-0.5"]),
		(vqubit.Rx(math.pi * 0.3)(Q0), ["X^0.3"]),
		(vqubit.Ry(math.pi / 2)(Q0), ["Y^0.5"]),
		(vqubit.Rz(-math.pi / 4)(Q0), ["Z^-0.25"]),
		(vqubit.I(Q0), [1]),
		(vqubit.CNOT(Q0, Q1), ["•", "X"]),
		(vqubit.CZ(Q0, Q1), ["•", "•"]),
		(vqubit.SWAP(Q0, Q1), ["Swap", "Swap"]),
		(vqubit.ControlledGate(vqubit.Y, num_controls=1).on(Q0, Q1), ["•", "Y"]),
		(vqubit.ControlledGate(vqubit.Z**0.5, num_controls=1).on(Q0, Q1), ["•", "Z^0.5"]),
		(vqubit.ControlledGate(vqubit.X**0.5, num_controls=1).on(Q0, Q1), ["•", "X^0.5"]),
	],
)
def test_quirk_url_writes_each_supported_gate(operation: vqubit.Operation, column: list[Any]) -> None:
	assert _columns(vqubit.quirk_url(vqubit.Circuit(operation))) == [column]


def test_quirk_url_writes_measurements_of_any_width() -> None:
	assert _columns(vqubit.quirk_url(vqubit.Circuit(vqubit.measure(Q0, key="m")))) == [["Measure"]]
	assert _columns(vqubit.quirk_url(vqubit.Circuit(vqubit.measure(Q0, Q1, key="m")))) == [["Measure", "Measure"]]


def test_quirk_url_puts_a_measurement_beside_an_unrelated_gate() -> None:
	circuit = vqubit.Circuit([vqubit.measure(Q0, key="m"), vqubit.X(Q1)])
	assert _columns(vqubit.quirk_url(circuit)) == [["Measure", "X"]]


@pytest.mark.parametrize(
	"operation",
	[
		(vqubit.H**0.5)(Q0),
		(vqubit.ISWAP(Q0, Q1)),
		(vqubit.CCX(Q0, Q1, Q2)),
		vqubit.ControlledGate(vqubit.X, num_controls=1, control_values=[0]).on(Q0, Q1),
		vqubit.ControlledGate(vqubit.MeasurementGate("m", num_qubits=1), num_controls=1).on(Q0, Q1),
	],
)
def test_quirk_url_rejects_gates_it_cannot_spell(operation: vqubit.Operation) -> None:
	with pytest.raises(TypeError, match="cannot be written as a Quirk URL"):
		vqubit.quirk_url(vqubit.Circuit(operation))


def test_a_quirk_url_round_trips_a_program() -> None:
	circuit = vqubit.Circuit(vqubit.H(Q0), vqubit.CNOT(Q0, Q1), vqubit.measure(Q0, Q1, key="m"))
	url = vqubit.quirk_url(circuit)
	restored = vqubit.quirk_url_to_circuit(url)

	alice, bob = vqubit.NamedQubit("q0"), vqubit.NamedQubit("q1")
	assert list(restored.all_qubits()) == [alice, bob]
	assert list(restored.all_operations()) == [
		vqubit.H(alice),
		vqubit.CNOT(alice, bob),
		vqubit.measure(alice, bob, key="m"),
	]
	# Exporting the imported circuit gives the same payload back.
	assert vqubit.quirk_url(restored) == url


def test_quirk_json_accepts_text_or_a_mapping() -> None:
	document = {"cols": [["H", 1], ["•", "X"]], "qubits": ["alice", "bob"]}
	from_text = vqubit.quirk_json_to_circuit(json.dumps(document))
	from_mapping = vqubit.quirk_json_to_circuit({"cols": [["H", 1], ["•", "X"]], "qubits": ["alice", "bob"]})
	assert from_text == from_mapping
	alice, bob = vqubit.NamedQubit("alice"), vqubit.NamedQubit("bob")
	assert list(from_text.all_operations()) == [vqubit.H(alice), vqubit.CNOT(alice, bob)]


def test_quirk_json_names_qubits_q0_q1_by_default() -> None:
	circuit = vqubit.quirk_json_to_circuit('{"cols": [["H", "X"]]}')
	assert list(circuit.all_qubits()) == [vqubit.NamedQubit("q0"), vqubit.NamedQubit("q1")]


@pytest.mark.parametrize(
	("entry", "expected"),
	[
		("X^½", 0.5),
		("X^-¼", -0.25),
		("Y^0.25", 0.25),
		("Z^-1", -1.0),
	],
)
def test_quirk_json_reads_power_shorthand(entry: str, expected: float) -> None:
	circuit = vqubit.quirk_json_to_circuit({"cols": [[entry]]})
	gate = next(circuit.all_operations()).gate
	assert isinstance(gate, (vqubit.XPowGate, vqubit.YPowGate, vqubit.ZPowGate))
	assert gate.exponent == pytest.approx(expected)


@pytest.mark.parametrize(
	("column", "gate"),
	[
		(["•", "X"], vqubit.CNOT),
		(["•", "•"], vqubit.CZ),
		(["•", "Z^1"], vqubit.CZ),
		(["•", "Z^0.5"], vqubit.CZ**0.5),
		(["•", "Y^1"], vqubit.ControlledGate(vqubit.Y, num_controls=1)),
		(["•", "X^0.5"], vqubit.ControlledGate(vqubit.X**0.5, num_controls=1)),
	],
)
def test_quirk_json_reads_controlled_columns(column: list[Any], gate: Any) -> None:
	circuit = vqubit.quirk_json_to_circuit({"cols": [column]})
	operations = list(circuit.all_operations())
	assert len(operations) == 1
	assert operations[0].gate == gate
	assert operations[0].qubits == (vqubit.NamedQubit("q0"), vqubit.NamedQubit("q1"))


def test_quirk_json_reads_swaps_and_wires() -> None:
	swapped = vqubit.quirk_json_to_circuit({"cols": [["Swap", "Swap"]]})
	assert list(swapped.all_operations()) == [vqubit.SWAP(vqubit.NamedQubit("q0"), vqubit.NamedQubit("q1"))]
	assert list(vqubit.quirk_json_to_circuit({"cols": [[1, "-"]]}).all_operations()) == []
	assert list(vqubit.quirk_json_to_circuit({"cols": [["H", 1]]}).all_operations()) == [
		vqubit.H(vqubit.NamedQubit("q0"))
	]


def test_quirk_json_reads_one_measurement_per_column() -> None:
	circuit = vqubit.quirk_json_to_circuit({"cols": [["Measure", "Measure"], ["Measure", 1]]})
	operations = list(circuit.all_operations())
	assert len(operations) == 2
	assert operations[0].gate == vqubit.MeasurementGate("m", num_qubits=2)
	assert operations[1].gate == vqubit.MeasurementGate("m", num_qubits=1)
	assert operations[0].qubits == (vqubit.NamedQubit("q0"), vqubit.NamedQubit("q1"))


def test_quirk_json_reports_unsupported_gates_and_dangling_controls() -> None:
	with pytest.raises(ValueError, match="unsupported Quirk gate 'Foo'"):
		vqubit.quirk_json_to_circuit({"cols": [["Foo"]]})
	with pytest.raises(ValueError, match="control has no target"):
		vqubit.quirk_json_to_circuit({"cols": [["•", 1]]})


def test_quirk_round_trips_a_controlled_y() -> None:
	circuit = vqubit.Circuit(vqubit.ControlledGate(vqubit.Y, num_controls=1).on(Q0, Q1))
	vqubit.quirk_url_to_circuit(vqubit.quirk_url(circuit))


def test_quirk_json_rejects_malformed_documents() -> None:
	with pytest.raises(json.JSONDecodeError):
		vqubit.quirk_json_to_circuit("not json")
	with pytest.raises(KeyError):
		vqubit.quirk_json_to_circuit({"qubits": ["q0"]})


@pytest.mark.parametrize(
	"gate",
	[
		vqubit.PhasedXPowGate(phase_exponent=0.25, exponent=0.5),
		vqubit.PhasedXZGate(x_exponent=0.5, z_exponent=0.25, axis_phase_exponent=0.125),
		vqubit.XXPowGate(exponent=0.5),
		vqubit.YYPowGate(exponent=0.5),
		vqubit.ZZPowGate(exponent=0.5),
		vqubit.ISwapPowGate(exponent=0.3),
		vqubit.SqrtISwapGate(),
		vqubit.CSWAP,
	],
)
def test_json_round_trips_the_less_common_gates(gate: vqubit.Gate) -> None:
	assert _roundtrip(gate) == gate


def test_json_round_trips_a_gate_with_several_controls() -> None:
	gate = vqubit.ControlledGate(vqubit.X, num_controls=2, control_values=[True, False])
	decoded = _roundtrip(gate)
	assert decoded == gate
	assert decoded.control_values == (True, False)


def test_json_round_trips_a_moment_with_several_operations() -> None:
	moment = vqubit.Moment([vqubit.H(Q0), vqubit.CNOT(Q1, Q2)])
	assert _roundtrip(moment) == moment


def test_json_round_trips_measurements_with_an_invert_mask() -> None:
	gate = vqubit.MeasurementGate("m", num_qubits=2, invert_mask=(True, False))
	decoded = _roundtrip(gate)
	assert decoded == gate
	assert decoded.invert_mask == (True, False)
	operation = gate.on(Q0, Q1)
	assert _roundtrip(operation) == operation


def test_json_round_trips_a_circuit_of_grid_and_named_qubits() -> None:
	grid = vqubit.GridQubit(1, 2)
	alice, bob = vqubit.NamedQubit("alice"), vqubit.NamedQubit("bob")
	circuit = vqubit.Circuit(vqubit.H(grid), vqubit.CNOT(alice, bob), vqubit.measure(alice, key="m"))
	assert _roundtrip(circuit) == circuit


def test_json_round_trips_a_gate_operation() -> None:
	operation = vqubit.CZ(Q0, Q1)
	assert _roundtrip(operation) == operation


def test_json_rejects_malformed_text_and_unknown_types() -> None:
	with pytest.raises(json.JSONDecodeError):
		vqubit.read_json("{")
	with pytest.raises(ValueError, match="no decoder registered"):
		vqubit.read_json('{"cirq_type": "Nope"}')
	with pytest.raises(KeyError):
		vqubit.read_json('{"cirq_type": "LineQubit"}')


def test_json_still_rejects_values_it_cannot_encode() -> None:
	with pytest.raises(TypeError, match="cannot serialize"):
		vqubit.to_json(object())


def test_the_cirq_bridge_round_trips_a_circuit() -> None:
	cirq = pytest.importorskip("cirq")
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="m"))

	converted = vqubit.to_cirq(circuit)
	assert isinstance(converted, cirq.Circuit)
	assert vqubit.from_cirq(converted) == circuit
