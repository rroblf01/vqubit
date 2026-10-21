"""Tests for the interop helpers (QasmArgs and quirk_url)."""

from __future__ import annotations

import base64
import json
import math

import pytest

import vqubit


def test_qasm_args_defaults() -> None:
	args = vqubit.QasmArgs()
	assert args.precision == 15
	assert args.version == "2.0"
	assert args.qubit_id_map == {}
	assert args.args == {}
	assert args.format("x q[0];") == "x q[0];"


def test_qasm_args_custom() -> None:
	q = vqubit.LineQubit(0)
	args = vqubit.QasmArgs(qubit_id_map={q: "a"}, precision=8, args={"theta": 1})
	assert args.qubit_id_map[q] == "a"
	assert args.precision == 8
	assert args.args == {"theta": 1}


def _decode(url: str) -> dict:
	payload = url.split("circuit=", 1)[1]
	payload += "=" * (-len(payload) % 4)
	return json.loads(base64.urlsafe_b64decode(payload))


def test_quirk_url_encodes_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	url = vqubit.quirk_url(circuit)
	assert url.startswith("https://algassert.com/quirk#circuit=")
	assert _decode(url) == {"cols": [["H", 1], ["•", "X"]], "qubits": ["q0", "q1"]}


def test_quirk_url_handles_controls_spanning_qubits() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit([vqubit.H(q1), vqubit.CNOT(q0, q2)])
	assert _decode(vqubit.quirk_url(circuit))["cols"] == [["•", "H", "X"]]


def test_quirk_url_writes_rotations_as_powers() -> None:
	# Quirk has no ``rx``; a rotation about X is ``X^t``, with ``t = rads/pi``.
	q0 = vqubit.LineQubit(0)
	url = vqubit.quirk_url(vqubit.Circuit(vqubit.Rx(math.pi / 2)(q0)))
	assert _decode(url)["cols"] == [["X^0.5"]]


def test_quirk_url_rejects_unknown_gate() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	with pytest.raises(TypeError):
		vqubit.quirk_url(vqubit.Circuit(vqubit.ISWAP(q0, q1)))
	with pytest.raises(TypeError):
		vqubit.quirk_url(vqubit.Circuit(vqubit.CCX(q0, q1, vqubit.LineQubit(2))))


def test_the_cirq_bridge_round_trips() -> None:
	"""A circuit goes to Cirq and back unchanged, through the shared JSON."""
	cirq = pytest.importorskip("cirq")
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="m"))
	converted = vqubit.to_cirq(circuit)
	assert isinstance(converted, cirq.Circuit)
	assert vqubit.from_cirq(converted) == circuit


def test_the_cirq_bridge_keeps_the_unitary() -> None:
	numpy = pytest.importorskip("numpy", exc_type=ImportError)
	cirq = pytest.importorskip("cirq")
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CZ(q0, q1), vqubit.rz(0.3)(q1))
	converted = vqubit.to_cirq(circuit)
	assert numpy.allclose(numpy.asarray(cirq.unitary(converted)), numpy.asarray(vqubit.unitary(circuit)), atol=1e-9)


def test_channels_and_the_other_gate_types_travel_between_the_two_libraries() -> None:
	"""The JSON of a channel, a wait, a matrix gate or a tagged op, both ways."""
	numpy = pytest.importorskip("numpy", exc_type=ImportError)
	cirq = pytest.importorskip("cirq", exc_type=ImportError)
	import vqubit as _vqubit

	q = _vqubit.LineQubit(0)
	cq = cirq.LineQubit(0)
	mine = [
		_vqubit.bit_flip(0.1),
		_vqubit.phase_flip(0.1),
		_vqubit.depolarize(0.1, 2),
		_vqubit.amplitude_damp(0.2),
		_vqubit.phase_damp(0.2),
		_vqubit.generalized_amplitude_damp(0.2, 0.3),
		_vqubit.ResetChannel(),
		_vqubit.KrausChannel([numpy.eye(2)]),
		_vqubit.WaitGate(0.3),
		_vqubit.MatrixGate(_vqubit.unitary(_vqubit.H)),
		_vqubit.FSimGate(theta=0.3, phi=0.6),
		_vqubit.H(q).with_tags("a"),
	]
	for gate in mine:
		text = _vqubit.to_json(gate)
		# vqubit reads its own JSON back...
		assert _vqubit.read_json(text) is not None
		# ...and Cirq reads it too, which is what the `cirq_type` format is for.
		assert cirq.read_json(json_text=text) is not None
	# And the other direction, from what Cirq writes.
	theirs = [
		cirq.bit_flip(0.1),
		cirq.depolarize(0.1, 2),
		cirq.amplitude_damp(0.2),
		cirq.generalized_amplitude_damp(0.2, 0.3),
		cirq.ResetChannel(),
		cirq.WaitGate(cirq.Duration(nanos=300)),
		cirq.MatrixGate(cirq.unitary(cirq.H)),
		cirq.KrausChannel([numpy.array([[1, 0], [0, 0.9]]), numpy.array([[0, 0.1], [0, 0]])]),
		cirq.FSimGate(theta=0.3, phi=0.6),
		cirq.H(cq).with_tags("a"),
	]
	for gate in theirs:
		restored = _vqubit.read_json(cirq.to_json(gate))
		assert restored is not None
