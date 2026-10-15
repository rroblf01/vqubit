"""Tests for the interop helpers (QasmArgs and quirk_url)."""

from __future__ import annotations

import base64
import json

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


def test_quirk_url_rejects_unknown_gate() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError):
		vqubit.quirk_url(vqubit.Circuit(vqubit.Rx(0.3)(q0)))
