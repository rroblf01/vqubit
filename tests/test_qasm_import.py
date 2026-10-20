"""Tests for the OpenQASM 2.0 importer.

The importer has to read what the exporter writes (so the two round-trip) plus a
useful part of what other tools emit, and it has to say *where* a program is
wrong when it is.
"""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def test_bell_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.circuit_from_qasm(
		HEADER
		+ """
qreg q[2];
creg c[2];
h q[0];
cx q[0],q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
"""
	)
	assert circuit == vqubit.Circuit(
		vqubit.H(q0),
		vqubit.CNOT(q0, q1),
		vqubit.measure(q0, q1, key="c"),
	)


def test_gates_of_qelib1() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.circuit_from_qasm(
		HEADER
		+ """
qreg q[3];
id q[0];
x q[0];
y q[0];
z q[0];
h q[0];
s q[0];
sdg q[0];
t q[0];
tdg q[0];
sx q[0];
sxdg q[0];
cx q[0],q[1];
cz q[0],q[1];
cy q[0],q[1];
ch q[0],q[1];
swap q[0],q[1];
ccx q[0],q[1],q[2];
"""
	)
	assert len(list(circuit.all_operations())) == 17
	assert set(circuit.all_qubits()) == {q0, q1, q2}


def test_parameter_expressions() -> None:
	circuit = vqubit.circuit_from_qasm(
		HEADER
		+ """
qreg q[1];
rx(-pi/2) q[0];
ry(2*pi/3) q[0];
rz(0.5) q[0];
u1(pi/4) q[0];
"""
	)
	rotations = [getattr(operation.gate, "radians", 0.0) for operation in circuit.all_operations()]
	assert [round(angle, 12) for angle in rotations] == [
		round(-math.pi / 2, 12),
		round(2 * math.pi / 3, 12),
		0.5,
		round(math.pi / 4, 12),
	]


def test_u3_matches_the_exact_matrix_up_to_phase() -> None:
	# ``u3`` is written as rotations, which equal the QASM matrix up to a phase.
	circuit = vqubit.circuit_from_qasm(HEADER + "qreg q[1];\nu3(0.3, 0.4, 0.5) q[0];")
	got = vqubit.unitary(circuit)
	cos, sin = math.cos(0.15), math.sin(0.15)
	want = [
		[complex(cos), -complex(math.cos(0.5), math.sin(0.5)) * sin],
		[complex(math.cos(0.4), math.sin(0.4)) * sin, complex(math.cos(0.9), math.sin(0.9)) * cos],
	]
	phase = got[0][0] / want[0][0]
	assert_matrix_close([[entry / phase for entry in row] for row in got], want)


def test_custom_gate_is_expanded() -> None:
	circuit = vqubit.circuit_from_qasm(
		HEADER
		+ """
gate pair(a) q0, q1 { h q0; rz(a) q1; cx q0, q1; }
qreg q[2];
pair(pi/2) q[0], q[1];
"""
	)
	q0, q1 = vqubit.LineQubit.range(2)
	assert list(circuit.all_operations()) == [
		vqubit.H(q0),
		vqubit.rz(math.pi / 2)(q1),
		vqubit.CNOT(q0, q1),
	]


def test_register_arguments_expand_elementwise() -> None:
	circuit = vqubit.circuit_from_qasm(HEADER + "qreg a[2];\nqreg b[2];\nh a;\ncx a,b;\n")
	# ``h a;`` is one gate per wire and ``cx a,b;`` one per pair of wires.
	assert len(list(circuit.all_operations())) == 4
	assert set(circuit.all_qubits()) == set(vqubit.LineQubit.range(4))


def test_reset_and_barrier() -> None:
	circuit = vqubit.circuit_from_qasm(HEADER + "qreg q[1];\nx q[0];\nbarrier q[0];\nreset q[0];\n")
	operations = list(circuit.all_operations())
	assert len(operations) == 2
	assert operations[1].gate == vqubit.ResetChannel()


def test_classical_control() -> None:
	circuit = vqubit.circuit_from_qasm(
		HEADER
		+ """
qreg q[1];
creg c[2];
x q[0];
measure q[0] -> c[1];
if(c==1) x q[0];
"""
	)
	operations = list(circuit.all_operations())
	assert isinstance(operations[-1], vqubit.ClassicallyControlledOperation)
	# ``c == 1`` is "bit 0 is 1 and bit 1 is 0", one condition per bit.
	conditions = operations[-1].conditions
	assert vqubit.KeyCondition("c", 0) in conditions
	assert ~vqubit.KeyCondition("c", 1) in conditions


def test_comments_and_includes_are_ignored() -> None:
	circuit = vqubit.circuit_from_qasm(
		"""// a comment
OPENQASM 2.0;
include "qelib1.inc";
/* a block
   comment */
qreg q[1];
h q[0];
"""
	)
	assert circuit == vqubit.Circuit(vqubit.H(vqubit.LineQubit(0)))


def test_the_importer_and_the_exporter_round_trip() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	corpus = [
		vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="m")),
		vqubit.Circuit(vqubit.X(q0), vqubit.Y(q0), vqubit.Z(q0), vqubit.S(q0), vqubit.T(q1), vqubit.SWAP(q1, q2)),
		vqubit.Circuit(vqubit.rx(0.3)(q0), vqubit.ry(1.2)(q1), vqubit.rz(-0.7)(q2)),
		vqubit.Circuit(vqubit.H(q0), vqubit.CZ(q0, q1), vqubit.CCX(q0, q1, q2), vqubit.measure(q2, key="out")),
		vqubit.Circuit(vqubit.I(q0), vqubit.CNOT(q1, q2), vqubit.measure(q0, q1, q2, key="all")),
	]
	for circuit in corpus:
		text = vqubit.qasm(circuit)
		imported = vqubit.circuit_from_qasm(text)
		assert vqubit.qasm(imported) == text
		assert list(imported.all_qubits()) == list(circuit.all_qubits())
		assert_matrix_close(
			vqubit.unitary(vqubit.Circuit(*[op for op in imported.all_operations() if not vqubit.is_measurement(op)])),
			vqubit.unitary(vqubit.Circuit(*[op for op in circuit.all_operations() if not vqubit.is_measurement(op)])),
		)


def test_a_malformed_program_says_where() -> None:
	with pytest.raises(vqubit.QasmParseError) as error:
		vqubit.circuit_from_qasm(HEADER + "qreg q[2];\nh q[5];\n")
	assert error.value.line == 4
	assert "out of range" in str(error.value)
	with pytest.raises(vqubit.QasmParseError):
		vqubit.circuit_from_qasm(HEADER + "qreg q[1];\nnot_a_gate q[0];\n")
	with pytest.raises(vqubit.QasmParseError):
		vqubit.circuit_from_qasm(HEADER + "qreg q[1];\nh q[0]\n")
	with pytest.raises(vqubit.QasmParseError):
		vqubit.circuit_from_qasm(HEADER + "opaque mystery q;\nmystery x;\nqreg x[1];\n")
	with pytest.raises(TypeError):
		vqubit.circuit_from_qasm(42)  # ty: ignore[invalid-argument-type]
