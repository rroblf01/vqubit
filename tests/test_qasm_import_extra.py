"""Extra tests for the OpenQASM 2.0 importer.

``tests/test_qasm_import.py`` covers the common path; this file documents the
rest of the grammar the parser accepts (the whole supported gate set, custom
gates, ``reset``/``barrier``, classical assignment and conditions) and the
diagnostics for malformed programs.
"""

from __future__ import annotations

import cmath
import math

import pytest

import vqubit
from tests.helpers import up_to_phase

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def _parse(body: str) -> vqubit.Circuit:
	"""Parse ``body`` after the usual header."""
	return vqubit.circuit_from_qasm(HEADER + body)


def _u3_matrix(theta: float, phi: float, lam: float) -> list[list[complex]]:
	"""The exact ``u3`` matrix, global phase included."""
	cosine, sine = math.cos(theta / 2), math.sin(theta / 2)
	return [
		[complex(cosine), -cmath.exp(1j * lam) * sine],
		[cmath.exp(1j * phi) * sine, cmath.exp(1j * (phi + lam)) * cosine],
	]


def test_programs_round_trip_through_the_exporter() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	corpus = [
		vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="c")),
		vqubit.Circuit(vqubit.rx(0.125)(q0), vqubit.ry(-0.5)(q1), vqubit.rz(math.pi / 8)(q2)),
		vqubit.Circuit(vqubit.S(q0), vqubit.T(q1), vqubit.SWAP(q0, q2), vqubit.CZ(q1, q2)),
		vqubit.Circuit(vqubit.X(q0), vqubit.CCX(q0, q1, q2), vqubit.measure(q2, key="out")),
	]
	for circuit in corpus:
		text = vqubit.qasm(circuit)
		again = vqubit.circuit_from_qasm(text)
		assert vqubit.qasm(again) == text
		assert list(again.all_qubits()) == list(circuit.all_qubits())


def test_adjacent_measurements_into_different_registers_round_trip() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="a"), vqubit.measure(q1, key="b"))
	vqubit.circuit_from_qasm(vqubit.qasm(circuit))


def test_every_plain_supported_gate_maps_to_the_same_gate() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	cases: list[tuple[str, vqubit.Operation]] = [
		("id q[0];", vqubit.I(q0)),
		("x q[0];", vqubit.X(q0)),
		("y q[0];", vqubit.Y(q0)),
		("z q[0];", vqubit.Z(q0)),
		("h q[0];", vqubit.H(q0)),
		("s q[0];", (vqubit.Z**0.5)(q0)),
		("sdg q[0];", (vqubit.Z**-0.5)(q0)),
		("t q[0];", (vqubit.Z**0.25)(q0)),
		("tdg q[0];", (vqubit.Z**-0.25)(q0)),
		("sx q[0];", (vqubit.X**0.5)(q0)),
		("sxdg q[0];", (vqubit.X**-0.5)(q0)),
		("cx q[0],q[1];", vqubit.CNOT(q0, q1)),
		("CX q[0],q[1];", vqubit.CNOT(q0, q1)),
		("cy q[0],q[1];", vqubit.ControlledGate(vqubit.Y, num_controls=1).on(q0, q1)),
		("cz q[0],q[1];", vqubit.CZ(q0, q1)),
		("ch q[0],q[1];", vqubit.ControlledGate(vqubit.H, num_controls=1).on(q0, q1)),
		("swap q[0],q[1];", vqubit.SWAP(q0, q1)),
		("ccx q[0],q[1],q[2];", vqubit.CCX(q0, q1, q2)),
	]
	for instruction, expected in cases:
		circuit = _parse("qreg q[3];\n" + instruction + "\n")
		assert list(circuit.all_operations()) == [expected], instruction


@pytest.mark.parametrize(
	("instruction", "matrix"),
	[
		("u1(0.4) q[0];", _u3_matrix(0.0, 0.0, 0.4)),
		("u2(0.2, 0.3) q[0];", _u3_matrix(math.pi / 2, 0.2, 0.3)),
		("u3(0.4, 0.5, 0.6) q[0];", _u3_matrix(0.4, 0.5, 0.6)),
		("U(0.4, 0.5, 0.6) q[0];", _u3_matrix(0.4, 0.5, 0.6)),
	],
)
def test_the_u_family_matches_the_qasm_matrix_up_to_a_global_phase(
	instruction: str, matrix: list[list[complex]]
) -> None:
	circuit = _parse("qreg q[1];\n" + instruction + "\n")
	assert up_to_phase(vqubit.unitary(circuit), matrix)


def test_u0_is_the_identity() -> None:
	circuit = _parse("qreg q[1];\nu0(0.9) q[0];\n")
	assert list(circuit.all_operations()) == [vqubit.I(vqubit.LineQubit(0))]


@pytest.mark.parametrize(
	("expression", "expected"),
	[
		("pi", math.pi),
		("-pi/2", -math.pi / 2),
		("+pi/4", math.pi / 4),
		("2*pi/3", 2 * math.pi / 3),
		("(pi/2)", math.pi / 2),
		("1+2", 3),
		("5-2", 3),
		("2^3", 8),
		("2^-1", 0.5),
		("6/2", 3),
		("2*-3", -6),
		(".5", 0.5),
		("1e-3", 1e-3),
	],
)
def test_angle_expressions(expression: str, expected: float) -> None:
	circuit = _parse(f"qreg q[1];\nrx({expression}) q[0];\n")
	gate = next(circuit.all_operations()).gate
	assert isinstance(gate, vqubit.RxGate)
	assert gate.radians == pytest.approx(expected)


def test_a_free_name_in_an_expression_is_a_parameter() -> None:
	circuit = _parse("qreg q[1];\nrx((theta + 1) * 2) q[0];\n")
	gate = next(circuit.all_operations()).gate
	assert isinstance(gate, vqubit.RxGate)
	assert str(gate.radians) == "((theta + 1) * 2)"
	# Symbolic circuits export and re-import unchanged.
	assert vqubit.qasm(vqubit.circuit_from_qasm(vqubit.qasm(circuit))) == vqubit.qasm(circuit)


def test_custom_gates_expand_inline_with_their_parameters() -> None:
	circuit = _parse(
		"""gate pair(a) q0, q1 {
  h q0;
  rz(2*a + pi/4) q1;
  cx q0, q1;
}
qreg q[2];
pair(pi/8) q[0], q[1];
"""
	)
	q0, q1 = vqubit.LineQubit.range(2)
	assert list(circuit.all_operations()) == [
		vqubit.H(q0),
		vqubit.Rz(2 * math.pi / 8 + math.pi / 4)(q1),
		vqubit.CNOT(q0, q1),
	]


def test_custom_gates_can_call_other_custom_gates() -> None:
	circuit = _parse(
		"""gate inner(a) q { rz(a) q; }
gate outer(b) q { inner(b/2) q; h q; }
qreg q[1];
outer(pi) q[0];
"""
	)
	q0 = vqubit.LineQubit(0)
	assert list(circuit.all_operations()) == [vqubit.Rz(math.pi / 2)(q0), vqubit.H(q0)]


def test_custom_gates_expand_over_register_pairs() -> None:
	circuit = _parse(
		"""gate pair(a) x, y { h x; rz(a) y; }
qreg a[2];
qreg b[2];
pair(0.5) a, b;
"""
	)
	operations = list(circuit.all_operations())
	assert [operation.qubits for operation in operations] == [
		(vqubit.LineQubit(0),),
		(vqubit.LineQubit(2),),
		(vqubit.LineQubit(1),),
		(vqubit.LineQubit(3),),
	]
	assert [operation.gate for operation in operations[1::2]] == [vqubit.Rz(0.5), vqubit.Rz(0.5)]


def test_custom_gate_argument_counts_are_checked() -> None:
	with pytest.raises(vqubit.QasmParseError, match=r"takes 1 parameter\(s\), got 2") as error:
		_parse("gate g(a) q { rz(a) q; }\nqreg r[1];\ng(1, 2) r[0];\n")
	assert error.value.line == 5
	with pytest.raises(vqubit.QasmParseError, match=r"expected 1 qubit argument\(s\), got 2"):
		_parse("gate g(a) q { rz(a) q; }\nqreg r[1];\ng(1) r[0], r[0];\n")


def test_measuring_a_register_produces_one_measurement_gate() -> None:
	circuit = _parse("qreg q[3];\ncreg c[3];\nmeasure q -> c;\n")
	operation = next(circuit.all_operations())
	assert operation.gate == vqubit.MeasurementGate("c", num_qubits=3)
	assert operation.qubits == vqubit.LineQubit.range(3)


def test_measurements_fuse_until_another_gate_interrupts_them() -> None:
	circuit = _parse(
		"""qreg q[2];
creg c[2];
measure q[0] -> c[0];
x q[1];
measure q[1] -> c[1];
"""
	)
	operations = list(circuit.all_operations())
	assert len(operations) == 3
	assert operations[0].gate == vqubit.MeasurementGate("c", num_qubits=1)
	assert operations[1] == vqubit.X(vqubit.LineQubit(1))
	assert operations[2].gate == vqubit.MeasurementGate("c", num_qubits=1)


def test_measurements_into_different_registers_are_separate_operations() -> None:
	# Consecutive measurements into different registers become one measurement
	# per register, in the order they were written.
	circuit = _parse("qreg q[1];\ncreg c[1];\ncreg d[1];\nmeasure q[0] -> c[0];\nmeasure q[0] -> d[0];\n")
	keys = [vqubit.measurement_key_name(operation) for operation in circuit.all_operations()]
	assert keys == ["c", "d"]


def test_classical_assignment_measurement_syntax() -> None:
	circuit = _parse("qreg q[2];\ncreg c[2];\nc = measure q;\n")
	operation = next(circuit.all_operations())
	assert operation.gate == vqubit.MeasurementGate("c", num_qubits=2)
	assert operation.qubits == vqubit.LineQubit.range(2)


def test_classical_assignment_can_target_one_bit() -> None:
	_parse("qreg q[1];\ncreg c[1];\nc[0] = measure q[0];\n")


def test_reset_and_barrier_accept_whole_registers() -> None:
	circuit = _parse("qreg q[3];\nreset q;\nbarrier q;\nbarrier;\n")
	operations = list(circuit.all_operations())
	assert len(operations) == 3
	assert [operation.qubits for operation in operations] == [
		(vqubit.LineQubit(0),),
		(vqubit.LineQubit(1),),
		(vqubit.LineQubit(2),),
	]
	assert all(operation.gate == vqubit.ResetChannel() for operation in operations)


def test_a_condition_on_one_classical_bit() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = _parse("qreg q[1];\ncreg c[2];\nif(c[0]==1) x q[0];\n")
	operation = next(circuit.all_operations())
	assert isinstance(operation, vqubit.ClassicallyControlledOperation)
	assert operation.conditions == (vqubit.KeyCondition("c", 0),)
	assert operation.sub_operation == vqubit.X(q0)


def test_a_zero_condition_is_written_as_its_negation() -> None:
	circuit = _parse("qreg q[1];\ncreg c[2];\nif(c[0]==0) x q[0];\n")
	operation = next(circuit.all_operations())
	assert isinstance(operation, vqubit.ClassicallyControlledOperation)
	assert operation.conditions == (~vqubit.KeyCondition("c", 0),)


def test_a_condition_on_a_whole_register_checks_every_bit() -> None:
	circuit = _parse("qreg q[1];\ncreg c[2];\nif(c==2) x q[0];\n")
	operation = next(circuit.all_operations())
	assert isinstance(operation, vqubit.ClassicallyControlledOperation)
	assert operation.conditions == (~vqubit.KeyCondition("c", 0), vqubit.KeyCondition("c", 1))


def test_conditions_can_be_combined_with_and() -> None:
	circuit = _parse("qreg q[1];\ncreg c[2];\ncreg d[2];\nif(c==2 && d[1]==0) x q[0];\n")
	operation = next(circuit.all_operations())
	assert isinstance(operation, vqubit.ClassicallyControlledOperation)
	assert operation.conditions == (
		~vqubit.KeyCondition("c", 0),
		vqubit.KeyCondition("c", 1),
		~vqubit.KeyCondition("d", 1),
	)


def test_a_condition_can_guard_a_multi_qubit_gate() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = _parse("qreg q[2];\ncreg c[1];\nif(c==1) cx q[0],q[1];\n")
	operations = list(circuit.all_operations())
	assert len(operations) == 1
	assert isinstance(operations[0], vqubit.ClassicallyControlledOperation)
	assert operations[0].sub_operation == vqubit.CNOT(q0, q1)


def test_qubits_are_numbered_across_registers_in_declaration_order() -> None:
	circuit = _parse("qreg a[2];\nqreg b[2];\nh a;\nh b[0];\ncx a[0], b[1];\n")
	assert set(circuit.all_qubits()) == set(vqubit.LineQubit.range(4))
	operation = list(circuit.all_operations())[-1]
	assert operation.qubits == (vqubit.LineQubit(0), vqubit.LineQubit(3))


def test_comments_and_blank_lines_can_appear_anywhere() -> None:
	circuit = vqubit.circuit_from_qasm(
		"""// a leading comment

OPENQASM 2.0;
/* a block
   comment */
include "qelib1.inc";
qreg q[1];

h /* between the name and its arguments */ q[0]; // trailing
"""
	)
	assert list(circuit.all_operations()) == [vqubit.H(vqubit.LineQubit(0))]


def test_an_empty_program_is_an_empty_circuit() -> None:
	# Nothing in the grammar requires a header, so an empty source is legal.
	assert list(vqubit.circuit_from_qasm("").all_operations()) == []
	assert list(vqubit.circuit_from_qasm(HEADER).all_operations()) == []


def test_opaque_declarations_parse_but_cannot_be_applied() -> None:
	assert list(_parse("opaque foo(a) q;\nopaque bar q;\nqreg q[1];\n").all_operations()) == []
	with pytest.raises(vqubit.QasmParseError, match="is opaque and cannot be expanded") as error:
		_parse("opaque foo(a) q;\nqreg q[1];\nfoo(0.5) q[0];\n")
	assert error.value.line == 5


def test_a_gate_definition_without_a_body_is_rejected() -> None:
	with pytest.raises(vqubit.QasmParseError, match=r"expected '\{'"):
		_parse("gate g q junk\n")


def test_the_importer_needs_text_not_a_circuit() -> None:
	with pytest.raises(TypeError, match="needs the QASM text"):
		vqubit.circuit_from_qasm(42)  # ty: ignore[invalid-argument-type]


def test_an_unterminated_gate_definition_is_rejected() -> None:
	with pytest.raises(vqubit.QasmParseError, match="unterminated definition of gate 'g'") as error:
		_parse("gate g q { h q;\n")
	assert error.value.line == 3


def test_an_unterminated_opaque_parenthesis_is_rejected() -> None:
	with pytest.raises(vqubit.QasmParseError, match="unterminated parenthesis") as error:
		_parse("opaque g(\n")
	assert error.value.line == 3


def test_an_unterminated_barrier_is_rejected() -> None:
	with pytest.raises(vqubit.QasmParseError, match="unterminated barrier"):
		_parse("qreg q[1];\nbarrier q[0]")


def test_an_unknown_gate_is_reported_with_its_line() -> None:
	with pytest.raises(vqubit.QasmParseError, match="unknown gate 'toffoli'") as error:
		_parse("qreg q[3];\ntoffoli q[0], q[1], q[2];\n")
	assert error.value.line == 4


def test_cswap_is_a_supported_qelib1_gate() -> None:
	_parse("qreg q[3];\ncswap q[0], q[1], q[2];\n")


def test_a_missing_semicolon_is_reported_at_the_end_of_the_statement() -> None:
	with pytest.raises(vqubit.QasmParseError, match=r"expected ';'") as error:
		_parse("qreg q[1];\nh q[0]")
	assert error.value.line == 4


def test_unbalanced_parentheses_in_a_gate_call_are_rejected() -> None:
	with pytest.raises(vqubit.QasmParseError, match=r"expected '\)'"):
		_parse("qreg q[1];\nrx(0.5 q[0];\n")


@pytest.mark.parametrize(
	("body", "message", "line"),
	[
		("qreg[1];", "expected a name", 3),
		("qreg 3[1];", "expected a name", 3),
		("qreg q[two];", "expected a number", 3),
		("qreg q[1.5];", "expected an integer", 3),
		("qreg q[1];\n;\n", "unexpected token ';'", 4),
	],
)
def test_declaration_and_statement_syntax_errors(body: str, message: str, line: int) -> None:
	with pytest.raises(vqubit.QasmParseError, match=message) as error:
		_parse(body)
	assert error.value.line == line


def test_duplicate_registers_are_rejected() -> None:
	with pytest.raises(vqubit.QasmParseError, match="quantum register 'q' is already declared") as error:
		_parse("qreg q[1];\nqreg q[1];\n")
	assert error.value.line == 4
	with pytest.raises(vqubit.QasmParseError, match="classical register 'c' is already declared") as error:
		_parse("creg c[1];\ncreg c[1];\n")
	assert error.value.line == 4


def test_undeclared_registers_and_out_of_range_indices_are_rejected() -> None:
	with pytest.raises(vqubit.QasmParseError, match="undefined quantum register 'r'") as error:
		_parse("qreg q[1];\nh r[0];\n")
	assert error.value.line == 4
	with pytest.raises(vqubit.QasmParseError, match="undefined classical register 'c'"):
		_parse("qreg q[1];\nmeasure q[0] -> c[0];\n")
	with pytest.raises(vqubit.QasmParseError, match=r"q\[2\] is out of range"):
		_parse("qreg q[1];\nh q[2];\n")
	with pytest.raises(vqubit.QasmParseError, match=r"c\[1\] is out of range"):
		_parse("qreg q[1];\ncreg c[1];\nmeasure q[0] -> c[1];\n")


def test_wrong_qubit_and_parameter_counts_are_rejected() -> None:
	with pytest.raises(vqubit.QasmParseError, match=r"expected 2 qubit argument\(s\), got 1"):
		_parse("qreg q[1];\ncx q[0];\n")
	with pytest.raises(vqubit.QasmParseError, match=r"expected 1 qubit argument\(s\), got 2"):
		_parse("qreg q[1];\nh q[0], q[0];\n")
	with pytest.raises(vqubit.QasmParseError, match="different register sizes"):
		_parse("qreg a[1];\nqreg b[2];\ncx a, b;\n")
	with pytest.raises(vqubit.QasmParseError, match=r"takes 3 parameter\(s\), got 2"):
		_parse("qreg q[1];\nu3(0.1, 0.2) q[0];\n")


def test_a_plain_gate_rejects_parameters() -> None:
	with pytest.raises(vqubit.QasmParseError, match="gate 'h' takes no parameters") as error:
		_parse("qreg q[1];\nh(0.4) q[0];\n")
	assert error.value.line == 4


def test_unsupported_expressions_and_tokens_are_rejected() -> None:
	with pytest.raises(vqubit.QasmParseError, match="in an expression"):
		_parse("qreg q[1];\nrx() q[0];\n")
	with pytest.raises(vqubit.QasmParseError, match="unexpected character '@'") as error:
		_parse("qreg q[1];\nrx(@) q[0];\n")
	assert error.value.line == 4


def test_the_header_and_include_are_validated() -> None:
	with pytest.raises(vqubit.QasmParseError, match=r"expected OpenQASM version 2\.0") as error:
		vqubit.circuit_from_qasm("OPENQASM 3.0;\n")
	assert error.value.line == 1
	with pytest.raises(vqubit.QasmParseError, match="expected a file name after include") as error:
		vqubit.circuit_from_qasm("include 123;\n")
	assert error.value.line == 1


def test_a_condition_needs_a_declared_classical_register() -> None:
	with pytest.raises(vqubit.QasmParseError, match="undefined classical register 'c'") as error:
		_parse("qreg q[1];\nif(c==1) x q[0];\n")
	assert error.value.line == 4


@pytest.mark.parametrize(
	"instruction",
	[
		"crz(0.3) q[0], q[1];",
		"cu1(0.3) q[0], q[1];",
		"cu3(0.1, 0.2, 0.3) q[0], q[1];",
	],
)
def test_the_controlled_u_gates_accept_two_qubits(instruction: str) -> None:
	_parse("qreg q[2];\n" + instruction + "\n")
