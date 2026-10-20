"""Tests for alignment, point optimizers, QASM output and the text drawer."""

from __future__ import annotations

from typing import Any, cast

import vqubit
from tests.helpers import assert_matrix_close


def test_alignment_has_the_cirq_strategies() -> None:
	assert [(member.name, member.value) for member in vqubit.Alignment] == [
		("EARLIEST", "earliest"),
		("NEW_THEN_INLINE", "new_then_inline"),
		("INLINE", "inline"),
		("APPEND", "append"),
		("LEFT", "left"),
		("RIGHT", "right"),
	]


def test_circuit_function_builds_a_circuit_from_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	calls: list[tuple[vqubit.Qid, ...]] = []

	def operations(qubits: tuple[vqubit.Qid, ...]) -> list[vqubit.Operation]:
		calls.append(qubits)
		return [vqubit.H(qubits[0]), vqubit.CNOT(qubits[0], qubits[1])]

	# CircuitFunction is annotated as taking a moment-producing function, but its
	# implementation also accepts operations; the cast records the actual behaviour.
	circuit = vqubit.CircuitFunction(cast(Any, operations), [q0, q1])

	# The function is handed exactly the declared qubits, once.
	assert calls == [(q0, q1)]
	assert circuit.qubits == (q0, q1)
	assert list(circuit.all_operations()) == [vqubit.H(q0), vqubit.CNOT(q0, q1)]


def test_circuit_function_packs_like_a_plain_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.CircuitFunction(cast(Any, lambda qubits: [vqubit.H(qubits[0]), vqubit.X(qubits[1])]), [q0, q1])

	# Operations on different qubits share a moment, exactly as in ``Circuit``.
	assert circuit == vqubit.Circuit(vqubit.H(q0), vqubit.X(q1))
	assert len(circuit) == 1


def test_circuit_function_qubits_are_declared_even_when_unused() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.CircuitFunction(cast(Any, lambda qubits: [vqubit.H(qubits[0])]), [q0, q1])

	assert circuit.qubits == (q0, q1)
	assert circuit.all_qubits() == (q0,)


def test_point_optimization_summary_is_a_value_object() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	summary = vqubit.PointOptimizationSummary(clear_span=2, clear_qubits=[q0, q1], new_operations=[vqubit.CNOT(q0, q1)])

	assert summary.clear_span == 2
	assert summary.clear_qubits == [q0, q1]
	assert summary.new_operations == [vqubit.CNOT(q0, q1)]
	assert summary == vqubit.PointOptimizationSummary(2, [q0, q1], [vqubit.CNOT(q0, q1)])
	assert summary != vqubit.PointOptimizationSummary(1, [q0, q1], [vqubit.CNOT(q0, q1)])


def test_point_optimizer_call_preserves_the_circuit_unitary() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.Z(q1))
	optimizer = vqubit.PointOptimizer()

	result = optimizer(circuit)

	assert isinstance(result, vqubit.Circuit)
	assert list(result.all_operations()) == list(circuit.all_operations())
	assert_matrix_close(vqubit.unitary(result), vqubit.unitary(circuit))


def test_point_optimizer_keeps_independent_operations_in_one_moment() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q1))
	result = vqubit.PointOptimizer().optimize_circuit(circuit)

	assert result == circuit
	assert len(result) == 1


def test_point_optimizer_default_hooks_do_nothing() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0))
	operation = vqubit.X(q0)
	optimizer = vqubit.PointOptimizer()

	assert optimizer.optimization_at(circuit, 0, operation) is None
	assert optimizer.convert(operation) == [operation]


def test_qasm_output_wraps_the_circuit_text() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	output = vqubit.QasmOutput(circuit)

	assert output.circuit is circuit
	assert output.qasm() == vqubit.qasm(circuit)
	assert str(output) == output.qasm()
	assert output.qasm().splitlines()[-1] == "cx q[0],q[1];"


def test_qasm_output_of_an_empty_circuit_still_declares_the_register() -> None:
	output = vqubit.QasmOutput(vqubit.Circuit())

	assert output.qasm().splitlines()[-1] == "qreg q[0];"


def test_text_diagram_drawer_starts_empty() -> None:
	drawer = vqubit.TextDiagramDrawer()

	assert drawer.width == 0
	assert drawer.height == 0
	assert drawer.render() == ""


def test_text_diagram_drawer_writes_cells_and_grows_its_extent() -> None:
	drawer = vqubit.TextDiagramDrawer()
	drawer.write(0, 0, "H")
	drawer.write(1, 0, "X")

	assert drawer.width == 2
	assert drawer.height == 1
	assert drawer.render() == " H X"


def test_text_diagram_drawer_renders_one_line_per_row() -> None:
	drawer = vqubit.TextDiagramDrawer()
	drawer.write(0, 0, "H")
	drawer.write(0, 1, "X")

	assert drawer.render() == " H\n X"


def test_text_diagram_drawer_overwrites_a_cell() -> None:
	drawer = vqubit.TextDiagramDrawer()
	drawer.write(0, 0, "H")
	drawer.write(0, 0, "Z")

	assert drawer.render() == " Z"
	assert drawer.width == 1


def test_text_diagram_drawer_draws_horizontal_wires_only_between_endpoints() -> None:
	drawer = vqubit.TextDiagramDrawer()
	drawer.write(0, 0, "a")
	drawer.write(2, 1, "b")
	drawer.horizontal_line(0, 3, 0)

	# The wire fills x=0,1,2 but not x=3, and leaves the second row connected to b.
	assert drawer.render() == " a - -\n     b"


def test_text_diagram_drawer_normalizes_reversed_horizontal_lines() -> None:
	q0 = vqubit.LineQubit(0)
	forward = vqubit.TextDiagramDrawer()
	forward.write(0, 0, str(q0))
	forward.horizontal_line(0, 2, 0)

	backward = vqubit.TextDiagramDrawer()
	backward.write(0, 0, str(q0))
	backward.horizontal_line(2, 0, 0)

	assert forward.render() == backward.render()


def test_text_diagram_drawer_vertical_lines_do_not_grow_the_grid() -> None:
	# Reserving a vertical wire between existing rows must not add rows or columns.
	drawer = vqubit.TextDiagramDrawer()
	drawer.write(0, 0, "H")
	drawer.write(0, 1, "X")
	drawer.vertical_line(0, 2, 0)

	assert (drawer.width, drawer.height) == (1, 2)


def test_text_diagram_drawer_spacing_adds_padding_and_blank_lines() -> None:
	drawer = vqubit.TextDiagramDrawer()
	drawer.write(0, 0, "H")
	drawer.write(0, 1, "X")

	assert drawer.render(horizontal_spacing=2) == "  H\n  X"
	assert drawer.render(vertical_spacing=1) == " H\n\n X"
	assert drawer.render(horizontal_spacing=2, vertical_spacing=1) == "  H\n\n  X"


def test_text_diagram_drawer_blank_cells_are_spaces() -> None:
	drawer = vqubit.TextDiagramDrawer()
	drawer.write(2, 0, "Z")

	assert drawer.width == 3
	assert drawer.render() == "     Z"
