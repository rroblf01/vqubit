"""Tests for the text diagrams."""

from __future__ import annotations

import vqubit


def test_empty_circuit_has_no_diagram() -> None:
	assert str(vqubit.Circuit()) == ""


def test_single_gate_diagram() -> None:
	q0 = vqubit.LineQubit(0)
	assert str(vqubit.Circuit(vqubit.H(q0))) == "0: ───H───"


def test_independent_gates_share_a_line() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert str(vqubit.Circuit(vqubit.H(q0), vqubit.X(q1))) == "0: ───H───\n1: ───X───"


def test_two_qubit_gate_diagram() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	expected = "0: ───H───@───\n          │\n1: ───────X───"
	assert str(vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))) == expected


def test_measurement_diagram() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert str(vqubit.Circuit(vqubit.measure(q0, q1, key="m"))) == "0: ───M───\n      │\n1: ───M───"


def test_wide_symbol_is_centered_in_the_column() -> None:
	q0 = vqubit.LineQubit(0)
	assert str(vqubit.Circuit(vqubit.Rx(0.3)(q0))) == "0: ───Rx───"


def test_swap_diagram_on_non_adjacent_qubits() -> None:
	q0, q2 = vqubit.LineQubit(0), vqubit.LineQubit(2)
	assert str(vqubit.Circuit(vqubit.SWAP(q0, q2))) == "0: ───×───\n      │\n2: ───×───"
