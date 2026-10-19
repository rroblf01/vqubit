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


def test_a_wide_symbol_sits_in_a_column_of_its_own() -> None:
	q0 = vqubit.LineQubit(0)
	assert str(vqubit.Circuit(vqubit.Rx(0.3)(q0))) == "0: ───Rx───"


def test_the_plain_gates_draw_without_an_exponent() -> None:
	# ``vq.CZ`` *is* a ``CZPowGate`` at exponent 1, and Cirq draws that as a bare
	# ``@``; only a fractional power gets the ``^`` annotation.
	q0, q1 = vqubit.LineQubit.range(2)
	assert str(vqubit.Circuit(vqubit.CZ(q0, q1))) == "0: ───@───\n      │\n1: ───@───"
	assert str(vqubit.Circuit(vqubit.CNOT(q0, q1))) == "0: ───@───\n      │\n1: ───X───"
	assert str(vqubit.Circuit(vqubit.ISWAP(q0, q1))) == "0: ───iSwap───\n      │\n1: ───iSwap───"


def test_a_power_gate_annotates_the_last_symbol() -> None:
	# The exponent rides on the last symbol only, as in Cirq, and the first symbol
	# stays left aligned: the rest of the column is wire, not padding.
	q0, q1 = vqubit.LineQubit.range(2)
	assert str(vqubit.Circuit((vqubit.CZ**0.5)(q0, q1))) == "0: ───@───────\n      │\n1: ───@^0.5───"
	assert str(vqubit.Circuit((vqubit.CNOT**0.5)(q0, q1))) == "0: ───@───────\n      │\n1: ───X^0.5───"
	assert str(vqubit.Circuit((vqubit.SWAP**0.5)(q0, q1))) == "0: ───×───────\n      │\n1: ───×^0.5───"
	assert str(vqubit.Circuit((vqubit.ISWAP**0.5)(q0, q1))) == "0: ───iSwap───────\n      │\n1: ───iSwap^0.5───"


def test_a_single_qubit_power_gate_annotates_its_symbol() -> None:
	q0 = vqubit.LineQubit(0)
	assert str(vqubit.Circuit((vqubit.X**0.5)(q0))) == "0: ───X^0.5───"
	assert str(vqubit.Circuit((vqubit.XX**0.5)(q0, vqubit.LineQubit(1)))) == "0: ───XX───────\n      │\n1: ───XX^0.5───"


def test_swap_diagram_on_non_adjacent_qubits() -> None:
	q0, q2 = vqubit.LineQubit(0), vqubit.LineQubit(2)
	assert str(vqubit.Circuit(vqubit.SWAP(q0, q2))) == "0: ───×───\n      │\n2: ───×───"
