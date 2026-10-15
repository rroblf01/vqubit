"""Tests for ArithmeticGate and BooleanHamiltonianGate."""

from __future__ import annotations

import cmath
from typing import Any

import pytest

import vqubit
from tests.helpers import assert_matrix_close


class _Add(vqubit.ArithmeticGate):
	def __init__(self, target_register: Any, input_register: Any) -> None:
		self.target_register = target_register
		self.input_register = input_register

	def registers(self) -> Any:
		return self.target_register, self.input_register

	def with_registers(self, *new_registers: Any) -> _Add:
		return _Add(*new_registers)

	def apply(self, *register_values: int) -> int:
		return sum(register_values)


def test_arithmetic_adder() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	unitary = vqubit.unitary(_Add([2, 2], 1).on(q0, q1))
	assert_matrix_close(unitary, [[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])


def test_arithmetic_gate_registers() -> None:
	gate = _Add([2, 2], 1)
	assert gate.num_qubits == 2
	assert gate.with_registers([2, 2], 2).registers() == ([2, 2], 2)


def test_boolean_hamiltonian_gate() -> None:
	gate = vqubit.BooleanHamiltonianGate(["x0", "x1"], ["x0 & x1"], 1.0)
	assert gate.num_qubits == 2
	unitary = vqubit.unitary(gate)
	assert vqubit.is_unitary(unitary)
	# Only |11> satisfies the expression, so only its phase differs from 1.
	assert unitary[0][0] == pytest.approx(1)
	assert abs(unitary[3][3] - cmath.exp(0.5j)) < 1e-9


def test_boolean_hamiltonian_rejects_bad_expression() -> None:
	with pytest.raises(ValueError):
		vqubit.BooleanHamiltonianGate(["x0"], ["__import__('os')"], 1.0)
