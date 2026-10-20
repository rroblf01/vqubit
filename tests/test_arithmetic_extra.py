"""Tests for the Boolean Hamiltonian gate and the arithmetic base class.

``BooleanHamiltonianGate`` is the diagonal gate a QAOA-style phase separator
needs: ``diag(exp(i theta/2 sum_k f_k(x)))`` over the bit strings ``x``. The
arithmetic base class turns a register permutation into a unitary; its
abstract pieces must refuse to run on their own.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import vqubit
from tests.helpers import assert_matrix_close


def test_arithmetic_gate_base_class_is_abstract() -> None:
	gate = vqubit.ArithmeticGate()
	with pytest.raises(NotImplementedError):
		gate.registers()
	with pytest.raises(NotImplementedError):
		gate.with_registers()
	with pytest.raises(NotImplementedError):
		gate.apply(0)


def test_boolean_hamiltonian_is_diagonal_and_unitary() -> None:
	gate = vqubit.BooleanHamiltonianGate(parameter_names=["a", "b"], boolean_strs=["a & b"], theta=0.5)
	assert gate.num_qubits == 2
	matrix = np.asarray(vqubit.unitary(gate))
	# Diagonal, with exp(i theta/2) on |11> and 1 elsewhere.
	assert_matrix_close(matrix, np.diag([1, 1, 1, complex(math.cos(0.25), math.sin(0.25))]))
	assert vqubit.has_unitary(gate)
	assert gate == vqubit.BooleanHamiltonianGate(parameter_names=["a", "b"], boolean_strs=["a & b"], theta=0.5)
	assert gate != vqubit.BooleanHamiltonianGate(parameter_names=["a", "b"], boolean_strs=["a | b"], theta=0.5)
	assert "BooleanHamiltonianGate" in str(gate)
	assert "boolean_strs" in repr(gate)


def test_boolean_expressions_support_the_operators() -> None:
	# One expression per operator, evaluated over the four bit strings.
	expressions = ["a & b", "a | b", "a ^ b", "~a"]
	gate = vqubit.BooleanHamiltonianGate(parameter_names=["a", "b"], boolean_strs=expressions, theta=0.4)
	diagonal = np.diag(np.asarray(vqubit.unitary(gate)))
	# |00>: only ~a; |01>: a|b, a^b and ~a; |10>: a|b and a^b; |11>: a&b and a|b.
	expected_counts = [1, 3, 2, 2]
	for index, count in enumerate(expected_counts):
		assert diagonal[index] == pytest.approx(complex(math.cos(0.2 * count), math.sin(0.2 * count))), index


def test_boolean_expressions_are_validated() -> None:
	with pytest.raises(ValueError, match="Unknown symbol"):
		vqubit.BooleanHamiltonianGate(parameter_names=["a"], boolean_strs=["c"], theta=0.5)
	with pytest.raises(ValueError, match="may only use"):
		vqubit.BooleanHamiltonianGate(parameter_names=["a"], boolean_strs=["a + 1"], theta=0.5)
