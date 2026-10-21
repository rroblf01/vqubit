"""Tests for tagged operations, circuit operations and the extra values."""

from __future__ import annotations

import cmath
import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close, up_to_phase
from vqubit._linalg import matmul


def test_tagged_operation() -> None:
	q0 = vqubit.LineQubit(0)
	tagged = vqubit.TaggedOperation(vqubit.H(q0), "hardware")
	assert tagged.tags == ("hardware",)
	assert tagged.qubits == (q0,)
	assert tagged.sub_operation == vqubit.H(q0)
	assert tagged.untagged == vqubit.H(q0)
	assert "hardware" in str(vqubit.Circuit(tagged))


def test_tagged_operation_validates() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(ValueError):
		vqubit.TaggedOperation(vqubit.H(q0))
	with pytest.raises(TypeError):
		vqubit.TaggedOperation("not an operation", "tag")  # ty: ignore[invalid-argument-type]


def test_circuit_operation() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	sub = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	operation = vqubit.CircuitOperation(sub, repetitions=2)
	assert operation.qubits == (q0, q1)
	assert vqubit.num_qubits(operation) == 2
	squared = matmul(sub.unitary(), sub.unitary())
	assert_matrix_close(vqubit.unitary(operation), squared)


def test_global_phase_gate() -> None:
	gate = vqubit.GlobalPhaseGate(0.5)
	assert gate.num_qubits == 0
	assert vqubit.unitary(gate)[0][0] == pytest.approx(cmath.exp(1j * math.pi * 0.5))
	assert vqubit.GlobalPhaseGate(0.5) == vqubit.global_phase(0.5)


def test_wait_gate() -> None:
	gate = vqubit.WaitGate(100.0, num_qubits=2)
	assert gate.num_qubits == 2
	assert_matrix_close(vqubit.unitary(gate), vqubit.unitary(vqubit.IdentityGate(2)))
	with pytest.raises(ValueError):
		vqubit.WaitGate(-1.0)


def test_qubit_order() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert vqubit.QubitOrder.DEFAULT.order_for([q1, q0]) == (q0, q1)
	assert vqubit.QubitOrder.explicit([q1, q0]).order_for([q0, q1]) == (q1, q0)


def test_dense_pauli_string() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	dense = vqubit.DensePauliString("XZ", coefficient=1j)
	assert str(dense) == "XZ"
	assert len(dense) == 2
	assert dense[0] is vqubit.Pauli.X
	assert dense.to_pauli_string([q0, q1]) == vqubit.PauliString(
		{q0: vqubit.Pauli.X, q1: vqubit.Pauli.Z}, coefficient=1j
	)
	assert dense * dense == vqubit.DensePauliString("II", coefficient=-1)
	assert dense.matrix()[0][0] == pytest.approx(0j)


def test_dense_pauli_string_identity_matrix() -> None:
	assert up_to_phase(vqubit.DensePauliString("X").matrix(), vqubit.unitary(vqubit.X))
