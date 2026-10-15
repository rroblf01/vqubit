"""Tests for the KAK decomposition and two-qubit synthesis."""

from __future__ import annotations

import math
import random

import pytest

import vqubit
from tests.helpers import random_unitary, up_to_phase

_QUARTER_PI = math.pi / 4


def _circuit_unitary(operations: list[vqubit.Operation]) -> list[list[complex]]:
	return vqubit.unitary(vqubit.Circuit(list(operations)))


def test_kak_angles_of_known_gates() -> None:
	cnot = vqubit.kak_decomposition(vqubit.unitary(vqubit.CNOT))
	assert sorted(abs(a) for a in cnot.interaction_coefficients) == pytest.approx([0, 0, _QUARTER_PI], abs=1e-6)
	swap = vqubit.kak_decomposition(vqubit.unitary(vqubit.SWAP))
	assert sorted(abs(a) for a in swap.interaction_coefficients) == pytest.approx(
		[_QUARTER_PI, _QUARTER_PI, _QUARTER_PI], abs=1e-6
	)


@pytest.mark.parametrize(
	"gate",
	[vqubit.CNOT, vqubit.CZ, vqubit.SWAP, vqubit.ISWAP, vqubit.XX**0.5, vqubit.XX, vqubit.YY, vqubit.ZZ],
)
def test_two_qubit_synthesis_of_known_gates(gate: vqubit.Gate) -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	matrix = vqubit.unitary(gate)
	cnot_ops = vqubit.two_qubit_matrix_to_operations(matrix, qubits=(q0, q1))
	cz_ops = vqubit.two_qubit_matrix_to_cz_operations(q0, q1, matrix)
	assert up_to_phase(_circuit_unitary(cnot_ops), matrix, tol=1e-5)
	assert up_to_phase(_circuit_unitary(cz_ops), matrix, tol=1e-5)


def test_two_qubit_synthesis_of_random_unitaries() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(1234)
	for _ in range(40):
		matrix = random_unitary(4, rng)
		cnot_ops = vqubit.two_qubit_matrix_to_operations(matrix, qubits=(q0, q1))
		cz_ops = vqubit.two_qubit_matrix_to_cz_operations(q0, q1, matrix)
		assert up_to_phase(_circuit_unitary(cnot_ops), matrix, tol=1e-4)
		assert up_to_phase(_circuit_unitary(cz_ops), matrix, tol=1e-4)


def test_cz_synthesis_uses_only_czs_and_single_qubit_gates() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(7)
	operations = vqubit.two_qubit_matrix_to_cz_operations(q0, q1, random_unitary(4, rng))
	for operation in operations:
		assert len(operation.qubits) == 1 or isinstance(operation.gate, vqubit.CZGate)


def test_target_gateset_synthesizes_an_arbitrary_two_qubit_gate() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rng = random.Random(99)
	matrix = random_unitary(4, rng)
	circuit = vqubit.Circuit(vqubit.MatrixGate(matrix).on(q0, q1))
	compiled = vqubit.optimize_for_target_gateset(circuit, vqubit.CZTargetGateset())
	for operation in compiled.all_operations():
		assert isinstance(operation.gate, (vqubit.CZGate, vqubit.PhasedXZGate))
	assert up_to_phase(vqubit.unitary(compiled), matrix, tol=1e-4)
