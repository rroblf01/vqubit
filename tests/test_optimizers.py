"""Tests for the optimizers."""

from __future__ import annotations

import pytest

import vqubit
from tests.helpers import assert_matrix_close


def _up_to_phase(actual: list[list[complex]], expected: list[list[complex]]) -> bool:
	phase = None
	for j in range(len(expected[0])):
		if abs(expected[0][j]) > 1e-9:
			phase = actual[0][j] / expected[0][j]
			break
	if phase is None:
		return False
	return all(
		abs(actual[i][j] - phase * expected[i][j]) < 1e-9 for i in range(len(expected)) for j in range(len(expected))
	)


@pytest.mark.parametrize(
	"gate",
	[vqubit.X, vqubit.Y, vqubit.Z, vqubit.H, vqubit.S, vqubit.T, vqubit.Rx(0.3), vqubit.Ry(1.2), vqubit.Rz(-0.4)],
)
def test_single_qubit_decomposition(gate: vqubit.Gate) -> None:
	q0 = vqubit.LineQubit(0)
	gates = vqubit.single_qubit_matrix_to_gates(vqubit.unitary(gate))
	assert all(isinstance(g, (vqubit.RzGate, vqubit.RyGate)) for g in gates)
	reconstructed = vqubit.unitary(vqubit.Circuit([g.on(q0) for g in gates]))
	assert _up_to_phase(reconstructed, vqubit.unitary(gate))


def test_single_qubit_decomposition_of_identity() -> None:
	q0 = vqubit.LineQubit(0)
	gates = vqubit.single_qubit_matrix_to_gates(vqubit.unitary(vqubit.I))
	reconstructed = vqubit.unitary(vqubit.Circuit([g.on(q0) for g in gates]))
	assert_matrix_close(reconstructed, vqubit.unitary(vqubit.I))


def test_single_qubit_matrix_to_operations() -> None:
	q0 = vqubit.LineQubit(0)
	operations = vqubit.single_qubit_matrix_to_operations(vqubit.unitary(vqubit.H), q0)
	assert all(op.qubits == (q0,) for op in operations)
	reconstructed = vqubit.unitary(vqubit.Circuit(list(operations)))
	assert _up_to_phase(reconstructed, vqubit.unitary(vqubit.H))
