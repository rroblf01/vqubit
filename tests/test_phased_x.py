"""Tests for the phased-X gates and the phased-X/Z optimizers."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close
from vqubit._linalg import matmul


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


def test_phased_x_is_conjugated_x() -> None:
	gate = vqubit.PhasedXPowGate(phase_exponent=0.25, exponent=1)
	expected = matmul(matmul(vqubit.unitary(vqubit.Z**0.25), vqubit.unitary(vqubit.X)), vqubit.unitary(vqubit.Z**-0.25))
	assert_matrix_close(vqubit.unitary(gate), expected)


def test_phased_x_with_half_phase_is_ry() -> None:
	gate = vqubit.PhasedXPowGate(phase_exponent=0.5, exponent=1)
	assert _up_to_phase(vqubit.unitary(gate), vqubit.unitary(vqubit.Ry(math.pi)))


def test_phased_x_power_and_inverse() -> None:
	gate = vqubit.PhasedXPowGate(phase_exponent=0.3, exponent=0.5)
	assert vqubit.pow(gate, 2) == vqubit.PhasedXPowGate(phase_exponent=0.3, exponent=1)
	inverse = vqubit.inverse(gate)
	assert _up_to_phase(matmul(vqubit.unitary(gate), vqubit.unitary(inverse)), [[1, 0], [0, 1]])


def test_phased_x_parameterization() -> None:
	t = vqubit.Parameter("t")
	gate = vqubit.PhasedXPowGate(phase_exponent=t)
	assert vqubit.is_parameterized(gate)
	resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({t: 0.5}))
	assert resolved == vqubit.PhasedXPowGate(phase_exponent=0.5)


def test_phased_xz_unitary() -> None:
	gate = vqubit.PhasedXZGate(x_exponent=0.5, z_exponent=0.25, axis_phase_exponent=0.1)
	expected = matmul(
		vqubit.unitary(vqubit.Z**0.25),
		vqubit.unitary(vqubit.PhasedXPowGate(phase_exponent=0.1, exponent=0.5)),
	)
	assert_matrix_close(vqubit.unitary(gate), expected)


def test_phased_xz_inverse() -> None:
	gate = vqubit.PhasedXZGate(x_exponent=0.3, z_exponent=0.2, axis_phase_exponent=0.1)
	inverse = vqubit.inverse(gate)
	assert _up_to_phase(matmul(vqubit.unitary(gate), vqubit.unitary(inverse)), [[1, 0], [0, 1]])


@pytest.mark.parametrize("gate", [vqubit.X, vqubit.H, vqubit.Rx(0.3), vqubit.Rz(-0.4), vqubit.S])
def test_single_qubit_matrix_to_phased_x_z(gate: vqubit.Gate) -> None:
	q0 = vqubit.LineQubit(0)
	gates = vqubit.single_qubit_matrix_to_phased_x_z(vqubit.unitary(gate))
	assert len(gates) == 2
	reconstructed = vqubit.unitary(vqubit.Circuit([g.on(q0) for g in gates]))
	assert _up_to_phase(reconstructed, vqubit.unitary(gate))


@pytest.mark.parametrize("gate", [vqubit.X, vqubit.H, vqubit.Rx(0.3), vqubit.Rz(-0.4), vqubit.S])
def test_single_qubit_matrix_to_phxz(gate: vqubit.Gate) -> None:
	q0 = vqubit.LineQubit(0)
	phxz = vqubit.single_qubit_matrix_to_phxz(vqubit.unitary(gate))
	assert isinstance(phxz, vqubit.PhasedXZGate)
	reconstructed = vqubit.unitary(vqubit.Circuit([phxz.on(q0)]))
	assert _up_to_phase(reconstructed, vqubit.unitary(gate))


def test_phased_gates_json_roundtrip() -> None:
	for gate in [
		vqubit.PhasedXPowGate(phase_exponent=0.25, exponent=0.5),
		vqubit.PhasedXZGate(x_exponent=0.3, z_exponent=0.2, axis_phase_exponent=0.1),
	]:
		assert vqubit.read_json(vqubit.to_json(gate)) == gate


def test_merge_single_qubit_moments_to_phxz() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.H(q0), vqubit.Rz(0.3)(q0))
	merged = vqubit.merge_single_qubit_moments_to_phxz(circuit)
	operations = list(merged.all_operations())
	assert len(operations) == 1
	assert isinstance(operations[0].gate, vqubit.PhasedXZGate)
	assert _up_to_phase(vqubit.unitary(merged), vqubit.unitary(circuit))


def test_merge_does_not_cross_two_qubit_gates() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.X(q0))
	merged = vqubit.merge_single_qubit_moments_to_phxz(circuit)
	# The CNOT on q0 splits the run, so q0 keeps two single-qubit gates.
	q0_gates = [op for op in merged.all_operations() if op.qubits == (q0,)]
	assert len(q0_gates) == 2


def test_merge_keeps_measurements() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.measure(q0, key="m"))
	merged = vqubit.merge_single_qubit_moments_to_phxz(circuit)
	assert vqubit.is_measurement(merged)
	operations = list(merged.all_operations())
	assert len(operations) == 2
	assert isinstance(operations[0].gate, vqubit.PhasedXZGate)
	assert operations[1].gate == vqubit.MeasurementGate("m")
