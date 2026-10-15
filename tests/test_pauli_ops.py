"""Tests for the Pauli-string operations and related helpers."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close, up_to_phase


def _bell_string() -> vqubit.PauliString:
	q0, q1 = vqubit.LineQubit.range(2)
	return vqubit.PauliString({q0: vqubit.Pauli.X, q1: vqubit.Pauli.Z})


def test_pauli_string_phasor() -> None:
	pauli_string = _bell_string()
	phasor = vqubit.PauliStringPhasor(pauli_string, half_turns=1.0)
	assert vqubit.is_unitary(vqubit.unitary(phasor))
	half = vqubit.PauliStringPhasor(pauli_string, half_turns=0.5)
	assert up_to_phase(vqubit.dot(vqubit.unitary(half), vqubit.unitary(half)), vqubit.unitary(phasor))
	assert phasor.gate.num_qubits == 2
	assert (phasor**0.5).half_turns == pytest.approx(0.5)


def test_pauli_string_gate_operation() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	pauli_string = _bell_string()
	operation = vqubit.PauliStringGateOperation(pauli_string, [q0, q1])
	assert operation.qubits == (q0, q1)
	assert operation.gate.num_qubits == 2
	assert up_to_phase(vqubit.unitary(operation), pauli_string.matrix())
	single = vqubit.SingleQubitPauliStringGateOperation(vqubit.PauliString({q0: vqubit.Pauli.X}), q0)
	assert single.gate == vqubit.X


def test_mutable_pauli_string() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	pauli_string = vqubit.MutablePauliString({q0: vqubit.Pauli.X})
	pauli_string[q1] = vqubit.Pauli.Z
	assert pauli_string[q1] == vqubit.Pauli.Z
	pauli_string.coefficient = 1j
	assert pauli_string.coefficient == 1j
	del pauli_string[q0]
	assert pauli_string[q0] == vqubit.Pauli.I
	assert str(pauli_string.to_pauli_string()) == "1j*Z(1)"


def test_mutable_dense_pauli_string() -> None:
	dense = vqubit.MutableDensePauliString("IX")
	dense[0] = vqubit.Pauli.Z
	assert str(dense) == "ZX"
	dense.coefficient = -1
	assert dense.coefficient == -1


def test_base_dense_pauli_string_alias() -> None:
	assert vqubit.BaseDensePauliString is vqubit.DensePauliString


def test_pauli_measurement_gate() -> None:
	pauli_string = _bell_string()
	gate = vqubit.PauliMeasurementGate(pauli_string, key="m")
	assert gate.num_qubits == 2
	assert gate._is_measurement_()
	assert gate._measurement_key_name_() == "m"
	operation = vqubit.measure_single_paulistring(pauli_string, key="k")
	assert operation.gate._measurement_key_name_() == "k"
	terms = vqubit.PauliSum({pauli_string: 1.0, vqubit.PauliString(): 0.5})
	measurements = vqubit.measure_paulistring_terms(terms, key="t")
	assert len(measurements) == 2


def test_pauli_sum_exponential() -> None:
	pauli_string = _bell_string()
	exponential = vqubit.PauliSumExponential(vqubit.PauliSum({pauli_string: 1.0}))
	assert vqubit.is_unitary(exponential.matrix())
	# exp(i π P) = cos(π) I + i sin(π) P = -I.
	assert_matrix_close(exponential.matrix(), [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])


def test_projector_string() -> None:
	pauli_string = _bell_string()
	projector = vqubit.ProjectorString(pauli_string)
	matrix = projector.matrix()
	assert matrix[0][0] == pytest.approx(1)
	projector_sum = vqubit.ProjectorSum({pauli_string: 1.0})
	assert len(projector_sum) == 1


def test_linear_combinations() -> None:
	q0 = vqubit.LineQubit(0)
	combination = vqubit.LinearCombinationOfGates({vqubit.X: 1.0, vqubit.Z: 1.0})
	matrix = combination.matrix()
	assert_matrix_close(matrix, [[1, 1], [1, -1]])
	operations = vqubit.LinearCombinationOfOperations({vqubit.X(q0): 1.0, vqubit.Z(q0): 1.0})
	assert_matrix_close(operations.matrix(), [[1, 1], [1, -1]])
