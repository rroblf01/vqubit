"""Tests for the Pauli-string operations behind ``vqubit.ops``.

``PauliString`` is immutable; ``MutablePauliString`` is the in-place editing
form. Both turn into operations: a gate operation that applies the string, a
phasor ``exp(i π half_turns P)``, and a measurement along the string's axes.
"""

from __future__ import annotations

import cmath
import math

import numpy as np
import pytest

import vqubit
from tests.helpers import assert_matrix_close

P = vqubit.Pauli


def test_mutable_pauli_string_edits_in_place() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	mutable = vqubit.MutablePauliString(coefficient=2.0)
	mutable[q0] = P.X
	mutable[q1] = P.Y
	assert mutable.qubits == (q0, q1)
	assert mutable[q0] is P.X
	assert mutable.coefficient == 2.0
	assert list(mutable.items()) == [(q0, P.X), (q1, P.Y)]
	# Setting I removes the qubit: the string no longer touches it.
	mutable[q1] = P.I
	assert mutable.qubits == (q0,)
	assert mutable[q1] is P.I
	# Deleting removes it too, and deleting an absent qubit is an error.
	del mutable[q0]
	assert mutable.qubits == ()
	with pytest.raises(KeyError):
		del mutable[q0]


def test_mutable_pauli_string_validation() -> None:
	q0 = vqubit.LineQubit(0)
	mutable = vqubit.MutablePauliString()
	with pytest.raises(TypeError, match="must be Pauli"):
		mutable[q0] = vqubit.X  # ty: ignore[invalid-assignment] - a gate is not a Pauli
	mutable[q0] = P.X
	mutable.coefficient = 3
	assert mutable.coefficient == 3
	assert mutable.to_pauli_string() == vqubit.PauliString({q0: P.X}, 3)
	assert str(mutable) == "3*X(0)"
	assert repr(mutable) == "vqubit.MutablePauliString(3*X(0))"


def test_mutable_dense_pauli_string_edits_in_place() -> None:
	dense = vqubit.MutableDensePauliString("XY", coefficient=2)
	assert dense.coefficient == 2
	dense[0] = P.Z
	dense.coefficient = 1.5
	assert str(dense) == "ZY"
	assert dense.coefficient == 1.5
	assert dense == vqubit.MutableDensePauliString("ZY", 1.5)
	# A mutable dense string is not equal to its immutable sibling, as in Cirq.
	assert dense != vqubit.DensePauliString("ZY", 1.5)


def test_pauli_string_gate_operation_is_the_pauli() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	string = vqubit.PauliString({q0: P.X, q1: P.Y})
	operation = vqubit.PauliStringGateOperation(string, [q0, q1])
	assert operation.pauli_string == string
	assert operation.qubits == (q0, q1)
	assert_matrix_close(vqubit.unitary(operation), np.kron(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.Y)))
	# The inverse of a Pauli string is the string itself, and its square is the
	# identity on the same qubits.
	inverse = operation**-1
	assert isinstance(inverse, vqubit.PauliStringGateOperation)
	assert inverse.pauli_string == string
	assert_matrix_close(vqubit.unitary(operation**2), np.eye(4))
	assert (operation**2).qubits == (q0, q1)
	# A remapped copy acts on the new qubits.
	remapped = operation.with_qubits(q1, q0)
	assert remapped.qubits == (q1, q0)
	assert "X(0)*Y(1)" in str(operation)
	assert "PauliStringGateOperation" in repr(operation)
	assert operation == vqubit.PauliStringGateOperation(string, [q0, q1])


def test_single_qubit_pauli_string_operation_uses_the_pauli_gate() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	for pauli, gate in ((P.X, vqubit.X), (P.Y, vqubit.Y), (P.Z, vqubit.Z)):
		operation = vqubit.SingleQubitPauliStringGateOperation(vqubit.PauliString({q0: pauli}), q0)
		assert operation.gate == gate
		assert_matrix_close(vqubit.unitary(operation), vqubit.unitary(gate))
		# The square of any Pauli is the identity on the same qubit.
		assert_matrix_close(vqubit.unitary(operation**2), np.eye(2))
		assert (operation**2).qubits == (q0,)
		assert operation.with_qubits(q1).qubits == (q1,)


def test_pauli_string_phasor_is_a_rotation_about_the_string() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	string = vqubit.PauliString({q0: P.X, q1: P.Y})
	phasor = vqubit.PauliStringPhasor(string, half_turns=0.5)
	assert phasor.half_turns == 0.5
	assert phasor.pauli_string == string
	assert phasor.qubits == (q0, q1)
	assert isinstance(phasor.gate, vqubit.PauliStringPhasorGate)
	unitary = np.asarray(vqubit.unitary(phasor))
	# exp(iπ P / 2) squares to exp(iπ P) = -I, because P² = I.
	assert_matrix_close(unitary @ unitary, -np.eye(4))
	# Squaring the operation doubles the angle.
	assert (phasor**2).half_turns == 1.0
	assert_matrix_close(vqubit.unitary(phasor**0), np.eye(4))
	# Remapping moves the Paulis onto the new qubits, which stay sorted.
	remapped = phasor.with_qubits(q1, q0)
	assert remapped.qubits == (q0, q1)
	assert remapped.pauli_string[q0] is P.Y
	assert remapped.pauli_string[q1] is P.X
	assert "half_turns=0.5" in str(phasor)
	assert "PauliStringPhasor" in repr(phasor)


def test_a_phasor_of_a_scaled_string_has_no_unitary() -> None:
	# The coefficient would break unitarity, so the gate refuses to build one.
	q0 = vqubit.LineQubit(0)
	phasor = vqubit.PauliStringPhasor(vqubit.PauliString({q0: P.X}, 2), half_turns=0.5)
	assert phasor.gate._unitary_() is None
	with pytest.raises(TypeError, match="has no unitary"):
		vqubit.unitary(phasor)


def test_pauli_measurement_gate() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gate = vqubit.PauliMeasurementGate(vqubit.PauliString({q0: P.X, q1: P.Z}), "m")
	assert gate.num_qubits == 2
	assert vqubit.is_measurement(gate)
	assert gate.pauli_string == vqubit.PauliString({q0: P.X, q1: P.Z})
	# Without a key the string itself names the measurement.
	anonymous = vqubit.PauliMeasurementGate(vqubit.PauliString({q0: P.X}))
	assert vqubit.measurement_key_name(anonymous) == "X(0)"
	assert "X(0)" in str(anonymous)
	assert "key='m'" in repr(gate)
	with pytest.raises(TypeError, match="has no unitary"):
		vqubit.unitary(gate)


def test_measuring_pauli_strings_expands_a_sum() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	total = vqubit.PauliString({q0: P.X}) + vqubit.PauliString({q1: P.Z})
	operations = vqubit.measure_paulistring_terms(total, key="t")
	assert [vqubit.measurement_key_name(operation) for operation in operations] == ["t0", "t1"]
	single = vqubit.measure_single_paulistring(vqubit.PauliString({q0: P.X}), key="s")
	assert single.qubits == (q0,)
	assert vqubit.measurement_key_name(single) == "s"


def test_pauli_sum_exponential() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	string = vqubit.PauliString({q0: P.Z})
	exponential = vqubit.PauliSumExponential(string * 0.5, exponent=0.5)
	assert exponential.pauli_sum == vqubit.PauliSum([(string * 0.5, 1.0)])
	assert exponential.exponent == 0.5
	# exp(iπ/2 · 0.5 Z) = diag(e^{iπ/4}, e^{-iπ/4}).
	assert_matrix_close(
		exponential.matrix(),
		[[cmath.exp(1j * math.pi / 4), 0], [0, cmath.exp(-1j * math.pi / 4)]],
	)
	# With exponent zero it is the identity, and it is unitary in general.
	assert_matrix_close(vqubit.PauliSumExponential(string * 0.5, 0).matrix(), np.eye(2))
	unitary = np.asarray(vqubit.PauliSumExponential(string, 0.25).matrix())
	assert_matrix_close(unitary.conj().T @ unitary, np.eye(2))
	assert "exponent=0.5" in str(exponential)
	assert "PauliSumExponential" in repr(exponential)
	assert exponential == vqubit.PauliSumExponential(string * 0.5, exponent=0.5)
	assert exponential != vqubit.PauliSumExponential(string * 0.5, exponent=1.0)
	del q1
