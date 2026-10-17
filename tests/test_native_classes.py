"""Tests for the native `@[vc_buffer]` state classes."""

from __future__ import annotations

import vqubit_native as native

X_GATE = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
H_GATE = [2**-0.5, 0.0, 2**-0.5, 0.0, 2**-0.5, 0.0, -(2**-0.5), 0.0]
ZERO_ONE = [0.0, 0.0, 1.0, 0.0]  # |1> of one qubit
ZERO_ZERO = [1.0, 0.0, 0.0, 0.0]  # |0> of one qubit


def test_native_state_apply_and_read() -> None:
	state = native.NativeState(1)
	state.set_amplitudes(ZERO_ZERO)
	state.apply_gate(X_GATE, [0])
	assert list(state.amplitudes()) == ZERO_ONE
	assert state.num_qubits == 1


def test_native_state_probabilities_and_measure() -> None:
	state = native.NativeState(1)
	state.set_amplitudes(ZERO_ONE)
	assert list(state.probabilities([0])) == [0.0, 1.0]
	# A basis state measures deterministically, whatever the draw.
	assert state.measure([0], 0.1) == 1
	assert state.measure([0], 0.9) == 1


def test_native_state_superposition() -> None:
	state = native.NativeState(1)
	state.set_amplitudes(ZERO_ZERO)
	state.apply_gate(H_GATE, [0])
	probabilities = list(state.probabilities([0]))
	assert abs(probabilities[0] - 0.5) < 1e-9
	assert abs(probabilities[1] - 0.5) < 1e-9


def test_native_state_memoryview_is_zero_copy() -> None:
	state = native.NativeState(1)
	state.set_amplitudes(ZERO_ZERO)
	view = memoryview(state)  # ty: ignore[invalid-argument-type]
	assert view.format == "d"
	view[0] = 0.25  # ty: ignore[invalid-assignment]
	assert state.amplitudes()[0] == 0.25


def test_native_density_channel() -> None:
	density = native.NativeDensity(1)
	density.set_from_amplitudes(ZERO_ZERO)
	density.apply_channel(X_GATE, 1, [0])
	matrix = density.density_matrix()
	assert matrix[2 * (1 * 2 + 1)] == 1.0  # |1><1| population
	assert density.num_qubits == 1


def test_native_density_measure() -> None:
	density = native.NativeDensity(1)
	density.set_from_amplitudes(ZERO_ONE)
	assert density.measure([0], 0.5) == 1


def test_native_tableau_bell_state() -> None:
	tableau = native.NativeTableau(2)
	tableau.apply(native.CliffordOp.h, 0, 0)  # H on qubit 0
	tableau.apply(native.CliffordOp.cnot, 0, 1)  # CNOT(0, 1)
	vector = tableau.state_vector()
	assert abs(vector[0] - 2**-0.5) < 1e-9
	assert abs(vector[6] - 2**-0.5) < 1e-9
	assert tableau.num_qubits == 2


def test_native_tableau_measurement_is_correlated() -> None:
	tableau = native.NativeTableau(2)
	tableau.apply(native.CliffordOp.h, 0, 0)
	tableau.apply(native.CliffordOp.cnot, 0, 1)
	first = tableau.measure(0, 0.9)
	second = tableau.measure(1, 0.1)
	assert first == second


def test_native_tableau_reset() -> None:
	tableau = native.NativeTableau(2)
	tableau.apply(native.CliffordOp.h, 0, 0)
	tableau.reset()
	vector = tableau.state_vector()
	assert abs(vector[0]) == 1.0
