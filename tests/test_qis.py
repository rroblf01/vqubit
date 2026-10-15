"""Tests for the quantum-state helpers (qis)."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import up_to_phase

SQRT_HALF = math.sqrt(0.5)


def test_to_valid_state_vector_basis() -> None:
	assert vqubit.to_valid_state_vector(1, (2,)) == [0j, 1 + 0j]
	assert vqubit.to_valid_state_vector(3, (2, 2)) == [0j, 0j, 0j, 1 + 0j]


def test_to_valid_state_vector_rejects_unnormalized() -> None:
	with pytest.raises(ValueError):
		vqubit.to_valid_state_vector([1, 1], (2,))
	with pytest.raises(ValueError):
		vqubit.to_valid_state_vector(4, (2, 2))


def test_state_vector_from_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	vector = vqubit.state_vector(vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1)))
	assert vector[0] == pytest.approx(SQRT_HALF)
	assert vector[3] == pytest.approx(SQRT_HALF)


def test_bloch_vector() -> None:
	assert vqubit.bloch_vector_from_state_vector([SQRT_HALF, SQRT_HALF], 0) == pytest.approx((1.0, 0.0, 0.0))
	assert vqubit.bloch_vector_from_state_vector([0j, 1 + 0j], 0) == pytest.approx((0.0, 0.0, -1.0))


def test_density_matrix_from_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rho = vqubit.density_matrix(vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1)))
	assert rho[0][0] == pytest.approx(0.5)
	assert rho[3][3] == pytest.approx(0.5)
	assert sum(rho[i][i].real for i in range(4)) == pytest.approx(1.0)


def test_partial_trace_of_bell_is_maximally_mixed() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	rho = vqubit.density_matrix(vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1)))
	reduced = vqubit.partial_trace(rho, [0])
	assert reduced[0][0] == pytest.approx(0.5)
	assert reduced[1][1] == pytest.approx(0.5)
	assert reduced[0][1] == pytest.approx(0.0)


def test_fidelity() -> None:
	assert vqubit.fidelity([1, 0], [1, 0]) == pytest.approx(1.0)
	assert vqubit.fidelity([1, 0], [0, 1]) == pytest.approx(0.0)
	q0, q1 = vqubit.LineQubit.range(2)
	rho = vqubit.density_matrix(vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1)))
	assert vqubit.fidelity(rho, rho) == pytest.approx(1.0)
	assert vqubit.fidelity([SQRT_HALF, 0, 0, SQRT_HALF], rho) == pytest.approx(1.0)


def test_von_neumann_entropy() -> None:
	assert vqubit.von_neumann_entropy([[1, 0], [0, 0]]) == pytest.approx(0.0)
	assert vqubit.von_neumann_entropy([[0.5, 0], [0, 0.5]]) == pytest.approx(1.0)


def test_random_states() -> None:
	vector = vqubit.random_state_vector(3, seed=0)
	assert sum(abs(x) ** 2 for x in vector) == pytest.approx(1.0)
	rho = vqubit.random_density_matrix(2, seed=0)
	assert sum(rho[i][i].real for i in range(4)) == pytest.approx(1.0)


def test_kraus_channel_validates() -> None:
	assert len(vqubit.kraus_channel([[[1, 0], [0, 1]]])) == 1
	with pytest.raises(ValueError):
		vqubit.kraus_channel([[[2, 0], [0, 2]]])


def test_deconstruct_angles_reconstruct() -> None:
	x, y, z = vqubit.deconstruct_single_qubit_matrix_into_angles(vqubit.unitary(vqubit.H))
	gate = vqubit.PhasedXZGate(x_exponent=x, z_exponent=z, axis_phase_exponent=y)
	assert up_to_phase(vqubit.unitary(gate), vqubit.unitary(vqubit.H))
