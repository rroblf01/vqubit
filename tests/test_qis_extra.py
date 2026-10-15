"""Tests for the extra qis helpers."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close, up_to_phase


def test_validate_qid_shape() -> None:
	vqubit.validate_qid_shape([2, 3], 2)
	with pytest.raises(ValueError):
		vqubit.validate_qid_shape([2, 3], 3)
	with pytest.raises(ValueError):
		vqubit.validate_qid_shape([0, 2])


def test_validate_state_helpers() -> None:
	vqubit.validate_normalized_state_vector([1, 0], [2])
	with pytest.raises(ValueError):
		vqubit.validate_normalized_state_vector([1, 1], [2])
	vqubit.validate_density_matrix([[1, 0], [0, 0]], [2])
	with pytest.raises(ValueError):
		vqubit.validate_density_matrix([[1, 0], [0, 1]], [2])
	vqubit.validate_indices([0, 1], [2, 2])
	with pytest.raises(ValueError):
		vqubit.validate_indices([2, 0], [2, 2])


def test_density_matrix_from_state_vector() -> None:
	bell = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	rho = vqubit.density_matrix_from_state_vector(bell)
	assert rho[0][0] == pytest.approx(0.5)
	assert rho[0][3] == pytest.approx(0.5)
	reduced = vqubit.density_matrix_from_state_vector(bell, [0])
	assert reduced[0][0] == pytest.approx(0.5)


def test_dirac_notation() -> None:
	assert vqubit.dirac_notation([1, 0]) == "1.00|0>"
	assert vqubit.dirac_notation([0, 0, 0, 1]) == "1.00|11>"
	assert vqubit.dirac_notation([2**-0.5, 0, 0, 2**-0.5]) == "0.71|00> + 0.71|11>"


def test_one_hot_and_eye_tensor() -> None:
	assert vqubit.one_hot([1, 0], [2, 2]) == [[0j, 0j], [1 + 0j, 0j]]
	assert vqubit.eye_tensor([2, 2]) == [[1 + 0j, 0j], [0j, 1 + 0j]]
	assert vqubit.eye_tensor([2, 2, 2, 2])[0][0][0][0] == 1
	assert vqubit.eye_tensor([2, 2, 2, 2])[0][1][0][1] == 0


def test_entanglement_fidelity() -> None:
	assert vqubit.entanglement_fidelity([1, 0], [1, 0]) == pytest.approx(1)
	assert vqubit.entanglement_fidelity([1, 0], [0, 1]) == pytest.approx(0)


def test_choi_kraus_superoperator_round_trips() -> None:
	unitary = vqubit.unitary(vqubit.X)
	choi = vqubit.kraus_to_choi([unitary])
	recovered = vqubit.choi_to_kraus(choi)
	assert len(recovered) == 1
	assert up_to_phase(recovered[0], unitary)
	superoperator = vqubit.kraus_to_superoperator([unitary])
	assert vqubit.choi_to_superoperator(choi) == superoperator
	assert vqubit.superoperator_to_choi(superoperator) == choi


def test_channel_round_trip() -> None:
	kraus_ops = vqubit.kraus(vqubit.bit_flip(0.3))
	superoperator = vqubit.kraus_to_superoperator(kraus_ops)
	recovered = vqubit.superoperator_to_kraus(superoperator)
	assert len(recovered) == 2
	assert_matrix_close(vqubit.kraus_to_superoperator(recovered), superoperator)


def test_operation_to_choi_and_superoperator() -> None:
	q0 = vqubit.LineQubit(0)
	choi = vqubit.operation_to_choi(vqubit.X(q0))
	assert len(choi) == 4
	superoperator = vqubit.operation_to_superoperator(vqubit.X(q0))
	assert len(superoperator) == 4


def test_quantum_state() -> None:
	state = vqubit.quantum_state([1, 0])
	assert state.state_vector() == [1 + 0j, 0j]
	assert_matrix_close(state.density_matrix(), [[1, 0], [0, 0]])
	density = vqubit.quantum_state([[1, 0], [0, 0]])
	assert density.density_matrix() == [[1 + 0j, 0j], [0j, 0j]]


def test_stabilizer_state() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	final_state = vqubit.CliffordSimulator().simulate(circuit).final_state
	state = vqubit.StabilizerState(final_state)
	assert len(state.state_vector()) == 4
	assert_matrix_close(state.density_matrix(), vqubit.density_matrix_from_state_vector(state.state_vector()))
