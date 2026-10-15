"""Tests for the extra noise channels."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close


def test_asymmetric_depolarize() -> None:
	channel = vqubit.asymmetric_depolarize(0.1, 0.2, 0.3)
	assert channel.num_qubits == 1
	assert channel.p_x == pytest.approx(0.1)
	assert channel.p_y == pytest.approx(0.2)
	assert channel.p_z == pytest.approx(0.3)
	assert sum(channel.error_probabilities.values()) == pytest.approx(1)
	assert len(vqubit.kraus(channel)) == 4
	assert vqubit.is_cptp(kraus_ops=vqubit.kraus(channel))


def test_asymmetric_depolarize_two_qubits() -> None:
	channel = vqubit.asymmetric_depolarize(error_probabilities={"XX": 0.2})
	assert channel.num_qubits == 2
	assert sum(channel.error_probabilities.values()) == pytest.approx(1)
	with pytest.raises(ValueError):
		vqubit.asymmetric_depolarize(error_probabilities={"AB": 0.2})


def test_generalized_amplitude_damp() -> None:
	channel = vqubit.generalized_amplitude_damp(0.3, 0.4)
	assert channel.num_qubits == 1
	assert channel.p == pytest.approx(0.3)
	assert channel.gamma == pytest.approx(0.4)
	assert len(vqubit.kraus(channel)) == 4
	assert vqubit.is_cptp(kraus_ops=vqubit.kraus(channel))
	with pytest.raises(ValueError):
		vqubit.generalized_amplitude_damp(1.5, 0.4)


def test_kraus_channel() -> None:
	channel = vqubit.KrausChannel([[[1, 0], [0, 0]], [[0, 1], [0, 0]]])
	assert channel.num_qubits == 1
	assert len(vqubit.kraus(channel)) == 2
	assert vqubit.is_cptp(kraus_ops=vqubit.kraus(channel))
	copied = vqubit.KrausChannel.from_channel(vqubit.bit_flip(0.2), key="k")
	assert copied.num_qubits == 1
	with pytest.raises(ValueError):
		vqubit.KrausChannel([])


def test_mixed_unitary_channel() -> None:
	channel = vqubit.MixedUnitaryChannel([(0.5, [[1, 0], [0, 1]]), (0.5, [[0, 1], [1, 0]])])
	assert channel.num_qubits == 1
	assert len(channel._mixture_()) == 2
	assert vqubit.is_cptp(kraus_ops=vqubit.kraus(channel))
	with pytest.raises(ValueError):
		vqubit.MixedUnitaryChannel([(0.4, [[1, 0], [0, 1]])])


def test_random_gate_channel() -> None:
	channel = vqubit.RandomGateChannel(sub_gate=vqubit.X, probability=0.2)
	assert channel.num_qubits == 1
	assert channel.sub_gate == vqubit.X
	assert channel.probability == pytest.approx(0.2)
	assert len(vqubit.kraus(channel)) == 2
	with pytest.raises(ValueError):
		vqubit.RandomGateChannel(sub_gate=vqubit.X, probability=1.5)


def test_state_preparation_channel() -> None:
	channel = vqubit.StatePreparationChannel([1, 0])
	assert channel.num_qubits == 1
	operators = vqubit.kraus(channel)
	assert len(operators) == 2
	# K_0 = |0><0|, K_1 = |0><1|.
	assert_matrix_close(operators[0], [[1, 0], [0, 0]])
	assert_matrix_close(operators[1], [[0, 1], [0, 0]])
	assert vqubit.is_cptp(kraus_ops=operators)


def test_state_preparation_channel_normalizes() -> None:
	channel = vqubit.StatePreparationChannel([1, 1])
	norm = math.sqrt(sum(abs(z) ** 2 for z in channel.state))
	assert norm == pytest.approx(1)
