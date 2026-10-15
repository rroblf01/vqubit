"""Tests for the noise channels."""

from __future__ import annotations

import pytest

import vqubit
from vqubit._linalg import dagger, matmul

IDENTITY = [[1 + 0j, 0j], [0j, 1 + 0j]]


def _trace_preserving(operators: list[list[list[complex]]]) -> bool:
	total = [[0j, 0j], [0j, 0j]]
	for operator in operators:
		product = matmul(dagger(operator), operator)
		for i in range(2):
			for j in range(2):
				total[i][j] += product[i][j]
	return all(abs(total[i][j] - IDENTITY[i][j]) < 1e-9 for i in range(2) for j in range(2))


def test_channel_kraus_counts() -> None:
	assert len(vqubit.kraus(vqubit.bit_flip(0.1))) == 2
	assert len(vqubit.kraus(vqubit.phase_flip(0.1))) == 2
	assert len(vqubit.kraus(vqubit.bit_phase_flip(0.1))) == 2
	assert len(vqubit.kraus(vqubit.depolarize(0.1))) == 4
	assert len(vqubit.kraus(vqubit.amplitude_damp(0.1))) == 2
	assert len(vqubit.kraus(vqubit.phase_damp(0.1))) == 2


def test_channels_are_trace_preserving() -> None:
	for channel in [
		vqubit.bit_flip(0.2),
		vqubit.phase_flip(0.3),
		vqubit.bit_phase_flip(0.4),
		vqubit.depolarize(0.15),
		vqubit.amplitude_damp(0.25),
		vqubit.phase_damp(0.35),
	]:
		assert _trace_preserving(vqubit.kraus(channel)), str(channel)


def test_unitary_gates_have_one_kraus_operator() -> None:
	operators = vqubit.kraus(vqubit.X)
	assert len(operators) == 1
	assert vqubit.allclose_up_to_global_phase(operators[0], vqubit.unitary(vqubit.X))


def test_channel_flags() -> None:
	channel = vqubit.depolarize(0.1)
	assert not vqubit.has_unitary(channel)
	assert vqubit.has_kraus(channel)
	assert vqubit.has_kraus(vqubit.X)


def test_channel_parameterization() -> None:
	t = vqubit.Parameter("t")
	channel = vqubit.bit_flip(t)
	assert vqubit.is_parameterized(channel)
	with pytest.raises(TypeError):
		vqubit.kraus(channel)
	resolved = vqubit.resolve_parameters(channel, vqubit.ParamResolver({t: 0.1}))
	assert resolved == vqubit.bit_flip(0.1)


def test_channel_equality() -> None:
	assert vqubit.bit_flip(0.1) == vqubit.bit_flip(0.1)
	assert vqubit.bit_flip(0.1) != vqubit.bit_flip(0.2)
	assert vqubit.depolarize(0.1, 2) != vqubit.depolarize(0.1, 1)


def test_channel_is_a_gate_operation() -> None:
	q0 = vqubit.LineQubit(0)
	operation = vqubit.bit_flip(0.1)(q0)
	assert operation.qubits == (q0,)
	assert operation.gate == vqubit.bit_flip(0.1)
