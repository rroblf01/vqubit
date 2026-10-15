"""Tests for the fSim, sqrt-iSWAP and MS gates, and the Xmon device."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close, up_to_phase
from vqubit._linalg import matmul


def test_fsim_gate() -> None:
	gate = vqubit.FSimGate(1.0, 0.5)
	assert gate.num_qubits == 2
	assert gate.theta == 1.0
	assert gate.phi == 0.5
	assert_matrix_close(
		matmul(vqubit.unitary(gate), vqubit.unitary(vqubit.inverse(gate))),
		[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
	)


def test_sycamore_is_fsim_pi_over_2_pi_over_6() -> None:
	assert up_to_phase(vqubit.unitary(vqubit.SycamoreGate()), vqubit.unitary(vqubit.FSimGate(math.pi / 2, math.pi / 6)))


def test_sqrt_iswap_squares_to_iswap() -> None:
	squared = matmul(vqubit.unitary(vqubit.SQRT_ISWAP), vqubit.unitary(vqubit.SQRT_ISWAP))
	assert_matrix_close(squared, vqubit.unitary(vqubit.ISWAP))


def test_ms_gate() -> None:
	gate = vqubit.MSGate(2, 0.5)
	assert gate.num_qubits == 2
	assert len(vqubit.unitary(gate)) == 4
	with pytest.raises(ValueError):
		vqubit.MSGate(3, 0.5)


def test_xmon_device_durations() -> None:
	device = vqubit.XmonDevice([vqubit.GridQubit(0, 0), vqubit.GridQubit(0, 1)], gate_durations={"HPowGate": 25.0})
	assert device.metadata.name == "XmonDevice"
	assert device.metadata.duration_of(vqubit.H(vqubit.GridQubit(0, 0))) == 25.0
	assert device.metadata.duration_of(vqubit.X(vqubit.GridQubit(0, 0))) is None
	device.validate_operation(vqubit.CNOT(vqubit.GridQubit(0, 0), vqubit.GridQubit(0, 1)))
