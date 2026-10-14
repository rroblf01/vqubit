"""Tests for the decompose protocol."""

from __future__ import annotations

import pytest

import vqubit
from tests.helpers import assert_matrix_close


def _unitary_of(operations: list[vqubit.Operation]) -> list[list[complex]]:
	return vqubit.unitary(vqubit.Circuit(list(operations)))


def test_ccx_decomposition_matches() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	operations = vqubit.decompose(vqubit.CCX(q0, q1, q2))
	assert len(operations) > 1
	assert all(not isinstance(op.gate, vqubit.CCXGate) for op in operations)
	assert_matrix_close(_unitary_of(operations), vqubit.unitary(vqubit.CCX))


def test_cswap_decomposition_matches() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	operations = vqubit.decompose(vqubit.CSWAP(q0, q1, q2))
	assert all(not isinstance(op.gate, vqubit.CSWAPGate) for op in operations)
	assert_matrix_close(_unitary_of(operations), vqubit.unitary(vqubit.CSWAP))


def test_decompose_once_is_a_single_level() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	level = vqubit.decompose_once(vqubit.CSWAP(q0, q1, q2))
	# CSWAP decomposes into CNOT, CCX, CNOT; the CCX survives one level.
	assert any(isinstance(op.gate, vqubit.CCXGate) for op in level)


def test_decompose_of_an_atomic_gate_returns_it() -> None:
	q0 = vqubit.LineQubit(0)
	assert vqubit.decompose(vqubit.H(q0)) == [vqubit.H(q0)]
	assert vqubit.decompose(vqubit.H) == [vqubit.H]


def test_decompose_of_a_circuit() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.CCX(q0, q1, q2), vqubit.H(q0))
	operations = vqubit.decompose(circuit)
	assert all(not isinstance(op.gate, vqubit.CCXGate) for op in operations)
	assert any(op == vqubit.H(q0) for op in operations)


def test_decompose_rejects_unknown() -> None:
	with pytest.raises(TypeError):
		vqubit.decompose(5)
