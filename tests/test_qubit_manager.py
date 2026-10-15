"""Tests for the qubit managers."""

from __future__ import annotations

import pytest

import vqubit


def test_simple_qubit_manager() -> None:
	manager = vqubit.SimpleQubitManager()
	first = manager.qalloc()
	second = manager.qalloc()
	assert first == vqubit.LineQubit(0)
	assert second == vqubit.LineQubit(1)
	assert manager.allocate(2) == [vqubit.LineQubit(2), vqubit.LineQubit(3)]


def test_simple_qubit_manager_max_size() -> None:
	manager = vqubit.SimpleQubitManager(max_size=1)
	manager.qalloc()
	with pytest.raises(ValueError):
		manager.qalloc()


def test_greedy_qubit_manager_reuses() -> None:
	manager = vqubit.GreedyQubitManager()
	qubit = manager.qalloc()
	manager.qfree(qubit)
	assert manager.qalloc() == qubit


def test_borrow_and_clean_qubits() -> None:
	borrowed = vqubit.BorrowableQubit(vqubit.LineQubit(0))
	clean = vqubit.CleanQubit(vqubit.LineQubit(1))
	assert borrowed.qubit == vqubit.LineQubit(0)
	assert clean.qubit == vqubit.LineQubit(1)
	assert isinstance(borrowed, vqubit.Qid)
	assert isinstance(clean, vqubit.Qid)
	assert vqubit.BorrowableQubit(vqubit.LineQubit(0)) == vqubit.BorrowableQubit(vqubit.LineQubit(0))
