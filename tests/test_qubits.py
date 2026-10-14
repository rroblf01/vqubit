"""Tests for qubit identifiers."""

from __future__ import annotations

import pytest

import vqubit


def test_line_qubit_basics() -> None:
	q = vqubit.LineQubit(3)
	assert q.x == 3
	assert str(q) == "3"
	assert repr(q) == "vqubit.LineQubit(3)"


def test_line_qubit_range() -> None:
	assert vqubit.LineQubit.range(3) == (vqubit.LineQubit(0), vqubit.LineQubit(1), vqubit.LineQubit(2))


def test_line_qubit_validates_position() -> None:
	with pytest.raises(TypeError):
		vqubit.LineQubit(1.5)  # ty: ignore[invalid-argument-type]
	with pytest.raises(ValueError):
		vqubit.LineQubit(-1)
	with pytest.raises(TypeError):
		vqubit.LineQubit(True)


def test_line_qubit_equality_and_hash() -> None:
	assert vqubit.LineQubit(1) == vqubit.LineQubit(1)
	assert vqubit.LineQubit(1) != vqubit.LineQubit(2)
	assert hash(vqubit.LineQubit(1)) == hash(vqubit.LineQubit(1))
	assert len({vqubit.LineQubit(1), vqubit.LineQubit(1)}) == 1


def test_qubits_of_different_types_are_not_equal() -> None:
	assert vqubit.LineQubit(0) != vqubit.NamedQubit("0")
	assert vqubit.NamedQubit("a") != vqubit.GridQubit(0, 0)


def test_qubit_ordering() -> None:
	qubits = [vqubit.LineQubit(2), vqubit.LineQubit(0), vqubit.LineQubit(1)]
	assert sorted(qubits) == [vqubit.LineQubit(0), vqubit.LineQubit(1), vqubit.LineQubit(2)]


def test_grid_qubit_basics() -> None:
	q = vqubit.GridQubit(2, 3)
	assert (q.row, q.col) == (2, 3)
	assert str(q) == "(2, 3)"
	assert repr(q) == "vqubit.GridQubit(2, 3)"


def test_grid_qubit_neighbors() -> None:
	neighbors = vqubit.GridQubit(0, 0).neighbors()
	assert set(neighbors) == {
		vqubit.GridQubit(-1, 0),
		vqubit.GridQubit(1, 0),
		vqubit.GridQubit(0, -1),
		vqubit.GridQubit(0, 1),
	}


def test_grid_qubit_neighbors_on_grid() -> None:
	grid = [vqubit.GridQubit(0, 0), vqubit.GridQubit(1, 0), vqubit.GridQubit(0, 1)]
	assert set(vqubit.GridQubit(0, 0).neighbors(grid)) == {vqubit.GridQubit(1, 0), vqubit.GridQubit(0, 1)}


def test_grid_qubit_adjacency() -> None:
	assert vqubit.GridQubit(0, 0).is_adjacent(vqubit.GridQubit(0, 1))
	assert not vqubit.GridQubit(0, 0).is_adjacent(vqubit.GridQubit(1, 1))


def test_grid_qubit_square() -> None:
	square = vqubit.GridQubit.square(2)
	assert set(square) == {
		vqubit.GridQubit(0, 0),
		vqubit.GridQubit(0, 1),
		vqubit.GridQubit(1, 0),
		vqubit.GridQubit(1, 1),
	}


def test_named_qubit_basics() -> None:
	q = vqubit.NamedQubit("alice")
	assert q.name == "alice"
	assert str(q) == "alice"
	assert repr(q) == "vqubit.NamedQubit('alice')"
