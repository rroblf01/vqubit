"""Tests for the value_equality helper."""

from __future__ import annotations

from typing import Any

import vqubit


@vqubit.value_equality
class _Point:
	def __init__(self, x: int, y: int) -> None:
		self.x = x
		self.y = y

	def _value_equality_values_(self) -> Any:
		return (self.x, self.y)


@vqubit.value_equality
class _Other:
	def __init__(self, x: int) -> None:
		self.x = x

	def _value_equality_values_(self) -> Any:
		return (self.x,)


def test_equal_values() -> None:
	assert _Point(1, 2) == _Point(1, 2)
	assert _Point(1, 2) != _Point(1, 3)


def test_not_equal_operator() -> None:
	assert _Point(1, 2) != _Point(1, 3)


def test_different_types_are_not_equal() -> None:
	assert _Point(1, 2) != _Other(1)


def test_hash_from_values() -> None:
	assert hash(_Point(1, 2)) == hash(_Point(1, 2))
	assert len({_Point(1, 2), _Point(1, 2), _Point(3, 4)}) == 2


def test_identity_is_equal() -> None:
	point = _Point(1, 2)
	assert point == point
