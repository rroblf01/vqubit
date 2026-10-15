"""Measurement keys.

A :class:`MeasurementKey` names a classical result. It is what a measurement
gate and a classical control refer to, and it carries an optional path so keys
can be nested (``outer.inner``).
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from ._value import value_equality


@value_equality
class MeasurementKey:
	"""A name for a classical measurement result."""

	def __init__(self, name: str, path: Sequence[str] = ()) -> None:
		if not isinstance(name, str):
			raise TypeError("MeasurementKey name must be a str")
		self._name = name
		self._path = tuple(str(part) for part in path)

	@property
	def name(self) -> str:
		return self._name

	@property
	def path(self) -> tuple[str, ...]:
		return self._path

	def _value_equality_values_(self) -> Any:
		return (self._name, self._path)

	def __str__(self) -> str:
		return ".".join((*self._path, self._name))

	def __repr__(self) -> str:
		return f"vqubit.MeasurementKey({self!s})"


__all__ = ["MeasurementKey"]
