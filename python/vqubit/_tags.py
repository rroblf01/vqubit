"""Tags attached to operations.

A tag is a marker value (with no effect on the operation's unitary) that a
transformer can look for, e.g. to mark a ``SWAP`` that only exists to move
qubits around during routing.
"""

from __future__ import annotations

from typing import Any

from ._value import value_equality


@value_equality
class VirtualTag:
	"""Marks an operation as virtual: it need not be implemented by hardware."""

	def _value_equality_values_(self) -> Any:
		return ()

	def __str__(self) -> str:
		return "virtual"

	def __repr__(self) -> str:
		return "vqubit.VirtualTag()"


@value_equality
class RoutingSwapTag:
	"""Marks a ``SWAP`` as a routing swap introduced by a mapping pass."""

	def _value_equality_values_(self) -> Any:
		return ()

	def __str__(self) -> str:
		return "routing-swap"

	def __repr__(self) -> str:
		return "vqubit.RoutingSwapTag()"


__all__ = ["RoutingSwapTag", "VirtualTag"]
