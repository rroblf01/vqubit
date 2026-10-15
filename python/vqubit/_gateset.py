"""Gatesets and gate families.

A :class:`GateFamily` recognizes operations built from a particular gate (or
satisfying a predicate), and a :class:`Gateset` is a named collection of gate
families. They are what target gatesets use to decide whether an operation is
already in the set of gates a device can run.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any

from ._gates import Gate
from ._ops import Operation


class GateFamily:
	"""A family of operations, matched by gate equality or a predicate."""

	def __init__(
		self,
		gate: Gate | None = None,
		*,
		name: str | None = None,
		description: str | None = None,
		predicate: Callable[[Gate], bool] | None = None,
	) -> None:
		if gate is None and predicate is None:
			raise ValueError("GateFamily needs a gate or a predicate")
		self.gate = gate
		self.name = name if name is not None else (str(gate) if gate is not None else "GateFamily")
		self.description = description
		self.predicate = predicate

	def _matches_gate(self, gate: Gate) -> bool:
		if self.gate is not None and gate == self.gate:
			return True
		return bool(self.predicate is not None and self.predicate(gate))

	def __contains__(self, item: Any) -> bool:
		gate = item.gate if isinstance(item, Operation) else item
		if not isinstance(gate, Gate):
			return False
		return self._matches_gate(gate)

	def validate(self, operation: Operation) -> None:
		"""Raise ``ValueError`` if ``operation`` is not in this family."""
		if operation not in self:
			raise ValueError(f"{operation!r} is not in {self.name}")

	def __repr__(self) -> str:
		return f"vqubit.GateFamily({self.gate!r}, name={self.name!r})"


class Gateset:
	"""A named collection of gates or gate families."""

	def __init__(
		self,
		gates: Iterable[Gate | GateFamily],
		*,
		name: str | None = None,
		description: str | None = None,
	) -> None:
		self.gates: tuple[Any, ...] = tuple(gates)
		self.name = name if name is not None else "Gateset"
		self.description = description

	def _families(self) -> list[GateFamily]:
		return [entry if isinstance(entry, GateFamily) else GateFamily(entry) for entry in self.gates]

	def __contains__(self, item: Any) -> bool:
		return any(item in family for family in self._families())

	def validate(self, operation: Operation) -> None:
		"""Raise ``ValueError`` if ``operation`` is not in this gateset."""
		if operation not in self:
			raise ValueError(f"{operation!r} is not in {self.name}")

	def with_added_gates(self, *gates: Gate | GateFamily, name: str | None = None) -> Gateset:
		"""Return a new gateset with the given gates added."""
		return Gateset([*self.gates, *gates], name=name or self.name, description=self.description)

	def __repr__(self) -> str:
		return f"vqubit.Gateset({list(self.gates)!r}, name={self.name!r})"


__all__ = ["GateFamily", "Gateset"]
