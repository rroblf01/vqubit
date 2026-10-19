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
from ._parameter import is_symbolic


def _gate_name(gate: Gate | type[Gate]) -> str:
	"""Name a gate type by its class and a gate instance by its ``str``."""
	return gate.__name__ if isinstance(gate, type) else str(gate)


class GateFamily:
	"""A family of operations, matched by gate equality, type, or a predicate.

	``gate`` is an instance (matched with ``==``) or a class (matched with
	``isinstance``), as in Cirq: ``GateFamily(vq.CZ)`` means *that* gate while
	``ParallelGateFamily(vq.ZPowGate)`` means *any* ``ZPowGate``.
	"""

	def __init__(
		self,
		gate: Gate | type[Gate] | None = None,
		*,
		name: str | None = None,
		description: str | None = None,
		predicate: Callable[[Gate], bool] | None = None,
	) -> None:
		if gate is None and predicate is None:
			raise ValueError("GateFamily needs a gate or a predicate")
		self.gate = gate
		self.name = name if name is not None else (_gate_name(gate) if gate is not None else "GateFamily")
		self.description = description
		self.predicate = predicate

	def _matches_gate(self, gate: Gate) -> bool:
		"""``True`` when ``gate`` belongs to this family.

		A family carrying a predicate is decided by it (that is how the
		subclasses below add their extra conditions); otherwise an instance is
		compared with ``==`` and a class with ``isinstance``.
		"""
		if self.predicate is not None:
			return bool(self.predicate(gate))
		if self.gate is None:
			return False
		if isinstance(self.gate, type):
			return isinstance(gate, self.gate)
		return bool(gate == self.gate)

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
		gates: Iterable[Gate | type[Gate] | GateFamily],
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

	def with_added_gates(self, *gates: Gate | type[Gate] | GateFamily, name: str | None = None) -> Gateset:
		"""Return a new gateset with the given gates added."""
		return Gateset([*self.gates, *gates], name=name or self.name, description=self.description)

	def __repr__(self) -> str:
		return f"vqubit.Gateset({list(self.gates)!r}, name={self.name!r})"


class AnyUnitaryGateFamily(GateFamily):
	"""A gate family matching every gate that has a unitary.

	``num_qubits``, when given, keeps only the gates of that width, as in Cirq.
	"""

	def __init__(self, num_qubits: int | None = None) -> None:
		if num_qubits is not None and num_qubits <= 0:
			raise ValueError("num_qubits must be a positive integer")
		self._num_qubits = num_qubits
		width = "" if num_qubits is None else f" on {num_qubits} qubits"
		super().__init__(
			predicate=self._match,
			name=f"{num_qubits if num_qubits is not None else 'Any'}-Qubit UnitaryGateFamily",
			description=f"Accepts any unitary gate{width}",
		)

	def _match(self, gate: Gate) -> bool:
		from ._protocols import has_unitary

		if self._num_qubits is not None and gate.num_qubits != self._num_qubits:
			return False
		return has_unitary(gate)


class AnyIntegerPowerGateFamily(GateFamily):
	"""A gate family matching integer powers of a base gate."""

	def __init__(
		self,
		gate: type[Gate] | Gate,
		name: str | None = None,
		description: str | None = None,
	) -> None:
		if not (isinstance(gate, Gate) or (isinstance(gate, type) and issubclass(gate, Gate))):
			raise TypeError("AnyIntegerPowerGateFamily needs a gate class or a gate instance")
		self._base_type = gate if isinstance(gate, type) else type(gate)
		super().__init__(
			gate,
			name=name or f"AnyIntegerPowerGateFamily({_gate_name(gate)})",
			description=description,
			predicate=self._match,
		)

	def _match(self, gate: Gate) -> bool:
		"""``True`` for an integer power of the base gate (a symbolic one is not)."""
		if not isinstance(gate, self._base_type):
			return False
		exponent = getattr(gate, "exponent", None)
		if exponent is None or isinstance(exponent, bool) or is_symbolic(exponent):
			return False
		return float(exponent).is_integer()


class ParallelGateFamily(GateFamily):
	"""A gate family matching a base gate applied in parallel to disjoint qubits.

	``gate`` may be a gate, a gate *type* (any instance of it is accepted, as in
	Cirq) or another family, and ``max_parallel_allowed`` caps the width of the
	gate the family accepts.
	"""

	def __init__(
		self,
		gate: Gate | type[Gate] | GateFamily,
		*,
		name: str | None = None,
		description: str | None = None,
		max_parallel_allowed: int | None = None,
	) -> None:
		from ._gates_more import ParallelGate

		if isinstance(gate, ParallelGate):
			if max_parallel_allowed is None:
				max_parallel_allowed = gate.num_qubits
			gate = gate.sub_gate
		self._inner = gate if isinstance(gate, GateFamily) else GateFamily(gate)
		self._max_parallel_allowed = max_parallel_allowed
		super().__init__(
			predicate=self._match,
			name=name or f"ParallelGateFamily({self._inner.name})",
			description=description,
		)

	def _match(self, gate: Gate) -> bool:
		from ._gates_more import ParallelGate

		if self._max_parallel_allowed is not None and gate.num_qubits > self._max_parallel_allowed:
			return False
		if self._inner._matches_gate(gate):
			return True
		return isinstance(gate, ParallelGate) and self._inner._matches_gate(gate.sub_gate)


__all__ = [
	"AnyIntegerPowerGateFamily",
	"AnyUnitaryGateFamily",
	"GateFamily",
	"Gateset",
	"ParallelGateFamily",
]
