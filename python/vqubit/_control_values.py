"""Boolean control-value specifications.

A control-value specification describes *when* a controlled gate fires. The
simple case is one condition per control qubit (all ``True``); these classes
cover the more general cases of a sum of products (``OR`` of ``AND`` clauses)
and a product of sums (``AND`` of ``OR`` clauses).
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from ._qubits import Qid
from ._value import value_equality


@value_equality
class AbstractControlValues:
	"""Base class for a boolean formula over the values of some qubits.

	``terms`` is a sequence of clauses; each clause is a sequence of
	``(qubit, value)`` pairs.
	"""

	def __init__(self, qubits: Sequence[Qid], terms: Sequence[Sequence[tuple[Qid, bool]]]) -> None:
		self._qubits = tuple(qubits)
		self._terms = tuple(tuple((q, bool(v)) for q, v in clause) for clause in terms)

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits

	@property
	def terms(self) -> tuple[tuple[tuple[Qid, bool], ...], ...]:
		return self._terms

	def _condition_index(self, conditions: Sequence[bool]) -> dict[Qid, bool]:
		if len(conditions) != len(self._qubits):
			raise ValueError("conditions length must match the number of control qubits")
		return dict(zip(self._qubits, conditions, strict=True))

	def __contains__(self, conditions: Sequence[bool]) -> bool:
		raise NotImplementedError

	def validate(self, conditions: Sequence[bool]) -> None:
		"""Raise ``ValueError`` if ``conditions`` do not satisfy the formula."""
		if conditions not in self:
			raise ValueError(f"{conditions!r} is not a valid control assignment for {self.to_str()}")

	def to_str(self) -> str:
		raise NotImplementedError

	def _value_equality_values_(self) -> Any:
		return (self._qubits, self._terms)

	def __str__(self) -> str:
		return self.to_str()

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}({self._qubits!r}, {self._terms!r})"


class SumOfProducts(AbstractControlValues):
	"""A control formula that is an ``OR`` of ``AND`` clauses."""

	def __contains__(self, conditions: Sequence[bool]) -> bool:
		assignment = self._condition_index(conditions)
		return any(all(assignment[q] == v for q, v in clause) for clause in self._terms)

	def to_str(self) -> str:
		clauses = ["&".join(("" if v else "!") + str(q) for q, v in clause) for clause in self._terms]
		return "|".join(f"({c})" for c in clauses) or "0"


class ProductOfSums(AbstractControlValues):
	"""A control formula that is an ``AND`` of ``OR`` clauses."""

	def __contains__(self, conditions: Sequence[bool]) -> bool:
		assignment = self._condition_index(conditions)
		return all(any(assignment[q] == v for q, v in clause) for clause in self._terms)

	def to_str(self) -> str:
		clauses = ["|".join(("" if v else "!") + str(q) for q, v in clause) for clause in self._terms]
		return "&".join(f"({c})" for c in clauses) or "1"


__all__ = ["AbstractControlValues", "ProductOfSums", "SumOfProducts"]
