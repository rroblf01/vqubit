"""Qubit ordering.

A :class:`QubitOrder` turns the set of qubits a circuit touches into an ordered
tuple, which is what a simulator needs. ``QubitOrder.DEFAULT`` sorts them;
``QubitOrder.explicit`` fixes an order; ``QubitOrder.sorted_by`` sorts by a key.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from typing import Any

from ._qubits import Qid


class QubitOrder:
	"""A callable that orders a collection of qubits."""

	DEFAULT: QubitOrder

	def __init__(self, order_func: Callable[[tuple[Qid, ...]], Iterable[Qid]]) -> None:
		self._order_func = order_func

	def order_for(self, qubits: Iterable[Qid]) -> tuple[Qid, ...]:
		"""Return ``qubits`` in this order."""
		return tuple(self._order_func(tuple(qubits)))

	def __call__(self, qubits: Iterable[Qid]) -> tuple[Qid, ...]:
		return self.order_for(qubits)

	def __repr__(self) -> str:
		return "vqubit.QubitOrder(...)"

	@staticmethod
	def explicit(qubits: Iterable[Qid]) -> QubitOrder:
		"""Order qubits as listed, appending any the list does not mention."""
		explicit = tuple(qubits)

		def order(qubits: tuple[Qid, ...]) -> list[Qid]:
			present = set(qubits)
			ordered = [q for q in explicit if q in present]
			ordered.extend(q for q in qubits if q not in set(explicit))
			return ordered

		return QubitOrder(order)

	@staticmethod
	def sorted_by(key: Callable[[Qid], Any]) -> QubitOrder:
		"""Order qubits by ``key``."""

		def order(qubits: tuple[Qid, ...]) -> list[Qid]:
			return sorted(qubits, key=key)

		return QubitOrder(order)


QubitOrder.DEFAULT = QubitOrder.sorted_by(lambda qubit: qubit)

QubitOrderOrList = QubitOrder | Sequence[Qid]

__all__ = ["QubitOrder", "QubitOrderOrList"]
