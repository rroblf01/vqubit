"""Qubit managers.

A qubit manager hands out and recycles qubits while a circuit is being built
or routed. :class:`SimpleQubitManager` always allocates fresh qubits;
:class:`GreedyQubitManager` reuses freed ones to keep the register small.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

from ._qubits import LineQubit, Qid
from ._value import value_equality


@value_equality
class BorrowableQubit(Qid):
	"""A qubit that is borrowed from a qubit manager and must be released."""

	def __init__(self, qubit: Qid) -> None:
		self._qubit = qubit

	@property
	def qubit(self) -> Qid:
		return self._qubit

	def _comparison_key(self) -> Any:
		return (self._qubit._comparison_key(),)

	def _value_equality_values_(self) -> Any:
		return (self._qubit,)

	def __str__(self) -> str:
		return f"borrow({self._qubit})"

	def __repr__(self) -> str:
		return f"vqubit.BorrowableQubit({self._qubit!r})"


@value_equality
class CleanQubit(Qid):
	"""A qubit known to start in the ``|0>`` state."""

	def __init__(self, qubit: Qid) -> None:
		self._qubit = qubit

	@property
	def qubit(self) -> Qid:
		return self._qubit

	def _comparison_key(self) -> Any:
		return (self._qubit._comparison_key(),)

	def _value_equality_values_(self) -> Any:
		return (self._qubit,)

	def __str__(self) -> str:
		return f"clean({self._qubit})"

	def __repr__(self) -> str:
		return f"vqubit.CleanQubit({self._qubit!r})"


class QubitManager:
	"""Base class for a manager of the qubits used while building a circuit."""

	def qalloc(self) -> Qid:
		"""Allocate a fresh qubit."""
		raise NotImplementedError

	def qfree(self, qubit: Qid) -> None:
		"""Return an allocated qubit to the manager."""
		raise NotImplementedError

	def qborrow(self) -> BorrowableQubit:
		"""Borrow a qubit known to start in an unknown state."""
		raise NotImplementedError

	def qrelease(self, qubit: Qid) -> None:
		"""Release a borrowed qubit."""
		raise NotImplementedError

	def allocate(self, n: int) -> list[Qid]:
		"""Allocate ``n`` fresh qubits."""
		return [self.qalloc() for _ in range(n)]

	def free(self, qubits: Iterable[Qid]) -> None:
		"""Return several qubits to the manager."""
		for qubit in qubits:
			self.qfree(qubit)


class SimpleQubitManager(QubitManager):
	"""A qubit manager that always allocates a fresh ``LineQubit``."""

	def __init__(self, qubits: Sequence[Qid] | None = None, max_size: int | None = None) -> None:
		self._qubits = list(qubits) if qubits is not None else []
		self._max_size = max_size
		self._next = len(self._qubits)

	def qalloc(self) -> Qid:
		if self._max_size is not None and self._next >= self._max_size:
			raise ValueError("SimpleQubitManager exceeded max_size")
		qubit = LineQubit(self._next)
		self._next += 1
		return qubit

	def qfree(self, qubit: Qid) -> None:
		return None

	def qborrow(self) -> BorrowableQubit:
		return BorrowableQubit(self.qalloc())

	def qrelease(self, qubit: Qid) -> None:
		return None

	@property
	def qubits(self) -> list[Qid]:
		return list(self._qubits)


class GreedyQubitManager(QubitManager):
	"""A qubit manager that reuses freed qubits when possible."""

	def __init__(
		self, qubits: Sequence[Qid] | None = None, max_size: int | None = None, maximize_reuse: bool = False
	) -> None:
		self._free_pool = list(qubits) if qubits is not None else []
		self._in_use: set[Qid] = set()
		self._max_size = max_size
		self._next = len(self._free_pool)
		self._maximize_reuse = maximize_reuse

	def qalloc(self) -> Qid:
		if self._free_pool:
			qubit = self._free_pool.pop()
		else:
			if self._max_size is not None and self._next >= self._max_size:
				raise ValueError("GreedyQubitManager exceeded max_size")
			qubit = LineQubit(self._next)
			self._next += 1
		self._in_use.add(qubit)
		return qubit

	def qfree(self, qubit: Qid) -> None:
		if qubit in self._in_use:
			self._in_use.discard(qubit)
			self._free_pool.append(qubit)

	def qborrow(self) -> BorrowableQubit:
		return BorrowableQubit(self.qalloc())

	def qrelease(self, qubit: Qid) -> None:
		self.qfree(qubit)


__all__ = ["BorrowableQubit", "CleanQubit", "GreedyQubitManager", "QubitManager", "SimpleQubitManager"]
