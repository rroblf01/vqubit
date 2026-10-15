"""Qubit identifiers.

A :class:`Qid` says *where* a qubit lives. vqubit offers three concrete kinds,
mirroring Cirq:

- :class:`LineQubit` numbers qubits on a line: ``0, 1, 2, ...``.
- :class:`GridQubit` places them on a two-dimensional grid, ``(row, col)``.
- :class:`NamedQubit` gives them a name.

They are immutable values, ordered so that a mixed collection can still be
sorted deterministically.
"""

from __future__ import annotations

from typing import Any

from ._value import value_equality


@value_equality
class Qid:
	"""A quantum identifier: the abstract location of a qubit."""

	def _comparison_key(self) -> Any:
		raise NotImplementedError

	def _value_equality_values_(self) -> Any:
		return self._comparison_key()

	def _type_order(self) -> int:
		return 0

	def __lt__(self, other: object) -> bool:
		if not isinstance(other, Qid):
			return NotImplemented
		return (self._type_order(), self._comparison_key()) < (other._type_order(), other._comparison_key())

	def __le__(self, other: object) -> bool:
		if not isinstance(other, Qid):
			return NotImplemented
		return (self._type_order(), self._comparison_key()) <= (other._type_order(), other._comparison_key())

	def __gt__(self, other: object) -> bool:
		if not isinstance(other, Qid):
			return NotImplemented
		return (self._type_order(), self._comparison_key()) > (other._type_order(), other._comparison_key())

	def __ge__(self, other: object) -> bool:
		if not isinstance(other, Qid):
			return NotImplemented
		return (self._type_order(), self._comparison_key()) >= (other._type_order(), other._comparison_key())


class Qubit(Qid):
	"""The base class of the concrete qubit kinds. Not instantiated directly."""


class LineQubit(Qubit):
	"""A qubit on a line, identified by a non-negative integer position."""

	def __init__(self, x: int) -> None:
		if not isinstance(x, int) or isinstance(x, bool):
			raise TypeError(f"LineQubit position must be an int, got {type(x).__name__}")
		if x < 0:
			raise ValueError(f"LineQubit position must be non-negative, got {x}")
		self._x = x

	@property
	def x(self) -> int:
		return self._x

	def _comparison_key(self) -> int:
		return self._x

	def __str__(self) -> str:
		return str(self._x)

	def __repr__(self) -> str:
		return f"vqubit.LineQubit({self._x})"

	@classmethod
	def range(cls, n: int) -> tuple[LineQubit, ...]:
		"""Return ``n`` qubits numbered ``0`` through ``n - 1``."""
		return tuple(cls(i) for i in range(n))


class GridQubit(Qubit):
	"""A qubit on a two-dimensional grid, identified by ``(row, col)``."""

	def __init__(self, row: int, col: int) -> None:
		if not isinstance(row, int) or not isinstance(col, int) or isinstance(row, bool) or isinstance(col, bool):
			raise TypeError("GridQubit row and col must be ints")
		self._row = row
		self._col = col

	@property
	def row(self) -> int:
		return self._row

	@property
	def col(self) -> int:
		return self._col

	def _comparison_key(self) -> tuple[int, int]:
		return (self._row, self._col)

	def __str__(self) -> str:
		return f"({self._row}, {self._col})"

	def __repr__(self) -> str:
		return f"vqubit.GridQubit({self._row}, {self._col})"

	def neighbors(self, on_grid: Any = None) -> list[GridQubit]:
		"""Return the grid qubits sharing an edge with this one."""
		candidates = [
			GridQubit(self._row - 1, self._col),
			GridQubit(self._row + 1, self._col),
			GridQubit(self._row, self._col - 1),
			GridQubit(self._row, self._col + 1),
		]
		if on_grid is None:
			return candidates
		allowed = set(on_grid)
		return [q for q in candidates if q in allowed]

	def is_adjacent(self, other: GridQubit) -> bool:
		"""Report whether this qubit shares an edge with ``other``."""
		return abs(self._row - other._row) + abs(self._col - other._col) == 1

	@classmethod
	def square(cls, diameter: int, top: int = 0, left: int = 0) -> list[GridQubit]:
		"""Return a ``diameter``-by-``diameter`` square of qubits."""
		return cls.rect(rows=range(top, top + diameter), cols=range(left, left + diameter))

	@classmethod
	def rect(cls, rows: Any, cols: Any) -> list[GridQubit]:
		"""Return every qubit in the product of ``rows`` and ``cols``."""
		return [cls(row, col) for row in rows for col in cols]


class NamedQubit(Qubit):
	"""A qubit identified by a name."""

	def __init__(self, name: str) -> None:
		if not isinstance(name, str):
			raise TypeError(f"NamedQubit name must be a str, got {type(name).__name__}")
		self._name = name

	@property
	def name(self) -> str:
		return self._name

	def _comparison_key(self) -> str:
		return self._name

	def __str__(self) -> str:
		return self._name

	def __repr__(self) -> str:
		return f"vqubit.NamedQubit({self._name!r})"


class LineQid(Qid):
	"""A qid on a line, identified by a non-negative integer position."""

	def __init__(self, x: int) -> None:
		if not isinstance(x, int) or isinstance(x, bool):
			raise TypeError(f"LineQid position must be an int, got {type(x).__name__}")
		if x < 0:
			raise ValueError(f"LineQid position must be non-negative, got {x}")
		self._x = x

	@property
	def x(self) -> int:
		return self._x

	def _comparison_key(self) -> int:
		return self._x

	def __str__(self) -> str:
		return str(self._x)

	def __repr__(self) -> str:
		return f"vqubit.LineQid({self._x})"

	@classmethod
	def range(cls, n: int) -> tuple[LineQid, ...]:
		"""Return ``n`` qids numbered ``0`` through ``n - 1``."""
		return tuple(cls(i) for i in range(n))


class GridQid(Qid):
	"""A qid on a two-dimensional grid, identified by ``(row, col)``."""

	def __init__(self, row: int, col: int) -> None:
		if not isinstance(row, int) or not isinstance(col, int) or isinstance(row, bool) or isinstance(col, bool):
			raise TypeError("GridQid row and col must be ints")
		self._row = row
		self._col = col

	@property
	def row(self) -> int:
		return self._row

	@property
	def col(self) -> int:
		return self._col

	def _comparison_key(self) -> tuple[int, int]:
		return (self._row, self._col)

	def __str__(self) -> str:
		return f"({self._row}, {self._col})"

	def __repr__(self) -> str:
		return f"vqubit.GridQid({self._row}, {self._col})"


class NamedQid(Qid):
	"""A qid identified by a name."""

	def __init__(self, name: str) -> None:
		if not isinstance(name, str):
			raise TypeError(f"NamedQid name must be a str, got {type(name).__name__}")
		self._name = name

	@property
	def name(self) -> str:
		return self._name

	def _comparison_key(self) -> str:
		return self._name

	def __str__(self) -> str:
		return self._name

	def __repr__(self) -> str:
		return f"vqubit.NamedQid({self._name!r})"
