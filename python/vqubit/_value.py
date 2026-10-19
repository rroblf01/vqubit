"""Value equality, the Cirq protocol behind ``==`` for immutable values.

Gates, qubits and operations are *values*: two of them are equal when they are
built from the same components, not when they happen to be the same object. A
class declares its components with ``_value_equality_values_`` and applies the
:func:`value_equality` decorator, which fills in ``__eq__``, ``__ne__`` and
``__hash__``.
"""

from __future__ import annotations

import contextlib
from collections.abc import Callable
from typing import Any, TypeVar, overload

T = TypeVar("T")


@overload
def value_equality(cls: T, *, unhashable: bool = ...) -> T: ...


@overload
def value_equality(*, unhashable: bool = ...) -> Callable[[T], T]: ...


def value_equality(cls: Any = None, *, unhashable: bool = False) -> Any:
	"""Define ``__eq__``, ``__ne__`` and ``__hash__`` from ``_value_equality_values_``.

	Two instances are equal when they are of the same concrete type and their
	``_value_equality_values_`` are equal, which is what makes an
	``XPowGate(exponent=0.5)`` equal to another one and not to a ``YPowGate``.

	Set ``unhashable`` for values that cannot be hashed (for instance because
	they hold a list), which sets ``__hash__`` to ``None`` the way an ordinary
	mutable class does.
	"""

	def wrap(cls: Any) -> Any:
		def _eq(self: Any, other: Any) -> Any:
			if self is other:
				return True
			if type(self) is not type(other):
				return NotImplemented
			return bool(self._value_equality_values_() == other._value_equality_values_())

		def _ne(self: Any, other: Any) -> Any:
			result = _eq(self, other)
			if result is NotImplemented:
				return result
			return not result

		def _hash(self: Any) -> int:
			cached = getattr(self, "_hash_cache", None)
			if cached is None:
				values = self._value_equality_values_()
				if not isinstance(values, tuple):
					values = (values,)
				cached = hash((type(self), *values))
				with contextlib.suppress(AttributeError):  # pragma: no cover - slotted classes
					object.__setattr__(self, "_hash_cache", cached)
			return cached

		cls.__eq__ = _eq
		cls.__ne__ = _ne
		cls.__hash__ = None if unhashable else _hash
		return cls

	if cls is None:
		return wrap
	return wrap(cls)


@value_equality
class PeriodicValue:
	"""A value that is periodic with a given period.

	Cirq's wrapper: the value is normalized into ``[0, period)`` on the way in,
	so ``PeriodicValue(3, 2) == PeriodicValue(1, 2)`` and both read ``1``. The
	power gates compare their eigenphases through it, which is what makes
	``cirq.X ** 3 == cirq.X`` and ``cirq.CNOT ** -1 == cirq.CNOT``.
	"""

	def __init__(self, value: float, period: float) -> None:
		self.period = float(period)
		self.value = float(value) % self.period

	def _value_equality_values_(self) -> Any:
		return (self.value, self.period)

	def __repr__(self) -> str:
		return f"vqubit.PeriodicValue({self.value!r}, {self.period!r})"
