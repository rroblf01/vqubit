"""Symbolic parameters, the dependency-free stand-in for Cirq's sympy symbols.

A :class:`Parameter` is a named variable. It supports arithmetic with numbers
and other parameters, producing small expression trees, so a gate angle can be
written as ``Parameter("theta") / 2`` and filled in later by a
:class:`ParamResolver`.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import Any

from ._value import value_equality

_NUMBER = (int, float, complex)


class _Expression:
	"""A node in a parameter expression tree."""

	def names(self) -> set[str]:
		raise NotImplementedError

	def substitute(self, resolver: ParamResolver) -> Any:
		raise NotImplementedError

	def __add__(self, other: Any) -> _Expression:
		return _BinaryOp("+", self, other)

	def __radd__(self, other: Any) -> _Expression:
		return _BinaryOp("+", other, self)

	def __sub__(self, other: Any) -> _Expression:
		return _BinaryOp("-", self, other)

	def __rsub__(self, other: Any) -> _Expression:
		return _BinaryOp("-", other, self)

	def __mul__(self, other: Any) -> _Expression:
		return _BinaryOp("*", self, other)

	def __rmul__(self, other: Any) -> _Expression:
		return _BinaryOp("*", other, self)

	def __truediv__(self, other: Any) -> _Expression:
		return _BinaryOp("/", self, other)

	def __rtruediv__(self, other: Any) -> _Expression:
		return _BinaryOp("/", other, self)

	def __pow__(self, other: Any) -> _Expression:
		return _BinaryOp("**", self, other)

	def __neg__(self) -> _Expression:
		return _UnaryOp("-", self)


@value_equality
class Parameter(_Expression):
	"""A named symbolic parameter."""

	def __init__(self, name: str) -> None:
		if not isinstance(name, str) or not name:
			raise ValueError("Parameter name must be a non-empty str")
		self._name = name

	@property
	def name(self) -> str:
		return self._name

	def names(self) -> set[str]:
		return {self._name}

	def substitute(self, resolver: ParamResolver) -> Any:
		return resolver.value_of(self)

	def _value_equality_values_(self) -> Any:
		return (self._name,)

	def __str__(self) -> str:
		return self._name

	def __repr__(self) -> str:
		return f"vqubit.Parameter({self._name!r})"


@value_equality
class _BinaryOp(_Expression):
	def __init__(self, op: str, left: Any, right: Any) -> None:
		self._op = op
		self._left = left
		self._right = right

	def names(self) -> set[str]:
		return _names(self._left) | _names(self._right)

	def substitute(self, resolver: ParamResolver) -> Any:
		left = _substitute(self._left, resolver)
		right = _substitute(self._right, resolver)
		if isinstance(left, _NUMBER) and isinstance(right, _NUMBER):
			return _apply(self._op, left, right)
		return _BinaryOp(self._op, left, right)

	def _value_equality_values_(self) -> Any:
		return (self._op, self._left, self._right)

	def __str__(self) -> str:
		return f"({self._left} {self._op} {self._right})"

	def __repr__(self) -> str:
		return f"({self._left!r} {self._op} {self._right!r})"


@value_equality
class _UnaryOp(_Expression):
	def __init__(self, op: str, operand: Any) -> None:
		self._op = op
		self._operand = operand

	def names(self) -> set[str]:
		return _names(self._operand)

	def substitute(self, resolver: ParamResolver) -> Any:
		operand = _substitute(self._operand, resolver)
		if isinstance(operand, _NUMBER):
			return -operand
		return _UnaryOp(self._op, operand)

	def _value_equality_values_(self) -> Any:
		return (self._op, self._operand)

	def __str__(self) -> str:
		return f"(-{self._operand})"

	def __repr__(self) -> str:
		return f"({self._op}{self._operand!r})"


def _apply(op: str, left: Any, right: Any) -> Any:
	"""Evaluate one binary operation.

	The operations are dispatched one at a time on purpose: building a dict of
	all of them evaluates ``left / right`` too, which raises ``ZeroDivisionError``
	for any expression whose other operand happens to be ``0``, so resolving a
	parameter *to zero* used to fail.
	"""
	if op == "+":
		return left + right
	if op == "-":
		return left - right
	if op == "*":
		return left * right
	if op == "/":
		return left / right
	if op == "**":
		return left**right
	raise ValueError(f"unknown operation {op!r}")


def _names(value: Any) -> set[str]:
	if isinstance(value, _Expression):
		return value.names()
	return set()


def _substitute(value: Any, resolver: ParamResolver) -> Any:
	if isinstance(value, _Expression):
		return value.substitute(resolver)
	return value


def is_symbolic(value: Any) -> bool:
	"""Report whether ``value`` is a parameter or an expression over parameters."""
	return isinstance(value, _Expression)


def parameter_names(value: Any) -> set[str]:
	"""Return the names of the parameters appearing in ``value``."""
	return _names(value)


def resolve_value(value: Any, resolver: ParamResolver) -> Any:
	"""Substitute parameters in ``value`` using ``resolver``."""
	return _substitute(value, resolver)


class ParamResolver(Mapping[str, Any]):
	"""A mapping from parameter names to values."""

	def __init__(self, param_dict: Mapping[Any, Any] | None = None, **kwargs: Any) -> None:
		self._map: dict[str, Any] = {}
		if param_dict is not None:
			for key, value in dict(param_dict).items():
				self._map[key.name if isinstance(key, Parameter) else str(key)] = value
		for key, value in kwargs.items():
			self._map[key] = value

	def value_of(self, value: Any) -> Any:
		"""Return the value of a parameter, or ``value`` itself if it is not one."""
		if isinstance(value, Parameter):
			return self._map.get(value.name, value)
		if isinstance(value, _Expression):
			return value.substitute(self)
		return value

	def __getitem__(self, key: str | Parameter) -> Any:
		name = key.name if isinstance(key, Parameter) else key
		return self._map[name]

	def __iter__(self) -> Iterator[str]:
		return iter(self._map)

	def __len__(self) -> int:
		return len(self._map)

	def __repr__(self) -> str:
		return f"vqubit.ParamResolver({self._map!r})"
