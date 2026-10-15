"""Extra ``cirq.value`` helpers: angle/bit utilities, durations, conditions,
classical data stores and state constants.
"""

from __future__ import annotations

import math
from abc import ABCMeta
from collections.abc import Iterable, Iterator, Mapping, MutableMapping, Sequence
from typing import Any

from ._value import value_equality

TParamKey = Any
TParamVal = Any
TParamValComplex = Any
DURATION_LIKE = Any
RANDOM_STATE_OR_SEED_LIKE = Any
PRNG_OR_SEED_LIKE = Any

MEASUREMENT_KEY_SEPARATOR = "."

_TWO_PI = 2 * math.pi


def big_endian_bits_to_int(bits: Iterable[int]) -> int:
	"""Return the integer whose big-endian bits are ``bits``."""
	value = 0
	for bit in bits:
		value = (value << 1) | (1 if bit else 0)
	return value


def big_endian_int_to_bits(value: int, *, bit_count: int) -> list[int]:
	"""Return the ``bit_count`` big-endian bits of ``value``."""
	return [(value >> (bit_count - 1 - i)) & 1 for i in range(bit_count)]


def big_endian_digits_to_int(digits: Iterable[int], *, base: int | Sequence[int]) -> int:
	"""Return the integer whose big-endian digits are ``digits`` in ``base``."""
	value = 0
	for index, digit in enumerate(digits):
		digit_base = base[index] if isinstance(base, Sequence) else base
		value = value * digit_base + digit
	return value


def big_endian_int_to_digits(value: int, *, base: int | Sequence[int]) -> list[int]:
	"""Return the big-endian digits of ``value`` in ``base``."""
	if isinstance(base, int):
		length = 0
		probe = value
		while probe:
			length += 1
			probe //= base
		length = max(length, 1)
		bases = [base] * length
	else:
		bases = list(base)
	digits = [0] * len(bases)
	remaining = value
	for index in range(len(bases) - 1, -1, -1):
		digits[index] = remaining % bases[index]
		remaining //= bases[index]
	return digits


def canonicalize_half_turns(half_turns: float) -> float:
	"""Shift ``half_turns`` into the interval ``(-1, 1]``."""
	half_turns = half_turns % 2
	if half_turns > 1:
		half_turns -= 2
	return half_turns


def chosen_angle_to_half_turns(angle: float, *, period: float = _TWO_PI) -> float:
	"""Convert an angle to half turns for the given period."""
	return 2 * angle / period


def chosen_angle_to_canonical_half_turns(angle: float, *, period: float = _TWO_PI) -> float:
	"""Convert an angle to canonical half turns for the given period."""
	return canonicalize_half_turns(chosen_angle_to_half_turns(angle, period=period))


def state_vector_to_probabilities(state_vector: Sequence[complex]) -> list[float]:
	"""Return the measurement probabilities of a state vector."""
	return [abs(complex(amplitude)) ** 2 for amplitude in state_vector]


def validate_probability(probability: float, prob_str: str | None = None) -> float:
	"""Return ``probability`` if it is in ``[0, 1]``, else raise ``ValueError``."""
	probability = float(probability)
	if probability < 0 or probability > 1:
		name = prob_str if prob_str is not None else str(probability)
		raise ValueError(f"The probability {name} was not in the range [0, 1].")
	return probability


def alternative(*args: Any, default: Any = None) -> Any:
	"""Return the first argument that is not ``None``."""
	for arg in args:
		if arg is not None:
			return arg
	return default


class ABCMetaImplementAnyOneOf(ABCMeta):
	"""A metaclass requiring exactly one of a set of abstract methods.

	This is a light version of Cirq's helper; it records the abstract method
	names but does not enforce the constraint at instantiation time.
	"""

	_abstract_one_of: tuple[str, ...] = ()


class LinearDict(MutableMapping):
	"""A mapping that keeps its entries in insertion order and drops zeros."""

	def __init__(self, values: Mapping[Any, Any] | Iterable[tuple[Any, Any]] | None = None) -> None:
		self._items: dict[Any, Any] = {}
		if values is not None:
			for key, value in values.items() if isinstance(values, Mapping) else values:
				if value:
					self._items[key] = value

	def __getitem__(self, key: Any) -> Any:
		return self._items[key]

	def __setitem__(self, key: Any, value: Any) -> None:
		if value:
			self._items[key] = value
		else:
			self._items.pop(key, None)

	def __delitem__(self, key: Any) -> None:
		del self._items[key]

	def __iter__(self) -> Iterator[Any]:
		return iter(self._items)

	def __len__(self) -> int:
		return len(self._items)

	def __add__(self, other: Mapping[Any, Any]) -> LinearDict:
		result = LinearDict(self._items)
		for key, value in other.items():
			result[key] = result.get(key, 0) + value
		return result

	def copy(self) -> LinearDict:
		return LinearDict(self._items)

	def __repr__(self) -> str:
		return f"vqubit.LinearDict({self._items!r})"


@value_equality
class Duration:
	"""A duration of time, stored as an integer number of picoseconds."""

	def __init__(self, *, picos: int) -> None:
		if picos < 0:
			raise ValueError("Duration cannot be negative")
		self._picos = int(picos)

	@property
	def picos(self) -> int:
		return self._picos

	def total_picos(self) -> int:
		return self._picos

	def total_nanos(self) -> float:
		return self._picos / 1000

	def __add__(self, other: Duration) -> Duration:
		return Duration(picos=self._picos + other._picos)

	def __sub__(self, other: Duration) -> Duration:
		return Duration(picos=self._picos - other._picos)

	def __mul__(self, factor: float) -> Duration:
		return Duration(picos=round(self._picos * factor))

	def __rmul__(self, factor: float) -> Duration:
		return self.__mul__(factor)

	def _value_equality_values_(self) -> Any:
		return (self._picos,)

	def __str__(self) -> str:
		return f"{self._picos}ps"

	def __repr__(self) -> str:
		return f"vqubit.Duration(picos={self._picos})"


@value_equality
class Timestamp:
	"""A point in time, stored as an integer number of picoseconds."""

	def __init__(self, *, picos: int) -> None:
		self._picos = int(picos)

	def raw_picos(self) -> int:
		return self._picos

	def __add__(self, other: Duration) -> Timestamp:
		return Timestamp(picos=self._picos + other.picos)

	def __sub__(self, other: Timestamp) -> Duration:
		return Duration(picos=self._picos - other._picos)

	def _value_equality_values_(self) -> Any:
		return (self._picos,)

	def __repr__(self) -> str:
		return f"vqubit.Timestamp(picos={self._picos})"


@value_equality
class PeriodicValue:
	"""A value that is periodic with a given period."""

	def __init__(self, value: float, period: float) -> None:
		self.value = float(value)
		self.period = float(period)

	def _value_equality_values_(self) -> Any:
		return (self.value % self.period, self.period)

	def __repr__(self) -> str:
		return f"vqubit.PeriodicValue({self.value!r}, {self.period!r})"


import enum  # noqa: E402


class MeasurementType(enum.Enum):
	"""Whether a measurement is a Pauli-Z measurement or a Pauli-product one."""

	Z = "Z"
	P = "P"


class BitMaskKeyCondition:
	"""A condition on a bitmask of a measurement result."""

	def __init__(self, key: Any, mask: int) -> None:
		self.key = str(key)
		self.mask = int(mask)

	def __invert__(self) -> BitMaskKeyCondition:
		return BitMaskKeyCondition(self.key, self.mask)

	def __repr__(self) -> str:
		return f"vqubit.BitMaskKeyCondition({self.key!r}, {self.mask})"


class SympyCondition:
	"""A condition on a symbolic expression."""

	def __init__(self, sympy_expr: Any, sympy_symbol_map: Mapping[Any, Any] | None = None) -> None:
		self.sympy_expr = sympy_expr
		self.sympy_symbol_map = dict(sympy_symbol_map) if sympy_symbol_map is not None else {}

	def __repr__(self) -> str:
		return f"vqubit.SympyCondition({self.sympy_expr!r})"


class ClassicalDataStoreReader:
	"""Base class for reading from a classical data store."""

	def keys(self) -> Iterable[str]:
		raise NotImplementedError


class ClassicalDataStore(ClassicalDataStoreReader):
	"""Base class for a store of classical measurement data."""

	def record(self, **kwargs: Any) -> None:
		raise NotImplementedError


class ClassicalDataDictionaryStore(ClassicalDataStore):
	"""A classical data store backed by a dictionary of key -> value list."""

	def __init__(
		self,
		measurements: Mapping[str, Sequence[Any]] | None = None,
		records: Mapping[str, Sequence[Any]] | None = None,
	) -> None:
		self._measurements: dict[str, list[Any]] = {key: list(values) for key, values in (measurements or {}).items()}
		self._records: dict[str, list[Any]] = {key: list(values) for key, values in (records or {}).items()}

	def keys(self) -> Iterable[str]:
		return self._measurements.keys()

	def records(self, key: str) -> list[Any]:
		return list(self._records.get(key, []))

	def measurements(self, key: str) -> list[Any]:
		return list(self._measurements.get(key, []))

	def record(self, **kwargs: Any) -> None:
		for key, value in kwargs.items():
			self._records.setdefault(key, []).append(value)

	def __repr__(self) -> str:
		return f"vqubit.ClassicalDataDictionaryStore(keys={list(self._measurements)!r})"


KET_ZERO: list[list[complex]] = [[1 + 0j], [0j]]
KET_ONE: list[list[complex]] = [[0j], [1 + 0j]]
KET_PLUS: list[list[complex]] = [[1 / math.sqrt(2)], [1 / math.sqrt(2)]]
KET_MINUS: list[list[complex]] = [[1 / math.sqrt(2)], [-1 / math.sqrt(2)]]
KET_IMAG: list[list[complex]] = [[1 / math.sqrt(2)], [1j / math.sqrt(2)]]
KET_MINUS_IMAG: list[list[complex]] = [[1 / math.sqrt(2)], [-1j / math.sqrt(2)]]

PAULI_STATES: dict[Any, list[list[complex]]] = {
	"I": KET_ZERO,
	"X": KET_PLUS,
	"Y": KET_IMAG,
	"Z": KET_ZERO,
}


def _flatten_state(state: Any) -> list[complex]:
	if isinstance(state, Sequence) and state and isinstance(state[0], Sequence):
		return [complex(z) for row in state for z in row]
	return [complex(z) for z in state]


@value_equality
class ProductState:
	"""A tensor product of single-qubit states."""

	def __init__(self, state_vector: Sequence[Any]) -> None:
		self._states = [_flatten_state(state) for state in state_vector]

	def __getitem__(self, index: int) -> list[complex]:
		return list(self._states[index])

	def __len__(self) -> int:
		return len(self._states)

	@property
	def states(self) -> list[list[complex]]:
		return [list(state) for state in self._states]

	def state_vector(self) -> list[complex]:
		vector = [1 + 0j]
		for state in self._states:
			vector = [a * b for a in vector for b in state]
		return vector

	@staticmethod
	def _flatten(state: Any) -> list[complex]:
		return _flatten_state(state)

	def _value_equality_values_(self) -> Any:
		return (tuple(tuple(state) for state in self._states),)

	def __repr__(self) -> str:
		return f"vqubit.ProductState({self._states!r})"


__all__ = [
	"DURATION_LIKE",
	"KET_IMAG",
	"KET_MINUS",
	"KET_MINUS_IMAG",
	"KET_ONE",
	"KET_PLUS",
	"KET_ZERO",
	"MEASUREMENT_KEY_SEPARATOR",
	"PAULI_STATES",
	"PRNG_OR_SEED_LIKE",
	"RANDOM_STATE_OR_SEED_LIKE",
	"ABCMetaImplementAnyOneOf",
	"BitMaskKeyCondition",
	"ClassicalDataDictionaryStore",
	"ClassicalDataStore",
	"ClassicalDataStoreReader",
	"Duration",
	"LinearDict",
	"MeasurementType",
	"PeriodicValue",
	"ProductState",
	"SympyCondition",
	"TParamKey",
	"TParamVal",
	"TParamValComplex",
	"Timestamp",
	"alternative",
	"big_endian_bits_to_int",
	"big_endian_digits_to_int",
	"big_endian_int_to_bits",
	"big_endian_int_to_digits",
	"canonicalize_half_turns",
	"chosen_angle_to_canonical_half_turns",
	"chosen_angle_to_half_turns",
	"state_vector_to_probabilities",
	"validate_probability",
]
