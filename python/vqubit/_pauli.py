"""Pauli operators and Pauli strings.

A :class:`Pauli` is one of ``I``, ``X``, ``Y`` or ``Z``. A :class:`PauliString`
attaches a Pauli to each of a collection of qubits, with an overall complex
coefficient -- the compact way to write operators like ``X(q0)*Z(q1)`` or the
stabilizers of a code.

Pauli strings are values: two are equal when they carry the same coefficient on
the same qubits, which is what makes them usable as dictionary keys.
"""

from __future__ import annotations

import enum
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import Any

from ._compat import asarray
from ._gates import I, X, Y, Z
from ._linalg import Matrix, kron
from ._qubits import Qid
from ._value import value_equality

_PAULI_MATRIX: dict[str, Matrix] = {
	"I": [[1 + 0j, 0j], [0j, 1 + 0j]],
	"X": [[0j, 1 + 0j], [1 + 0j, 0j]],
	"Y": [[0j, -1j], [1j, 0j]],
	"Z": [[1 + 0j, 0j], [0j, -1 + 0j]],
}

_PAULI_GATE: dict[str, Any] = {"I": I, "X": X, "Y": Y, "Z": Z}

# Multiplication of two Paulis: the resulting Pauli and the phase. The identity
# row and column are filled in programmatically below.
_PRODUCTS: dict[tuple[str, str], tuple[str, complex]] = {
	("X", "Y"): ("Z", 1j),
	("Y", "X"): ("Z", -1j),
	("Y", "Z"): ("X", 1j),
	("Z", "Y"): ("X", -1j),
	("Z", "X"): ("Y", 1j),
	("X", "Z"): ("Y", -1j),
}


def _product(a: str, b: str) -> tuple[str, complex]:
	if a == "I":
		return (b, 1 + 0j)
	if b == "I":
		return (a, 1 + 0j)
	if a == b:
		return ("I", 1 + 0j)
	return _PRODUCTS[(a, b)]


class Pauli(enum.Enum):
	"""A single-qubit Pauli operator."""

	I = "I"  # noqa: E741 - the name matches the Pauli I notation
	X = "X"
	Y = "Y"
	Z = "Z"

	def matrix(self) -> Matrix:
		"""Return the 2-by-2 matrix of this Pauli."""
		return [list(row) for row in _PAULI_MATRIX[self.value]]

	def on(self, qubit: Qid) -> Any:
		"""Return this Pauli applied to ``qubit`` as an operation."""
		return _PAULI_GATE[self.value].on(qubit)

	def __str__(self) -> str:
		return self.value


@value_equality
class PauliString:
	"""A Pauli on each of a set of qubits, with a complex coefficient."""

	def __init__(
		self,
		qubit_pauli_map: Mapping[Qid, Pauli] | Iterable[tuple[Qid, Pauli]] | None = None,
		coefficient: complex = 1,
	) -> None:
		items = dict(qubit_pauli_map) if qubit_pauli_map is not None else {}
		self._map: dict[Qid, Pauli] = {}
		for qubit, pauli in items.items():
			if not isinstance(qubit, Qid):
				raise TypeError(f"PauliString keys must be Qid, got {type(qubit).__name__}")
			if not isinstance(pauli, Pauli):
				raise TypeError(f"PauliString values must be Pauli, got {type(pauli).__name__}")
			if pauli is not Pauli.I:
				self._map[qubit] = pauli
		self._coefficient = complex(coefficient)

	@property
	def coefficient(self) -> complex:
		return self._coefficient

	def items(self) -> Iterator[tuple[Qid, Pauli]]:
		return iter(sorted(self._map.items()))

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return tuple(sorted(self._map))

	def __getitem__(self, qubit: Qid) -> Pauli:
		return self._map.get(qubit, Pauli.I)

	def __mul__(self, other: Any) -> PauliString:
		if isinstance(other, complex | int | float):
			return PauliString(self._map, self._coefficient * other)
		if not isinstance(other, PauliString):
			return NotImplemented
		phase = self._coefficient * other._coefficient
		result: dict[Qid, Pauli] = {}
		for qubit in set(self._map) | set(other._map):
			product, factor = _product(self[qubit].value, other[qubit].value)
			phase *= factor
			if product != "I":
				result[qubit] = Pauli(product)
		return PauliString(result, phase)

	def __rmul__(self, other: Any) -> PauliString:
		if isinstance(other, complex | int | float):
			return PauliString(self._map, other * self._coefficient)
		return NotImplemented

	def __neg__(self) -> PauliString:
		return PauliString(self._map, -self._coefficient)

	def __pow__(self, exponent: int) -> PauliString:
		coefficient = self._coefficient**exponent
		result: dict[Qid, Pauli] = {}
		for qubit, pauli in self._map.items():
			product, factor = "I", 1 + 0j
			for _ in range(exponent):
				product, step = _product(product, pauli.value)
				factor *= step
			coefficient *= factor
			if product != "I":
				result[qubit] = Pauli(product)
		return PauliString(result, coefficient)

	def matrix(self) -> Matrix:
		"""Return the full matrix over this string's qubits, coefficient included."""
		qubits = self.qubits
		if not qubits:
			return asarray([[self._coefficient]])
		result = [[self._coefficient]]
		for qubit in qubits:
			result = kron(result, self[qubit].matrix())
		return asarray(result)

	def unitary(self) -> Matrix:
		"""Alias of :meth:`matrix`."""
		return self.matrix()

	def commutes_with(self, other: PauliString) -> bool:
		"""Report whether this string commutes with ``other``."""
		return bool((self * other).coefficient == (other * self).coefficient)

	def _value_equality_values_(self) -> Any:
		return (self._coefficient, tuple(sorted(self._map.items())))

	def __str__(self) -> str:
		if not self._map:
			return "I" if self._coefficient == 1 else str(_format_coefficient(self._coefficient))
		body = "*".join(f"{pauli.value}({qubit})" for qubit, pauli in self.items())
		if self._coefficient == 1:
			return body
		return f"{_format_coefficient(self._coefficient)}*{body}"

	def __repr__(self) -> str:
		return f"vqubit.PauliString({self!s})"


def _format_coefficient(coefficient: complex) -> str:
	if coefficient.imag == 0:
		return str(int(coefficient.real)) if coefficient.real == int(coefficient.real) else str(coefficient.real)
	return str(coefficient)


class PauliSum:
	"""A weighted sum of Pauli strings, the output of a Pauli expansion."""

	def __init__(self, terms: Mapping[PauliString, complex] | None = None) -> None:
		self._terms: dict[PauliString, complex] = {}
		if terms is not None:
			for pauli_string, coefficient in dict(terms).items():
				if not isinstance(pauli_string, PauliString):
					raise TypeError("PauliSum keys must be PauliString instances")
				if abs(coefficient) > 0:
					self._terms[pauli_string] = complex(coefficient)

	def __getitem__(self, pauli_string: PauliString) -> complex:
		return self._terms.get(pauli_string, 0j)

	def items(self) -> Iterator[tuple[PauliString, complex]]:
		return iter(self._terms.items())

	def __iter__(self) -> Iterator[PauliString]:
		return iter(self._terms)

	def __len__(self) -> int:
		return len(self._terms)

	def __str__(self) -> str:
		if not self._terms:
			return "0"
		return " + ".join(f"{_format_coefficient(c)}*{ps}" for ps, c in self._terms.items())

	def __repr__(self) -> str:
		return f"vqubit.PauliSum({self})"


@value_equality
class DensePauliString:
	"""A Pauli on every qubit, written densely as a string like ``"IXYZ"``."""

	def __init__(self, pauli_string: str, coefficient: complex = 1) -> None:
		if not isinstance(pauli_string, str):
			raise TypeError("DensePauliString needs a string like 'IXYZ'")
		self._paulis = tuple(Pauli(letter) for letter in pauli_string)
		self._coefficient = complex(coefficient)

	@property
	def paulis(self) -> tuple[Pauli, ...]:
		return self._paulis

	@property
	def coefficient(self) -> complex:
		return self._coefficient

	def __len__(self) -> int:
		return len(self._paulis)

	def __getitem__(self, index: int) -> Pauli:
		return self._paulis[index]

	def _value_equality_values_(self) -> Any:
		return (self._paulis, self._coefficient)

	def to_pauli_string(self, qubits: Sequence[Qid]) -> PauliString:
		"""Return the sparse :class:`~vqubit.PauliString` over ``qubits``."""
		qubit_pauli_map: dict[Qid, Pauli] = {}
		for qubit, pauli in zip(qubits, self._paulis, strict=True):
			if pauli is not Pauli.I:
				qubit_pauli_map[qubit] = pauli
		return PauliString(qubit_pauli_map, coefficient=self._coefficient)

	def matrix(self) -> Matrix:
		"""Return the full matrix, coefficient included."""
		result: Matrix = [[self._coefficient]]
		for pauli in self._paulis:
			result = kron(result, pauli.matrix())
		return asarray(result)

	def __mul__(self, other: Any) -> DensePauliString:
		if isinstance(other, DensePauliString):
			if len(other) != len(self):
				raise ValueError("DensePauliString lengths must match")
			coefficient = self._coefficient * other._coefficient
			paulis: list[Pauli] = []
			for left, right in zip(self._paulis, other._paulis, strict=True):
				product, factor = _product(left.value, right.value)
				coefficient *= factor
				paulis.append(Pauli(product))
			return DensePauliString("".join(p.value for p in paulis), coefficient)
		if isinstance(other, complex | int | float):
			return DensePauliString("".join(p.value for p in self._paulis), self._coefficient * other)
		return NotImplemented

	def __rmul__(self, other: Any) -> DensePauliString:
		if isinstance(other, complex | int | float):
			return DensePauliString("".join(p.value for p in self._paulis), other * self._coefficient)
		return NotImplemented

	def __str__(self) -> str:
		return "".join(pauli.value for pauli in self._paulis)

	def __repr__(self) -> str:
		return f"vqubit.DensePauliString({self!s})"


PAULI_GATE_LIKE = Any
PAULI_STRING_LIKE = Any
PauliSumLike = Any
