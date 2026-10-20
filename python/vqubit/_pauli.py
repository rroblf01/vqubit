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
import math
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import Any, TypeVar

from ._compat import asarray
from ._gates import Gate, I, X, Y, Z
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


#: What a Pauli string's values may be, as Cirq accepts them: the enum, its gate
#: (``vq.X``, and any power that equals it) or its letter.
PauliLike = Pauli | Gate | str

#: The key of a qubit map, left open as Cirq leaves it: ``Mapping`` is invariant
#: in its key type, so ``Mapping[Qid, int]`` would refuse a ``dict[LineQubit, int]``.
TKey = TypeVar("TKey")

# A product table keyed by enum members, so multiplying strings avoids string
# comparisons in the hot loop.
_PRODUCT_MAP: dict[tuple[Pauli, Pauli], tuple[Pauli, complex]] = {}
for _left in Pauli:
	for _right in Pauli:
		_result, _phase = _product(_left.value, _right.value)
		_PRODUCT_MAP[(_left, _right)] = (Pauli(_result), _phase)
del _left, _right, _result, _phase


def _flatten_state(values: Any) -> list[complex]:
	"""Flatten a flat or nested state vector (list or numpy array) to amplitudes."""
	if hasattr(values, "ravel"):
		return [complex(entry) for entry in values.ravel()]
	flat: list[complex] = []
	for value in values:
		if isinstance(value, (list, tuple)):
			flat.extend(_flatten_state(value))
		else:
			flat.append(complex(value))
	return flat


def _num_qubits_of(size: int) -> int:
	"""``log2(size)`` for a power of two, the qubit count of a state of that size."""
	if size == 0 or size & (size - 1):
		raise ValueError(f"the state has {size} entries, which is not a power of two")
	return size.bit_length() - 1


def _validate_qubit_map(qubit_map: Mapping[Any, int], qubits: Sequence[Qid], num_qubits: int) -> None:
	"""Refuse a map that is out of range or does not cover every qubit, as Cirq does."""
	for index in qubit_map.values():
		if not 0 <= index < num_qubits:
			raise ValueError(f"the qubit map indices must be valid for a state over {num_qubits} qubits")
	missing = set(qubits) - set(qubit_map)
	if missing:
		raise ValueError(f"the qubit map must cover every qubit of the Pauli string: {sorted(missing)!r} is missing")


def _check_hermitian(coefficient: complex, pauli_string: PauliString) -> None:
	"""Refuse a non-Hermitian string, whose ``<P>`` is not what a user means.

	Expectation values are defined for Hermitian operators, which is what Cirq
	checks (and the message it gives).
	"""
	if abs(coefficient.imag) > 1e-4:
		raise NotImplementedError(
			f"cannot compute the expectation value of the non-Hermitian PauliString <{pauli_string}>; "
			"the coefficient must be real"
		)


def _check_normalized(amplitudes: Sequence[complex], atol: float) -> None:
	"""Refuse a state vector that is not normalized, as Cirq does."""
	norm = sum(abs(amplitude) ** 2 for amplitude in amplitudes) ** 0.5
	if abs(norm - 1) > atol:
		raise ValueError(f"the state vector is not normalized: its norm is {norm}")


def _check_density_matrix(trace: complex, atol: float) -> None:
	"""Refuse a density matrix that does not have unit trace, as Cirq refuses one."""
	if abs(trace - 1) > atol:
		raise ValueError(f"the density matrix is not normalized: its trace is {trace}")


def _as_pauli(value: PauliLike) -> Pauli:
	"""Accept what Cirq accepts for a Pauli string's values.

	That is a :class:`Pauli`, one of its gates (``vq.X``, ``vq.Y``, ``vq.Z``,
	``vq.I``, and any power that *equals* one of them) or its letter, so the Cirq
	idiom ``vq.PauliString({q0: vq.Z})`` works.
	"""
	if isinstance(value, Pauli):
		return value
	if isinstance(value, str):
		try:
			return Pauli(value.upper())
		except ValueError as exc:
			raise TypeError(f"{value!r} is not a Pauli") from exc
	for letter, gate in (("X", X), ("Y", Y), ("Z", Z), ("I", I)):
		if value == gate:
			return Pauli(letter)
	raise TypeError(f"PauliString values must be a Pauli or a Pauli gate, got {type(value).__name__}")


@value_equality
class PauliString:
	"""A Pauli on each of a set of qubits, with a complex coefficient."""

	def __init__(
		self,
		qubit_pauli_map: Mapping[Any, PauliLike] | Iterable[tuple[Any, PauliLike]] | None = None,
		coefficient: complex = 1,
	) -> None:
		items = dict(qubit_pauli_map) if qubit_pauli_map is not None else {}
		self._map: dict[Qid, Pauli] = {}
		for qubit, pauli in items.items():
			if not isinstance(qubit, Qid):
				raise TypeError(f"PauliString keys must be Qid, got {type(qubit).__name__}")
			value = _as_pauli(pauli)
			if value is not Pauli.I:
				self._map[qubit] = value
		self._coefficient = complex(coefficient)

	@property
	def coefficient(self) -> complex:
		"""The complex coefficient multiplying the string."""
		return self._coefficient

	def items(self) -> Iterator[tuple[Qid, Pauli]]:
		"""Iterate ``(qubit, pauli)`` pairs in qubit order."""
		return iter(sorted(self._map.items()))

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits the string touches, sorted."""
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
		other_map = other._map
		for qubit, left in self._map.items():
			right = other_map.get(qubit)
			if right is None:
				result[qubit] = left
				continue
			product, factor = _PRODUCT_MAP[(left, right)]
			phase *= factor
			if product is not Pauli.I:
				result[qubit] = product
		for qubit, right in other_map.items():
			if qubit not in self._map:
				result[qubit] = right
		return PauliString(result, phase)

	def __rmul__(self, other: Any) -> PauliString:
		if isinstance(other, complex | int | float):
			return PauliString(self._map, other * self._coefficient)
		return NotImplemented

	def __neg__(self) -> PauliString:
		return PauliString(self._map, -self._coefficient)

	def __pow__(self, exponent: int) -> PauliString:
		# Every Pauli squares to the identity, so only the parity of the
		# exponent matters: ``X ** -1`` is ``X``, not the identity.
		coefficient = self._coefficient**exponent
		if exponent % 2 == 0:
			return PauliString(coefficient=coefficient)
		return PauliString(self._map, coefficient)

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

	def __bool__(self) -> bool:
		return bool(self._map)

	def keys(self) -> Any:
		"""The qubits the string touches."""
		return self._map.keys()

	def values(self) -> Iterator[Pauli]:
		"""The Pauli on each qubit the string touches."""
		return iter([pauli for _, pauli in self.items()])

	def get(self, key: Qid, default: Any = None) -> Any:
		"""The Pauli on ``key``, or ``default`` when the string does not touch it."""
		return self._map.get(key, Pauli.I if default is None else default)

	def __iter__(self) -> Iterator[Qid]:
		"""Iterate the qubits the string touches, as Cirq does."""
		return iter(self.qubits)

	def __len__(self) -> int:
		return len(self._map)

	def __truediv__(self, other: Any) -> PauliString:
		if isinstance(other, complex | int | float):
			return PauliString(self._map, self._coefficient / other)
		return NotImplemented

	def __add__(self, other: Any) -> Any:
		"""A sum, as in ``cirq.PauliString + ...``: ``ps + ps`` is a ``PauliSum``."""
		if not isinstance(other, PauliSum | PauliString | complex | int | float):
			return NotImplemented
		return PauliSum() + self + other

	def __radd__(self, other: Any) -> Any:
		if not isinstance(other, PauliSum | PauliString | complex | int | float):
			return NotImplemented
		return PauliSum() + other + self

	def __sub__(self, other: Any) -> Any:
		if not isinstance(other, PauliSum | PauliString | complex | int | float):
			return NotImplemented
		return PauliSum() + self - other

	def __rsub__(self, other: Any) -> Any:
		if not isinstance(other, PauliSum | PauliString | complex | int | float):
			return NotImplemented
		return PauliSum() + other - self

	def to_z_basis_ops(self) -> Iterator[Any]:
		"""Return single-qubit operations that rotate each Pauli to the Z basis."""
		from ._clifford_gates import SingleQubitCliffordGate

		for qubit, pauli in self.items():
			yield SingleQubitCliffordGate.from_single_map({pauli: (Pauli.Z, False)}).on(qubit)

	def _expectation_masks(self, qubit_map: Mapping[Any, int], num_qubits: int) -> tuple[int, int, int]:
		"""``(bits X flips, bits that flip the sign, number of Y)`` for this string.

		The Paulis act on a computational basis state ``|x>`` as
		``P|x> = phase(x) |x XOR flipped>``, with ``phase(x) = i^(#Y) (-1)^(bits
		of x in the Z and Y positions)``, which is all the expectation value needs.

		Cirq's map gives each qubit its *axis* in the state, with axis ``0`` the
		most significant bit and the last qubit the least significant one, so the
		axis is turned into a bit position here. Getting that backwards is
		invisible on a symmetric state such as a Bell pair, so it is exactly what
		the harness compares with a non-symmetric one.
		"""
		flipped = 0
		signs = 0
		num_y = 0
		for qubit, pauli in self._map.items():
			index = num_qubits - 1 - qubit_map[qubit]
			if pauli is not Pauli.Z:
				flipped |= 1 << index
			if pauli is not Pauli.X:
				signs |= 1 << index
			if pauli is Pauli.Y:
				num_y += 1
		return flipped, signs, num_y

	def _expectation_from_amplitudes(
		self, amplitudes: Sequence[complex], qubit_map: Mapping[Any, int], num_qubits: int
	) -> complex:
		"""``<psi|P|psi>``, walking the amplitudes ``P`` pairs up."""
		flipped, signs, num_y = self._expectation_masks(qubit_map, num_qubits)
		phase = 1j**num_y
		total = 0j
		for index, amplitude in enumerate(amplitudes):
			sign = -1 if (index & signs).bit_count() % 2 else 1
			total += amplitudes[index ^ flipped].conjugate() * phase * sign * amplitude
		return self._coefficient * total

	def _expectation_from_matrix(
		self, matrix: Sequence[complex], dim: int, qubit_map: Mapping[Any, int], num_qubits: int
	) -> complex:
		"""``Tr(rho P)``, the same walk down a flattened density matrix."""
		flipped, signs, num_y = self._expectation_masks(qubit_map, num_qubits)
		phase = (1j**num_y).conjugate()
		total = 0j
		for row in range(dim):
			sign = -1 if (row & signs).bit_count() % 2 else 1
			total += matrix[(row ^ flipped) * dim + row] * phase * sign
		return self._coefficient * total

	def expectation_from_state_vector(
		self,
		state_vector: Any,
		qubit_map: Mapping[TKey, int],
		*,
		atol: float = 1e-7,
		check_preconditions: bool = True,
	) -> complex:
		"""Return ``<psi|P|psi>`` for a state vector expressed over ``qubit_map``.

		``qubit_map`` says which index of the state each qubit of this string
		holds, so the same state can be read in any qubit order, as in Cirq. The
		expectation value is defined for Hermitian operators, so a non-Hermitian
		coefficient raises, and the state must be normalized unless
		``check_preconditions`` is ``False``.
		"""
		_check_hermitian(self._coefficient, self)
		amplitudes = _flatten_state(state_vector)
		num_qubits = _num_qubits_of(len(amplitudes))
		_validate_qubit_map(qubit_map, self.qubits, num_qubits)
		if check_preconditions:
			_check_normalized(amplitudes, atol)
		return self._expectation_from_amplitudes(amplitudes, qubit_map, num_qubits)

	def expectation_from_density_matrix(
		self,
		state: Any,
		qubit_map: Mapping[TKey, int],
		*,
		atol: float = 1e-7,
		check_preconditions: bool = True,
	) -> complex:
		"""Return ``Tr(rho P)`` for a density matrix expressed over ``qubit_map``."""
		_check_hermitian(self._coefficient, self)
		matrix = _flatten_state(state)
		dim = math.isqrt(len(matrix))
		if dim * dim != len(matrix):
			raise ValueError(f"the array has {len(matrix)} entries, which is not 2**n by 2**n")
		num_qubits = _num_qubits_of(dim)
		_validate_qubit_map(qubit_map, self.qubits, num_qubits)
		if check_preconditions:
			_check_density_matrix(sum(matrix[step * dim + step] for step in range(dim)), atol)
		return self._expectation_from_matrix(matrix, dim, qubit_map, num_qubits)

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

	def __init__(self, terms: Any = None) -> None:
		self._terms: dict[PauliString, complex] = {}
		if terms is not None:
			for pauli_string, coefficient in dict(terms).items():
				if not isinstance(pauli_string, PauliString):
					raise TypeError("PauliSum keys must be PauliString instances")
				if abs(coefficient) > 0:
					self._terms[pauli_string] = complex(coefficient)

	def __eq__(self, other: object) -> bool:
		if not isinstance(other, PauliSum):
			return NotImplemented
		return self._terms == other._terms

	def __ne__(self, other: object) -> bool:
		result = self.__eq__(other)
		return result if result is NotImplemented else not result

	def __hash__(self) -> int:
		return hash(frozenset(self._terms.items()))

	@classmethod
	def wrap(cls, value: Any) -> PauliSum:
		"""Convert a ``PauliSumLike`` (number, PauliString or PauliSum) to a ``PauliSum``."""
		if isinstance(value, PauliSum):
			return value
		return cls() + value

	@classmethod
	def from_pauli_strings(cls, terms: Any) -> PauliSum:
		"""Build a ``PauliSum`` from one or more ``PauliString`` terms."""
		if isinstance(terms, PauliString):
			terms = [terms]
		result = cls()
		for term in terms:
			result = result + term
		return result

	def __add__(self, other: Any) -> PauliSum:
		terms = dict(self._terms)
		if isinstance(other, PauliSum):
			items: Any = other._terms.items()
		elif isinstance(other, PauliString):
			items = [(other, other.coefficient)]
		elif isinstance(other, complex | int | float):
			items = [(PauliString(), complex(other))]
		else:
			return NotImplemented
		for pauli_string, coefficient in items:
			unit = PauliString(pauli_string._map)
			terms[unit] = terms.get(unit, 0j) + coefficient
		return PauliSum(terms)

	def __radd__(self, other: Any) -> PauliSum:
		return self + other

	def __neg__(self) -> PauliSum:
		return PauliSum({pauli_string: -coefficient for pauli_string, coefficient in self._terms.items()})

	def __sub__(self, other: Any) -> PauliSum:
		if isinstance(other, PauliSum | PauliString):
			return self + (-other)
		if isinstance(other, complex | int | float):
			return self + (-complex(other))
		return NotImplemented

	def __rsub__(self, other: Any) -> PauliSum:
		return (-self) + other

	def __mul__(self, other: Any) -> PauliSum:
		"""Multiply by a number (scaling every term) or by a Pauli string."""
		if isinstance(other, complex | int | float):
			return PauliSum({pauli_string: coefficient * other for pauli_string, coefficient in self._terms.items()})
		if isinstance(other, PauliString):
			result = PauliSum()
			for pauli_string, coefficient in self._terms.items():
				result = result + PauliString(pauli_string._map, coefficient) * other
			return result
		if isinstance(other, PauliSum):
			result = PauliSum()
			for left, left_coefficient in self._terms.items():
				for right, right_coefficient in other._terms.items():
					product = PauliString(left._map, left_coefficient) * PauliString(right._map, right_coefficient)
					result = result + product
			return result
		return NotImplemented

	def __rmul__(self, other: Any) -> PauliSum:
		return self * other

	def __truediv__(self, other: Any) -> PauliSum:
		if isinstance(other, complex | int | float):
			return PauliSum({pauli_string: coefficient / other for pauli_string, coefficient in self._terms.items()})
		return NotImplemented

	def __getitem__(self, pauli_string: PauliString) -> complex:
		return self._terms.get(pauli_string, 0j)

	def items(self) -> Iterator[tuple[PauliString, complex]]:
		"""Iterate ``(pauli_string, coefficient)`` pairs."""
		return iter(self._terms.items())

	def __iter__(self) -> Iterator[PauliString]:
		for pauli_string, coefficient in self._terms.items():
			yield PauliString(pauli_string._map, coefficient)

	def expectation_from_state_vector(
		self,
		state_vector: Any,
		qubit_map: Mapping[TKey, int],
		*,
		atol: float = 1e-7,
		check_preconditions: bool = True,
	) -> complex:
		"""Return ``<psi|S|psi>`` for the whole sum, one term at a time."""
		amplitudes = _flatten_state(state_vector)
		num_qubits = _num_qubits_of(len(amplitudes))
		if check_preconditions:
			_check_normalized(amplitudes, atol)
		total = 0j
		for pauli_string, coefficient in self._terms.items():
			_check_hermitian(coefficient * pauli_string.coefficient, pauli_string)
			_validate_qubit_map(qubit_map, pauli_string.qubits, num_qubits)
			total += coefficient * pauli_string._expectation_from_amplitudes(amplitudes, qubit_map, num_qubits)
		return total

	def expectation_from_density_matrix(
		self,
		state: Any,
		qubit_map: Mapping[TKey, int],
		*,
		atol: float = 1e-7,
		check_preconditions: bool = True,
	) -> complex:
		"""Return ``Tr(rho S)`` for the whole sum, one term at a time."""
		matrix = _flatten_state(state)
		dim = math.isqrt(len(matrix))
		if dim * dim != len(matrix):
			raise ValueError(f"the array has {len(matrix)} entries, which is not 2**n by 2**n")
		num_qubits = _num_qubits_of(dim)
		if check_preconditions:
			_check_density_matrix(sum(matrix[step * dim + step] for step in range(dim)), atol)
		total = 0j
		for pauli_string, coefficient in self._terms.items():
			_check_hermitian(coefficient * pauli_string.coefficient, pauli_string)
			_validate_qubit_map(qubit_map, pauli_string.qubits, num_qubits)
			total += coefficient * pauli_string._expectation_from_matrix(matrix, dim, qubit_map, num_qubits)
		return total

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
		"""The Pauli on every qubit, including the identities."""
		return self._paulis

	@property
	def coefficient(self) -> complex:
		"""The complex coefficient multiplying the dense string."""
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
