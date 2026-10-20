"""Pauli-string operations, phasors, measurement and linear combinations.

These build on the values in :mod:`vqubit._pauli` (``Pauli``, ``PauliString``,
``DensePauliString``) and provide the operation/gate wrappers that Cirq exposes
under ``cirq.ops``.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from ._gates import Gate, I, X, Y, Z
from ._linalg import Matrix
from ._ops import GateOperation, Operation
from ._pauli import DensePauliString, Pauli, PauliString, PauliSum
from ._qubits import Qid
from ._value import value_equality

_PAULI_GATE: dict[Pauli, Any] = {Pauli.I: I, Pauli.X: X, Pauli.Y: Y, Pauli.Z: Z}

BaseDensePauliString = DensePauliString


class MutablePauliString:
	"""A mutable :class:`~vqubit.PauliString`."""

	def __init__(
		self,
		qubit_pauli_map: Mapping[Qid, Pauli] | Iterable[tuple[Qid, Pauli]] | None = None,
		coefficient: complex = 1,
	) -> None:
		items = dict(qubit_pauli_map) if qubit_pauli_map is not None else {}
		self._map: dict[Qid, Pauli] = {q: p for q, p in items.items() if p is not Pauli.I}
		self._coefficient = complex(coefficient)

	@property
	def coefficient(self) -> complex:
		"""The coefficient, which can be reassigned."""
		return self._coefficient

	@coefficient.setter
	def coefficient(self, value: complex) -> None:
		self._coefficient = complex(value)

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits the string currently touches, sorted."""
		return tuple(sorted(self._map))

	def __getitem__(self, qubit: Qid) -> Pauli:
		return self._map.get(qubit, Pauli.I)

	def __setitem__(self, qubit: Qid, pauli: Pauli) -> None:
		if not isinstance(pauli, Pauli):
			raise TypeError("MutablePauliString values must be Pauli")
		if pauli is Pauli.I:
			self._map.pop(qubit, None)
		else:
			self._map[qubit] = pauli

	def __delitem__(self, qubit: Qid) -> None:
		del self._map[qubit]

	def items(self) -> Iterable[tuple[Qid, Pauli]]:
		"""Iterate ``(qubit, pauli)`` pairs in qubit order."""
		return iter(sorted(self._map.items()))

	def to_pauli_string(self) -> PauliString:
		"""Return an immutable copy of the string."""
		return PauliString(self._map, self._coefficient)

	def __str__(self) -> str:
		return str(self.to_pauli_string())

	def __repr__(self) -> str:
		return f"vqubit.MutablePauliString({self})"


class MutableDensePauliString(DensePauliString):
	"""A mutable :class:`~vqubit.DensePauliString`."""

	def __setitem__(self, index: int, pauli: Pauli) -> None:
		paulis = list(self._paulis)
		paulis[index] = pauli
		self._paulis = tuple(paulis)

	@property
	def coefficient(self) -> complex:
		"""The coefficient, which can be reassigned."""
		return self._coefficient

	@coefficient.setter
	def coefficient(self, value: complex) -> None:
		self._coefficient = complex(value)


@value_equality
class PauliStringGateOperation(Operation):
	"""An operation that applies a Pauli string."""

	def __init__(self, pauli_string: PauliString, qubits: Sequence[Qid]) -> None:
		self._pauli_string = pauli_string
		self._qubits = tuple(qubits)

	@property
	def pauli_string(self) -> PauliString:
		"""The Pauli string the operation applies."""
		return self._pauli_string

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits the string acts on."""
		return self._qubits

	@property
	def gate(self) -> Gate:
		"""The gate that applies the string."""
		from ._gates import MatrixGate

		return MatrixGate(self._pauli_string.matrix())

	def _unitary_(self) -> Matrix:
		return self._pauli_string.matrix()

	def with_qubits(self, *new_qubits: Qid) -> PauliStringGateOperation:
		"""Return the same string applied to ``new_qubits``."""
		return type(self)(self._pauli_string, new_qubits)

	def _value_equality_values_(self) -> Any:
		return (self._pauli_string, self._qubits)

	def __pow__(self, exponent: int) -> Any:
		powered = self._pauli_string**exponent
		if not powered.qubits:
			# An even power of a Pauli string is the identity on its qubits.
			from ._gates import IdentityGate

			return IdentityGate(len(self._qubits)).on(*self._qubits)
		return PauliStringGateOperation(powered, self._qubits)

	def __str__(self) -> str:
		return f"{self._pauli_string}{list(self._qubits)}"

	def __repr__(self) -> str:
		return f"vqubit.PauliStringGateOperation({self._pauli_string!r}, {self._qubits!r})"


@value_equality
class SingleQubitPauliStringGateOperation(PauliStringGateOperation):
	"""A single-qubit Pauli operation."""

	def __init__(self, pauli_string: PauliString, qubit: Qid) -> None:
		if len(pauli_string.qubits) > 1:
			raise ValueError("SingleQubitPauliStringGateOperation needs at most one qubit")
		super().__init__(pauli_string, (qubit,))

	@property
	def gate(self) -> Gate:
		"""The Pauli gate for this qubit."""
		from ._gates import MatrixGate

		paulis = [self._pauli_string[q] for q in (self._qubits[0],)]
		pauli = paulis[0]
		if self._pauli_string.coefficient == 1:
			return _PAULI_GATE[pauli]
		return MatrixGate(self._pauli_string.matrix())

	def with_qubits(self, *new_qubits: Qid) -> SingleQubitPauliStringGateOperation:
		"""Return the same string applied to ``new_qubit``."""
		(qubit,) = new_qubits
		return SingleQubitPauliStringGateOperation(self._pauli_string, qubit)


@value_equality
class PauliStringPhasorGate(Gate):
	"""A gate applying ``exp(i π half_turns P)`` for a Pauli string ``P``."""

	def __init__(self, pauli_string: PauliString, *, half_turns: Any = 1.0) -> None:
		self._pauli_string = pauli_string
		self._half_turns = half_turns

	@property
	def pauli_string(self) -> PauliString:
		"""The Pauli string in the exponent."""
		return self._pauli_string

	@property
	def half_turns(self) -> Any:
		"""The angle of the phasor, in half turns."""
		return self._half_turns

	def _num_qubits_(self) -> int:
		return len(self._pauli_string.qubits)

	def _unitary_(self) -> Matrix | None:
		if self._pauli_string.coefficient != 1:
			return None
		matrix = self._pauli_string.matrix()
		dim = len(matrix)
		c = math.cos(math.pi * float(self._half_turns))
		s = math.sin(math.pi * float(self._half_turns))
		return [[(c if i == j else 0j) + 1j * s * matrix[i][j] for j in range(dim)] for i in range(dim)]

	def _value_equality_values_(self) -> Any:
		return (self._pauli_string, self._half_turns)

	def __pow__(self, power: float) -> PauliStringPhasorGate:
		return PauliStringPhasorGate(self._pauli_string, half_turns=self._half_turns * power)

	def __str__(self) -> str:
		return f"PauliStringPhasorGate({self._pauli_string}, half_turns={self._half_turns})"

	def __repr__(self) -> str:
		return f"vqubit.PauliStringPhasorGate({self._pauli_string!r}, half_turns={self._half_turns!r})"


@value_equality
class PauliStringPhasor(Operation):
	"""The operation ``exp(i π half_turns P)`` on a Pauli string's qubits."""

	def __init__(self, pauli_string: PauliString, *, half_turns: Any = 1.0) -> None:
		self._gate = PauliStringPhasorGate(pauli_string, half_turns=half_turns)

	@property
	def gate(self) -> PauliStringPhasorGate:
		"""The gate that applies the phasor."""
		return self._gate

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits the phasor acts on."""
		return self._gate.pauli_string.qubits

	@property
	def pauli_string(self) -> PauliString:
		"""The Pauli string in the exponent."""
		return self._gate.pauli_string

	@property
	def half_turns(self) -> Any:
		"""The angle of the phasor, in half turns."""
		return self._gate.half_turns

	def _unitary_(self) -> Matrix | None:
		return self._gate._unitary_()

	def with_qubits(self, *new_qubits: Qid) -> PauliStringPhasor:
		"""Return the same phasor on ``new_qubits``."""
		mapping = dict(zip(self.qubits, new_qubits, strict=True))
		remapped = PauliString({mapping[q]: p for q, p in self.pauli_string.items()}, self.pauli_string.coefficient)
		return PauliStringPhasor(remapped, half_turns=self.half_turns)

	def _value_equality_values_(self) -> Any:
		return (self.pauli_string, self.half_turns)

	def __pow__(self, power: float) -> PauliStringPhasor:
		return PauliStringPhasor(self.pauli_string, half_turns=self.half_turns * power)

	def __str__(self) -> str:
		return f"PauliStringPhasor({self.pauli_string}, half_turns={self.half_turns})"

	def __repr__(self) -> str:
		return f"vqubit.PauliStringPhasor({self.pauli_string!r}, half_turns={self.half_turns!r})"


@value_equality
class PauliMeasurementGate(Gate):
	"""Measures the qubits of a Pauli string along its Pauli axes."""

	def __init__(self, pauli_string: PauliString, key: str | None = None, invert_mask: Sequence[bool] = ()) -> None:
		self._pauli_string = pauli_string
		self._key = key
		self._invert_mask = tuple(bool(v) for v in invert_mask)

	@property
	def pauli_string(self) -> PauliString:
		"""The string whose axes are measured."""
		return self._pauli_string

	def _num_qubits_(self) -> int:
		return len(self._pauli_string.qubits)

	def _is_measurement_(self) -> bool:
		return True

	def _measurement_key_name_(self) -> str:
		if self._key is not None:
			return self._key
		return str(self._pauli_string)

	def _value_equality_values_(self) -> Any:
		return (self._pauli_string, self._key, self._invert_mask)

	def __str__(self) -> str:
		return f"PauliMeasurementGate({self._pauli_string})"

	def __repr__(self) -> str:
		return f"vqubit.PauliMeasurementGate({self._pauli_string!r}, key={self._key!r})"


def measure_single_paulistring(pauli_string: PauliString, key: str | None = None) -> GateOperation:
	"""Return a measurement of a single Pauli string."""
	return PauliMeasurementGate(pauli_string, key=key).on(*pauli_string.qubits)


def measure_paulistring_terms(pauli_sum: PauliSum, key: str | None = None) -> list[GateOperation]:
	"""Return one measurement operation per Pauli string in a sum."""
	return [measure_single_paulistring(ps, key=None if key is None else f"{key}{i}") for i, ps in enumerate(pauli_sum)]


@value_equality
class PauliSumExponential:
	"""The operator ``exp(i π exponent H)`` for a Hermitian :class:`~vqubit.PauliSum`."""

	def __init__(self, pauli_sum: Any, exponent: Any = 1.0) -> None:
		# A lone Pauli string is the one-term sum, so ``PauliString(...) * 0.5``
		# can be handed straight to the exponential.
		self._pauli_sum = pauli_sum if isinstance(pauli_sum, PauliSum) else PauliSum([(pauli_sum, 1.0)])
		self._exponent = exponent

	@property
	def pauli_sum(self) -> PauliSum:
		"""The Hermitian sum in the exponent."""
		return self._pauli_sum

	@property
	def exponent(self) -> Any:
		"""The exponent multiplying the sum, in half turns."""
		return self._exponent

	def matrix(self) -> Matrix:
		"""The unitary ``exp(iπ·exponent·H)``."""
		from ._linalg_public import map_eigenvalues

		dim = 1 << len({q for ps in self._pauli_sum for q in ps.qubits})
		total = [[0j] * dim for _ in range(dim)]
		for pauli_string, coefficient in self._pauli_sum.items():
			matrix = pauli_string.matrix()
			for i in range(dim):
				for j in range(dim):
					total[i][j] += coefficient * matrix[i][j]
		return map_eigenvalues(total, lambda v: cmath.exp(1j * math.pi * float(self._exponent) * v))

	_unitary_ = matrix

	def _value_equality_values_(self) -> Any:
		return (tuple(self._pauli_sum.items()), self._exponent)

	def __str__(self) -> str:
		return f"PauliSumExponential({self._pauli_sum}, exponent={self._exponent})"

	def __repr__(self) -> str:
		return f"vqubit.PauliSumExponential({self._pauli_sum!r}, exponent={self._exponent!r})"


@value_equality
class ProjectorString:
	"""The projector ``(I + coefficient * P) / 2`` onto a Pauli string's ``+1`` eigenspace."""

	def __init__(self, pauli_string: PauliString, coefficient: complex = 1) -> None:
		self._pauli_string = pauli_string
		self._coefficient = complex(coefficient)

	@property
	def pauli_string(self) -> PauliString:
		"""The string projected onto its ``+1`` eigenspace."""
		return self._pauli_string

	@property
	def coefficient(self) -> complex:
		"""The coefficient multiplying the projector."""
		return self._coefficient

	def matrix(self) -> Matrix:
		"""The matrix of the projector."""
		matrix = self._pauli_string.matrix()
		dim = len(matrix)
		return [[(1 + 0j if i == j else 0j) + self._coefficient * matrix[i][j] for j in range(dim)] for i in range(dim)]

	def _value_equality_values_(self) -> Any:
		return (self._pauli_string, self._coefficient)

	def __str__(self) -> str:
		return f"ProjectorString({self._pauli_string})"

	def __repr__(self) -> str:
		return f"vqubit.ProjectorString({self._pauli_string!r}, coefficient={self._coefficient!r})"


@value_equality
class ProjectorSum:
	"""A sum of :class:`ProjectorString` terms."""

	def __init__(self, terms: Mapping[PauliString, complex] | None = None) -> None:
		self._terms: dict[PauliString, complex] = dict(terms) if terms is not None else {}

	def items(self) -> Iterable[tuple[PauliString, complex]]:
		"""Iterate ``(projector, coefficient)`` pairs."""
		return iter(self._terms.items())

	def __iter__(self) -> Iterable[PauliString]:
		return iter(self._terms)

	def __len__(self) -> int:
		return len(self._terms)

	def _value_equality_values_(self) -> Any:
		return (tuple(sorted(self._terms.items())),)

	def __str__(self) -> str:
		return " + ".join(f"{c}*P({ps})" for ps, c in self._terms.items()) or "0"

	def __repr__(self) -> str:
		return f"vqubit.ProjectorSum({self._terms!r})"


@value_equality
class LinearCombinationOfGates:
	"""A linear combination of gates."""

	def __init__(self, terms: Mapping[Gate, complex]) -> None:
		self._terms = {gate: complex(c) for gate, c in terms.items()}

	@property
	def terms(self) -> dict[Gate, complex]:
		"""The ``(gate, coefficient)`` pairs of the combination."""
		return dict(self._terms)

	def matrix(self) -> Matrix:
		"""The matrix of the linear combination."""
		from ._protocols import unitary

		first = next(iter(self._terms))
		dim = 1 << first.num_qubits
		total = [[0j] * dim for _ in range(dim)]
		for gate, coefficient in self._terms.items():
			matrix = unitary(gate)
			for i in range(dim):
				for j in range(dim):
					total[i][j] += coefficient * matrix[i][j]
		return total

	def _value_equality_values_(self) -> Any:
		return (tuple(sorted(self._terms.items(), key=lambda kv: str(kv[0]))),)

	def __str__(self) -> str:
		return " + ".join(f"{c}*{g}" for g, c in self._terms.items()) or "0"

	def __repr__(self) -> str:
		return f"vqubit.LinearCombinationOfGates({self._terms!r})"


@value_equality
class LinearCombinationOfOperations:
	"""A linear combination of operations."""

	def __init__(self, terms: Mapping[Operation, complex]) -> None:
		self._terms = {op: complex(c) for op, c in terms.items()}

	@property
	def terms(self) -> dict[Operation, complex]:
		"""The ``(operation, coefficient)`` pairs of the combination."""
		return dict(self._terms)

	def matrix(self) -> Matrix:
		"""The matrix of the linear combination."""
		from ._protocols import unitary

		first = next(iter(self._terms))
		dim = 1 << len(first.qubits)
		total = [[0j] * dim for _ in range(dim)]
		for op, coefficient in self._terms.items():
			matrix = unitary(op)
			for i in range(dim):
				for j in range(dim):
					total[i][j] += coefficient * matrix[i][j]
		return total

	def _value_equality_values_(self) -> Any:
		return (tuple(sorted(self._terms.items(), key=lambda kv: str(kv[0]))),)

	def __str__(self) -> str:
		return " + ".join(f"{c}*{op}" for op, c in self._terms.items()) or "0"

	def __repr__(self) -> str:
		return f"vqubit.LinearCombinationOfOperations({self._terms!r})"


@value_equality
class PauliInteractionGate(Gate):
	"""Applies a pair of Paulis to two qubits, optionally inverted."""

	def __init__(self, pauli0: Pauli, invert0: bool, pauli1: Pauli, invert1: bool) -> None:
		self._pauli0 = pauli0
		self._invert0 = bool(invert0)
		self._pauli1 = pauli1
		self._invert1 = bool(invert1)

	@property
	def pauli0(self) -> Pauli:
		"""The Pauli applied on the first qubit."""
		return self._pauli0

	@property
	def pauli1(self) -> Pauli:
		"""The Pauli applied on the second qubit."""
		return self._pauli1

	def _num_qubits_(self) -> int:
		return 2

	def _unitary_(self) -> Matrix:
		from ._linalg import kron

		m0 = self._pauli0.matrix()
		m1 = self._pauli1.matrix()
		if self._invert0:
			m0 = [[-z for z in row] for row in m0]
		if self._invert1:
			m1 = [[-z for z in row] for row in m1]
		return kron(m0, m1)

	def _value_equality_values_(self) -> Any:
		return (self._pauli0, self._invert0, self._pauli1, self._invert1)

	def __str__(self) -> str:
		return f"PauliInteractionGate({self._pauli0}{'^-1' if self._invert0 else ''}, {self._pauli1}{'^-1' if self._invert1 else ''})"

	def __repr__(self) -> str:
		return f"vqubit.PauliInteractionGate({self._pauli0!r}, {self._invert0!r}, {self._pauli1!r}, {self._invert1!r})"


__all__ = [
	"BaseDensePauliString",
	"LinearCombinationOfGates",
	"LinearCombinationOfOperations",
	"MutableDensePauliString",
	"MutablePauliString",
	"PauliInteractionGate",
	"PauliMeasurementGate",
	"PauliStringGateOperation",
	"PauliStringPhasor",
	"PauliStringPhasorGate",
	"PauliSumExponential",
	"ProjectorString",
	"ProjectorSum",
	"SingleQubitPauliStringGateOperation",
	"measure_paulistring_terms",
	"measure_single_paulistring",
]
