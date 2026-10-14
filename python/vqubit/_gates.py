"""Quantum gates.

A :class:`Gate` is a fixed quantum operation on a number of qubits. The gates
here are the ones you reach for first: the Pauli and phase gates, the
Hadamard, the rotation gates and the standard two-qubit gates.

Most single-qubit gates are *power* gates (:class:`EigenGate`): ``X`` is
``XPowGate(exponent=1)`` and ``X**0.5`` is the square root of ``X``. Each
power gate is defined by its eigencomponents, so powers, inverses and unitaries
all follow from one definition.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Iterable, Sequence
from typing import Any

from ._linalg import Matrix
from ._ops import GateOperation
from ._qubits import Qid
from ._value import value_equality

_EPS = 1e-12


class Gate:
	"""The base class of all gates."""

	def _num_qubits_(self) -> int:
		raise NotImplementedError

	@property
	def num_qubits(self) -> int:
		return self._num_qubits_()

	def _qid_shape_(self) -> tuple[int, ...]:
		return (2,) * self._num_qubits_()

	@property
	def qid_shape(self) -> tuple[int, ...]:
		return self._qid_shape_()

	def _unitary_(self) -> Matrix | None:
		return None

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return (str(self),) * self.num_qubits

	def _circuit_diagram_info_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return self._diagram_symbols_(qubits)

	def on(self, *qubits: Qid | Iterable[Qid]) -> GateOperation:
		"""Return this gate applied to ``qubits``."""
		flat: list[Qid] = []
		for qubit in qubits:
			if isinstance(qubit, Qid):
				flat.append(qubit)
			else:
				flat.extend(qubit)
		return GateOperation(self, flat)

	def on_each(self, *targets: Qid | Iterable[Qid]) -> tuple[GateOperation, ...]:
		"""Return this gate applied to each consecutive group of qubits."""
		operations: list[GateOperation] = []
		for target in targets:
			if isinstance(target, Qid):
				target = [target]
			target = list(target)
			if len(target) % self.num_qubits != 0:
				raise ValueError(f"{len(target)} qubits do not divide into groups of {self.num_qubits}")
			for k in range(0, len(target), self.num_qubits):
				operations.append(self.on(*target[k : k + self.num_qubits]))
		return tuple(operations)

	def __pow__(self, power: float) -> Gate:
		raise NotImplementedError(f"{type(self).__name__} does not support powers")

	def __call__(self, *qubits: Qid | Iterable[Qid]) -> GateOperation:
		return self.on(*qubits)

	def __str__(self) -> str:
		return type(self).__name__


@value_equality
class EigenGate(Gate):
	"""A single-qubit gate defined by its eigencomponents and an exponent."""

	def __init__(self, *, exponent: float = 1.0, global_shift: float = 0.0) -> None:
		self._exponent = float(exponent)
		self._global_shift = float(global_shift)

	@property
	def exponent(self) -> float:
		return self._exponent

	@property
	def global_shift(self) -> float:
		return self._global_shift

	def _eigen_components(self) -> list[tuple[float, Matrix]]:
		raise NotImplementedError

	def _unitary_(self) -> Matrix:
		result: Matrix = [[0j, 0j], [0j, 0j]]
		for eigenvalue, projector in self._eigen_components():
			if eigenvalue == 0:
				coefficient = 1.0 if self._exponent == 0 else 0.0
			else:
				coefficient = complex(eigenvalue) ** self._exponent
			coefficient *= cmath.exp(1j * math.pi * self._global_shift * eigenvalue)
			for i in range(2):
				for j in range(2):
					result[i][j] += coefficient * projector[i][j]
		return result

	def _value_equality_values_(self) -> Any:
		return (self._exponent, self._global_shift)

	def __pow__(self, power: float) -> EigenGate:
		return type(self)(exponent=self._exponent * power, global_shift=self._global_shift)

	def _inverse_(self) -> EigenGate:
		return type(self)(exponent=-self._exponent, global_shift=self._global_shift)

	def __str__(self) -> str:
		return f"{type(self).__name__}(exponent={_clean(self._exponent)})"


def _clean(value: float) -> float | int:
	"""Return ``value`` as an int when it is a whole number, else as a float."""
	return int(value) if abs(value - round(value)) < _EPS else float(value)


class XPowGate(EigenGate):
	"""The power of the Pauli ``X`` gate."""

	def _num_qubits_(self) -> int:
		return 1

	def _eigen_components(self) -> list[tuple[float, Matrix]]:
		return [
			(1.0, [[0.5, 0.5], [0.5, 0.5]]),
			(-1.0, [[0.5, -0.5], [-0.5, 0.5]]),
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("X",)


class YPowGate(EigenGate):
	"""The power of the Pauli ``Y`` gate."""

	def _num_qubits_(self) -> int:
		return 1

	def _eigen_components(self) -> list[tuple[float, Matrix]]:
		return [
			(1.0, [[0.5, -0.5j], [0.5j, 0.5]]),
			(-1.0, [[0.5, 0.5j], [-0.5j, 0.5]]),
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Y",)


class ZPowGate(EigenGate):
	"""The power of the Pauli ``Z`` gate."""

	def _num_qubits_(self) -> int:
		return 1

	def _eigen_components(self) -> list[tuple[float, Matrix]]:
		return [
			(1.0, [[1.0, 0.0], [0.0, 0.0]]),
			(-1.0, [[0.0, 0.0], [0.0, 1.0]]),
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		if abs(self._exponent - 0.5) < _EPS:
			return ("S",)
		if abs(self._exponent - 0.25) < _EPS:
			return ("T",)
		return ("Z",)


class HPowGate(EigenGate):
	"""The power of the Hadamard gate."""

	def _num_qubits_(self) -> int:
		return 1

	def _eigen_components(self) -> list[tuple[float, Matrix]]:
		s = math.sqrt(0.5)
		return [
			(1.0, [[(1 + s) / 2, s / 2], [s / 2, (1 - s) / 2]]),
			(-1.0, [[(1 - s) / 2, -s / 2], [-s / 2, (1 + s) / 2]]),
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("H",)


@value_equality
class RxGate(Gate):
	"""A rotation of ``radians`` around the X axis."""

	def __init__(self, radians: float) -> None:
		self._radians = float(radians)

	@property
	def radians(self) -> float:
		return self._radians

	def _num_qubits_(self) -> int:
		return 1

	def _unitary_(self) -> Matrix:
		c = math.cos(self._radians / 2)
		s = math.sin(self._radians / 2)
		return [[c + 0j, -1j * s], [-1j * s, c + 0j]]

	def _value_equality_values_(self) -> Any:
		return (self._radians,)

	def _inverse_(self) -> RxGate:
		return RxGate(-self._radians)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Rx",)

	def __str__(self) -> str:
		return f"Rx({_clean(self._radians)})"


@value_equality
class RyGate(Gate):
	"""A rotation of ``radians`` around the Y axis."""

	def __init__(self, radians: float) -> None:
		self._radians = float(radians)

	@property
	def radians(self) -> float:
		return self._radians

	def _num_qubits_(self) -> int:
		return 1

	def _unitary_(self) -> Matrix:
		c = math.cos(self._radians / 2)
		s = math.sin(self._radians / 2)
		return [[c + 0j, -s + 0j], [s + 0j, c + 0j]]

	def _value_equality_values_(self) -> Any:
		return (self._radians,)

	def _inverse_(self) -> RyGate:
		return RyGate(-self._radians)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Ry",)

	def __str__(self) -> str:
		return f"Ry({_clean(self._radians)})"


@value_equality
class RzGate(Gate):
	"""A rotation of ``radians`` around the Z axis."""

	def __init__(self, radians: float) -> None:
		self._radians = float(radians)

	@property
	def radians(self) -> float:
		return self._radians

	def _num_qubits_(self) -> int:
		return 1

	def _unitary_(self) -> Matrix:
		return [[cmath.exp(-0.5j * self._radians), 0j], [0j, cmath.exp(0.5j * self._radians)]]

	def _value_equality_values_(self) -> Any:
		return (self._radians,)

	def _inverse_(self) -> RzGate:
		return RzGate(-self._radians)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Rz",)

	def __str__(self) -> str:
		return f"Rz({_clean(self._radians)})"


@value_equality
class IdentityGate(Gate):
	"""The identity on ``num_qubits`` qubits."""

	def __init__(self, num_qubits: int = 1) -> None:
		self._num_qubits = int(num_qubits)

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _unitary_(self) -> Matrix | None:
		if self._num_qubits == 1:
			return [[1.0 + 0j if i == j else 0j for j in range(2)] for i in range(2)]
		return None

	def _value_equality_values_(self) -> Any:
		return (self._num_qubits,)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("I",) * self._num_qubits


@value_equality
class MeasurementGate(Gate):
	"""A measurement of one or more qubits into a classical key."""

	def __init__(self, key: str, num_qubits: int = 1, invert_mask: Sequence[bool] = ()) -> None:
		if not isinstance(key, str):
			raise TypeError("MeasurementGate key must be a str")
		invert_mask = tuple(bool(b) for b in invert_mask)
		if len(invert_mask) > num_qubits:
			raise ValueError("invert_mask is longer than the number of qubits")
		self._key = key
		self._num_qubits = int(num_qubits)
		self._invert_mask = invert_mask

	@property
	def key(self) -> str:
		return self._key

	@property
	def invert_mask(self) -> tuple[bool, ...]:
		return self._invert_mask

	@property
	def full_invert_mask(self) -> tuple[bool, ...]:
		return self._invert_mask + (False,) * (self._num_qubits - len(self._invert_mask))

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _value_equality_values_(self) -> Any:
		return (self._key, self._num_qubits, self._invert_mask)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("M",) * self._num_qubits

	def __str__(self) -> str:
		return f"measure({self._key!r})"


@value_equality
class _FixedGate(Gate):
	"""Shared machinery for the fixed two-qubit gates.

	The matrix and the diagram symbols are produced by static methods rather
	than stored as class attributes, so a gate instance never hands out a
	shared, mutable matrix.
	"""

	@staticmethod
	def _base_matrix() -> Matrix:
		raise NotImplementedError

	@staticmethod
	def _base_symbols() -> tuple[str, ...]:
		raise NotImplementedError

	def _num_qubits_(self) -> int:
		return 2

	def _unitary_(self) -> Matrix:
		return [list(row) for row in self._base_matrix()]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return self._base_symbols()

	def _value_equality_values_(self) -> Any:
		return ()

	def _inverse_(self) -> _FixedGate:
		return type(self)()

	def __str__(self) -> str:
		return type(self).__name__


class CNOTGate(_FixedGate):
	"""The controlled-NOT gate: ``X`` on the target when the control is ``1``."""

	@staticmethod
	def _base_symbols() -> tuple[str, ...]:
		return ("@", "X")

	@staticmethod
	def _base_matrix() -> Matrix:
		return [
			[1.0 + 0j, 0j, 0j, 0j],
			[0j, 1.0 + 0j, 0j, 0j],
			[0j, 0j, 0j, 1.0 + 0j],
			[0j, 0j, 1.0 + 0j, 0j],
		]


class CZGate(_FixedGate):
	"""The controlled-Z gate: ``Z`` on the target when the control is ``1``."""

	@staticmethod
	def _base_symbols() -> tuple[str, ...]:
		return ("@", "@")

	@staticmethod
	def _base_matrix() -> Matrix:
		return [
			[1.0 + 0j, 0j, 0j, 0j],
			[0j, 1.0 + 0j, 0j, 0j],
			[0j, 0j, 1.0 + 0j, 0j],
			[0j, 0j, 0j, -1.0 + 0j],
		]


class SWAPGate(_FixedGate):
	"""The gate that exchanges two qubits."""

	@staticmethod
	def _base_symbols() -> tuple[str, ...]:
		return ("×", "×")

	@staticmethod
	def _base_matrix() -> Matrix:
		return [
			[1.0 + 0j, 0j, 0j, 0j],
			[0j, 0j, 1.0 + 0j, 0j],
			[0j, 1.0 + 0j, 0j, 0j],
			[0j, 0j, 0j, 1.0 + 0j],
		]


# Instances, as in Cirq.
X: XPowGate = XPowGate()
Y: YPowGate = YPowGate()
Z: ZPowGate = ZPowGate()
H: HPowGate = HPowGate()
S: ZPowGate = ZPowGate(exponent=0.5)
T: ZPowGate = ZPowGate(exponent=0.25)
I: IdentityGate = IdentityGate()  # noqa: E741 - the name matches Cirq's identity gate
CNOT: CNOTGate = CNOTGate()
CX: CNOTGate = CNOTGate()
CZ: CZGate = CZGate()
SWAP: SWAPGate = SWAPGate()


def Rx(radians: float) -> RxGate:
	"""Return a rotation of ``radians`` around the X axis."""
	return RxGate(radians)


def Ry(radians: float) -> RyGate:
	"""Return a rotation of ``radians`` around the Y axis."""
	return RyGate(radians)


def Rz(radians: float) -> RzGate:
	"""Return a rotation of ``radians`` around the Z axis."""
	return RzGate(radians)


def measure(*qubits: Qid | Iterable[Qid], key: str | None = None) -> GateOperation:
	"""Return a measurement of ``qubits`` into the classical key ``key``."""
	flat: list[Qid] = []
	for q in qubits:
		if isinstance(q, Qid):
			flat.append(q)
		else:
			flat.extend(q)
	if not flat:
		raise ValueError("measure() needs at least one qubit")
	if key is None:
		key = "".join(str(q) for q in flat)
	return MeasurementGate(key, num_qubits=len(flat)).on(*flat)


def measure_each(*qubits: Qid, key_func: Any = None) -> tuple[GateOperation, ...]:
	"""Return one single-qubit measurement per qubit."""
	if key_func is None:
		key_func = str
	return tuple(MeasurementGate(key_func(q), num_qubits=1).on(q) for q in qubits)
