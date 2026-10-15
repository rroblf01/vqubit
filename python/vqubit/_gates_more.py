"""Additional gates, matching the rest of Cirq's ``cirq.ops`` surface.

This module adds the controlled power gates (``CXPowGate`` and friends), the
swap-family power gates, the diagonal gates, the permutation and Fourier gates,
and the lowercase factory helpers (``q``, ``rx``, ``ms`` ...). The base classes
and the elementary gates live in :mod:`vqubit._gates`.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Sequence
from typing import Any

from ._gates import (
	EigenGate,
	Gate,
	GlobalPhaseGate,
	MatrixGate,
	MSGate,
	WaitGate,
	XPowGate,
	YPowGate,
	ZPowGate,
	identity,
	is_symbolic,
	resolve_value,
	value_equality,
)
from ._linalg import Matrix
from ._parameter import ParamResolver
from ._qubits import Qid


class _PowGate(Gate):
	"""Shared machinery for the multi-qubit power gates."""

	def __init__(self, *, exponent: Any = 1.0, global_shift: Any = 0.0) -> None:
		self._exponent = exponent if is_symbolic(exponent) else float(exponent)
		self._global_shift = global_shift if is_symbolic(global_shift) else float(global_shift)

	@property
	def exponent(self) -> Any:
		return self._exponent

	@property
	def global_shift(self) -> Any:
		return self._global_shift

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._exponent) or is_symbolic(self._global_shift)

	def _value_equality_values_(self) -> Any:
		return (self._exponent, self._global_shift)

	def __pow__(self, power: float) -> _PowGate:
		return type(self)(exponent=self._exponent * power, global_shift=self._global_shift)

	def _inverse_(self) -> _PowGate:
		return type(self)(exponent=-self._exponent, global_shift=self._global_shift)

	def _resolve_parameters_(self, resolver: ParamResolver) -> _PowGate:
		if not self._is_parameterized_():
			return self
		return type(self)(
			exponent=resolve_value(self._exponent, resolver),
			global_shift=resolve_value(self._global_shift, resolver),
		)

	def __str__(self) -> str:
		return f"{type(self).__name__}(exponent={self._exponent})"

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}(exponent={self._exponent!r}, global_shift={self._global_shift!r})"


@value_equality
class _ControlledPowGate(_PowGate):
	"""A power of a controlled single-qubit power gate."""

	_SubGate: type[EigenGate] = XPowGate
	_NumControls = 1
	_Symbol = "X"

	def _num_qubits_(self) -> int:
		return self._NumControls + 1

	def _unitary_(self) -> Matrix | None:
		sub = self._SubGate(exponent=self._exponent, global_shift=self._global_shift)._unitary_()
		if sub is None:
			return None
		block = len(sub)
		dim = (1 << self._NumControls) * block
		matrix = identity(dim)
		offset = dim - block
		for i in range(block):
			for j in range(block):
				matrix[offset + i][offset + j] = sub[i][j]
		return matrix

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("@",) * self._NumControls + (self._Symbol,)


class CXPowGate(_ControlledPowGate):
	"""A power of the controlled-``X`` gate."""

	_SubGate = XPowGate
	_Symbol = "X"


class CYPowGate(_ControlledPowGate):
	"""A power of the controlled-``Y`` gate."""

	_SubGate = YPowGate
	_Symbol = "Y"


class CZPowGate(_ControlledPowGate):
	"""A power of the controlled-``Z`` gate."""

	_SubGate = ZPowGate
	_Symbol = "Z"


class CCXPowGate(_ControlledPowGate):
	"""A power of the doubly-controlled ``X`` gate."""

	_SubGate = XPowGate
	_NumControls = 2
	_Symbol = "X"


class CCYPowGate(_ControlledPowGate):
	"""A power of the doubly-controlled ``Y`` gate."""

	_SubGate = YPowGate
	_NumControls = 2
	_Symbol = "Y"


class CCZPowGate(_ControlledPowGate):
	"""A power of the doubly-controlled ``Z`` gate."""

	_SubGate = ZPowGate
	_NumControls = 2
	_Symbol = "Z"


CNotPowGate = CXPowGate
CCNotPowGate = CCXPowGate


@value_equality
class SwapPowGate(_PowGate):
	"""A power of the ``SWAP`` gate."""

	def _num_qubits_(self) -> int:
		return 2

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		phase = cmath.exp(1j * math.pi * (self._exponent - self._global_shift))
		a = (1 + phase) / 2
		b = (1 - phase) / 2
		return [
			[a + b, 0j, 0j, 0j],
			[0j, a, b, 0j],
			[0j, b, a, 0j],
			[0j, 0j, 0j, a + b],
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Swap", "Swap")


@value_equality
class ISwapPowGate(_PowGate):
	"""A power of the ``iSWAP`` gate."""

	def _num_qubits_(self) -> int:
		return 2

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		# The {|01>, |10>} block has eigenvalues +/- i; a global shift scales
		# each eigencomponent by exp(i pi g * eigenvalue).
		plus = cmath.exp(1j * math.pi * self._exponent / 2) * cmath.exp(1j * math.pi * self._global_shift * 1j)
		minus = cmath.exp(-1j * math.pi * self._exponent / 2) * cmath.exp(-1j * math.pi * self._global_shift * 1j)
		a = (plus + minus) / 2
		b = (plus - minus) / 2
		return [
			[1 + 0j, 0j, 0j, 0j],
			[0j, a, b, 0j],
			[0j, b, a, 0j],
			[0j, 0j, 0j, 1 + 0j],
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("iSwap", "iSwap")


@value_equality
class PhasedISwapPowGate(_PowGate):
	"""A fractional iSWAP conjugated by ``Z`` rotations.

	With phase exponent ``p`` and exponent ``t`` it is
	``(Z^-p ⊗ Z^p) iSWAP^t (Z^p ⊗ Z^-p)``.
	"""

	def __init__(self, *, phase_exponent: Any = 0.25, exponent: Any = 1.0, global_shift: Any = 0.0) -> None:
		super().__init__(exponent=exponent, global_shift=global_shift)
		self._phase_exponent = phase_exponent if is_symbolic(phase_exponent) else float(phase_exponent)

	@property
	def phase_exponent(self) -> Any:
		return self._phase_exponent

	def _num_qubits_(self) -> int:
		return 2

	def _is_parameterized_(self) -> bool:
		return super()._is_parameterized_() or is_symbolic(self._phase_exponent)

	def _value_equality_values_(self) -> Any:
		return (self._phase_exponent, self._exponent, self._global_shift)

	def __pow__(self, power: float) -> PhasedISwapPowGate:
		return PhasedISwapPowGate(
			phase_exponent=self._phase_exponent,
			exponent=self._exponent * power,
			global_shift=self._global_shift,
		)

	def _inverse_(self) -> PhasedISwapPowGate:
		return self.__pow__(-1)

	def _resolve_parameters_(self, resolver: ParamResolver) -> PhasedISwapPowGate:
		if not self._is_parameterized_():
			return self
		return PhasedISwapPowGate(
			phase_exponent=resolve_value(self._phase_exponent, resolver),
			exponent=resolve_value(self._exponent, resolver),
			global_shift=resolve_value(self._global_shift, resolver),
		)

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		c = math.cos(math.pi * self._exponent / 2)
		s = math.sin(math.pi * self._exponent / 2)
		f = cmath.exp(2j * math.pi * self._phase_exponent)
		phase = cmath.exp(1j * math.pi * self._exponent * self._global_shift)
		return [
			[phase, 0j, 0j, 0j],
			[0j, phase * c, phase * 1j * s * f.conjugate(), 0j],
			[0j, phase * 1j * s * f, phase * c, 0j],
			[0j, 0j, 0j, phase],
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("PhISwap", "PhISwap")

	def __str__(self) -> str:
		return f"PhasedISwapPowGate(phase_exponent={self._phase_exponent}, exponent={self._exponent})"

	def __repr__(self) -> str:
		return (
			f"vqubit.PhasedISwapPowGate(phase_exponent={self._phase_exponent!r}, "
			f"exponent={self._exponent!r}, global_shift={self._global_shift!r})"
		)


def givens(angle_rads: Any) -> PhasedISwapPowGate:
	"""Return ``exp(-i angle (Y⊗X - X⊗Y) / 2)``."""
	return PhasedISwapPowGate(exponent=2 * angle_rads / math.pi)


@value_equality
class DiagonalGate(Gate):
	"""A gate whose unitary is a given diagonal of phases."""

	def __init__(self, diag: Sequence[complex]) -> None:
		entries = [complex(z) for z in diag]
		dim = len(entries)
		if dim == 0 or dim & (dim - 1):
			raise ValueError("DiagonalGate needs a power-of-two number of entries")
		self._diag = entries

	def _num_qubits_(self) -> int:
		return len(self._diag).bit_length() - 1

	def _unitary_(self) -> Matrix:
		n = len(self._diag)
		return [[self._diag[i] if i == j else 0j for j in range(n)] for i in range(n)]

	def _value_equality_values_(self) -> Any:
		return (tuple(self._diag),)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("diag",) * self.num_qubits

	def __str__(self) -> str:
		return f"DiagonalGate({self._diag})"

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}({self._diag!r})"


class TwoQubitDiagonalGate(DiagonalGate):
	"""A two-qubit diagonal gate."""

	def __init__(self, diag: Sequence[complex]) -> None:
		super().__init__(diag)
		if len(diag) != 4:
			raise ValueError("TwoQubitDiagonalGate needs four entries")


class ThreeQubitDiagonalGate(DiagonalGate):
	"""A three-qubit diagonal gate."""

	def __init__(self, diag: Sequence[complex]) -> None:
		super().__init__(diag)
		if len(diag) != 8:
			raise ValueError("ThreeQubitDiagonalGate needs eight entries")


@value_equality
class QubitPermutationGate(Gate):
	"""A gate that permutes the basis states of its qubits."""

	def __init__(self, permutation: Sequence[int]) -> None:
		self._permutation = tuple(int(p) for p in permutation)
		dim = len(self._permutation)
		if dim == 0 or dim & (dim - 1) or sorted(self._permutation) != list(range(dim)):
			raise ValueError("QubitPermutationGate needs a permutation of a power-of-two range")

	def _num_qubits_(self) -> int:
		return len(self._permutation).bit_length() - 1

	def _unitary_(self) -> Matrix:
		dim = len(self._permutation)
		return [[1 + 0j if self._permutation[j] == i else 0j for j in range(dim)] for i in range(dim)]

	def _value_equality_values_(self) -> Any:
		return (self._permutation,)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("perm",) * self.num_qubits

	def __str__(self) -> str:
		return f"QubitPermutationGate({self._permutation})"

	def __repr__(self) -> str:
		return f"vqubit.QubitPermutationGate({self._permutation!r})"


@value_equality
class QuantumFourierTransformGate(Gate):
	"""The quantum Fourier transform on ``num_qubits`` qubits."""

	def __init__(self, num_qubits: int, *, without_reverse: bool = False) -> None:
		self._num_qubits = int(num_qubits)
		self._without_reverse = bool(without_reverse)

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _unitary_(self) -> Matrix:
		dim = 1 << self._num_qubits
		root = cmath.exp(2j * math.pi / dim)
		matrix = [[root ** (i * j) / math.sqrt(dim) for j in range(dim)] for i in range(dim)]
		if not self._without_reverse:
			matrix = [matrix[dim - 1 - i] for i in range(dim)]
		return matrix

	def _value_equality_values_(self) -> Any:
		return (self._num_qubits, self._without_reverse)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("QFT",) * self._num_qubits

	def __str__(self) -> str:
		return f"QuantumFourierTransformGate({self._num_qubits})"

	def __repr__(self) -> str:
		return f"vqubit.QuantumFourierTransformGate({self._num_qubits}, without_reverse={self._without_reverse!r})"


# Fixed-gate instances and aliases, as in Cirq.
CY: CYPowGate = CYPowGate()
CCZ: CCZPowGate = CCZPowGate()
CCY: CCYPowGate = CCYPowGate()
CCNOT: CCXPowGate = CCXPowGate()
ISWAP_INV: ISwapPowGate = ISwapPowGate(exponent=-1)
SQRT_ISWAP_INV: ISwapPowGate = ISwapPowGate(exponent=-0.5)
SWAP_INV: SwapPowGate = SwapPowGate(exponent=-1)


def rx(radians: Any) -> Gate:
	"""Return an ``X`` rotation by ``radians``."""
	from ._gates import Rx

	return Rx(radians)


def ry(radians: Any) -> Gate:
	"""Return a ``Y`` rotation by ``radians``."""
	from ._gates import Ry

	return Ry(radians)


def rz(radians: Any) -> Gate:
	"""Return a ``Z`` rotation by ``radians``."""
	from ._gates import Rz

	return Rz(radians)


def cphase(radians: Any) -> CZPowGate:
	"""Return a controlled-phase gate of angle ``radians``."""
	return CZPowGate(exponent=radians / math.pi)


def ms(theta: Any) -> MSGate:
	"""Return a Mølmer-Sørensen gate of angle ``theta``."""
	return MSGate(2, theta)


def wait(duration_nanos: float, num_qubits: int = 1) -> WaitGate:
	"""Return a wait gate of the given duration."""
	return WaitGate(duration_nanos, num_qubits)


def global_phase_operation(coefficient: complex) -> Any:
	"""Return an operation implementing a global phase of ``coefficient``."""
	return GlobalPhaseGate(cmath.phase(coefficient) / math.pi).on()


def reset_each(*qubits: Qid) -> tuple[Any, ...]:
	"""Return a reset operation on each of the given qubits."""
	from ._channels import reset

	return tuple(reset(qubit) for qubit in qubits)


def identity_each(*qubits: Qid) -> tuple[Any, ...]:
	"""Return an identity operation on each of the given qubits."""
	from ._gates import I

	return tuple(I.on(qubit) for qubit in qubits)


def qft(*qubits: Qid, without_reverse: bool = False, inverse: bool = False) -> list[Any]:
	"""Return operations implementing the quantum Fourier transform on ``qubits``."""
	gate: Gate = QuantumFourierTransformGate(len(qubits), without_reverse=without_reverse)
	if inverse:
		from ._linalg import dagger

		gate = MatrixGate(dagger(gate._unitary_()))
	return [gate.on(*qubits)]


def q(index: int) -> Qid:
	"""Return the line qubit at ``index`` (Cirq's short ``q`` helper)."""
	from ._qubits import LineQubit

	return LineQubit(index)


@value_equality
class ParallelGate(Gate):
	"""Applies a sub gate to two disjoint groups of qubits in parallel."""

	def __init__(self, sub_gate: Gate) -> None:
		if sub_gate.num_qubits == 0:
			raise ValueError("ParallelGate needs a sub gate on at least one qubit")
		self._sub_gate = sub_gate

	@property
	def sub_gate(self) -> Gate:
		return self._sub_gate

	def _num_qubits_(self) -> int:
		return 2 * self._sub_gate.num_qubits

	def _unitary_(self) -> Matrix | None:
		sub = self._sub_gate._unitary_()
		if sub is None:
			return None
		from ._linalg import kron

		return kron(sub, sub)

	def _value_equality_values_(self) -> Any:
		return (self._sub_gate,)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		half = self._sub_gate.num_qubits
		return self._sub_gate._diagram_symbols_(qubits[:half]) + self._sub_gate._diagram_symbols_(qubits[half:])

	def __str__(self) -> str:
		return f"ParallelGate({self._sub_gate})"

	def __repr__(self) -> str:
		return f"vqubit.ParallelGate({self._sub_gate!r})"


def parallel_gate_op(sub_gate: Gate, *qubits: Qid) -> Any:
	"""Return a parallel application of ``sub_gate`` to disjoint qubit groups."""
	return ParallelGate(sub_gate).on(*qubits)


__all__ = [
	"CCNOT",
	"CCY",
	"CCZ",
	"CY",
	"ISWAP_INV",
	"SQRT_ISWAP_INV",
	"SWAP_INV",
	"CCNotPowGate",
	"CCXPowGate",
	"CCYPowGate",
	"CCZPowGate",
	"CNotPowGate",
	"CXPowGate",
	"CYPowGate",
	"CZPowGate",
	"DiagonalGate",
	"ISwapPowGate",
	"PhasedISwapPowGate",
	"QuantumFourierTransformGate",
	"QubitPermutationGate",
	"SwapPowGate",
	"ThreeQubitDiagonalGate",
	"TwoQubitDiagonalGate",
	"cphase",
	"givens",
	"global_phase_operation",
	"identity_each",
	"ms",
	"q",
	"qft",
	"reset_each",
	"rx",
	"ry",
	"rz",
	"wait",
]
