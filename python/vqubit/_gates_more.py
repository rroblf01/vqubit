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
		phase = cmath.exp(1j * math.pi * self._exponent)
		a = (1 + phase) / 2
		b = (1 - phase) / 2
		shift = cmath.exp(1j * math.pi * self._exponent * self._global_shift)
		return [
			[shift * (a + b), 0j, 0j, 0j],
			[0j, shift * a, shift * b, 0j],
			[0j, shift * b, shift * a, 0j],
			[0j, 0j, 0j, shift * (a + b)],
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
		c = math.cos(math.pi * self._exponent / 2)
		s = math.sin(math.pi * self._exponent / 2)
		shift = cmath.exp(1j * math.pi * self._exponent * self._global_shift)
		return [
			[shift, 0j, 0j, 0j],
			[0j, shift * c, shift * 1j * s, 0j],
			[0j, shift * 1j * s, shift * c, 0j],
			[0j, 0j, 0j, shift],
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


def riswap(rads: Any) -> ISwapPowGate:
	"""Return ``exp(+i rads (X⊗X + Y⊗Y) / 2)``."""
	return ISwapPowGate(exponent=2 * rads / math.pi)


def _canonicalize_angle(value: Any) -> Any:
	"""Shift a 2π-periodic angle into ``[-π, π)``."""
	if is_symbolic(value):
		return value
	period = 2 * math.pi
	return float(value) - period * math.floor((float(value) + math.pi) / period)


@value_equality
class PhaseGradientGate(Gate):
	"""Phases each basis state ``|x>`` by ``omega**x`` with ``omega = e^{2πi/2^n}``."""

	def __init__(self, *, num_qubits: int, exponent: Any) -> None:
		self._num_qubits = int(num_qubits)
		self._exponent = exponent if is_symbolic(exponent) else float(exponent)

	@property
	def exponent(self) -> Any:
		return self._exponent

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _unitary_(self) -> Matrix | None:
		if is_symbolic(self._exponent):
			return None
		size = 1 << self._num_qubits
		return [
			[cmath.exp(1j * math.pi * 2 * i * self._exponent / size) if i == j else 0j for j in range(size)]
			for i in range(size)
		]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._exponent)

	def _resolve_parameters_(self, resolver: ParamResolver) -> PhaseGradientGate:
		if not self._is_parameterized_():
			return self
		return PhaseGradientGate(num_qubits=self._num_qubits, exponent=resolve_value(self._exponent, resolver))

	def __pow__(self, power: float) -> PhaseGradientGate:
		return PhaseGradientGate(num_qubits=self._num_qubits, exponent=self._exponent * power)

	def _value_equality_values_(self) -> Any:
		return (self._num_qubits, self._exponent)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Grad", *(f"#{k + 1}" for k in range(1, self._num_qubits)))

	def __str__(self) -> str:
		return f"PhaseGradientGate(num_qubits={self._num_qubits}, exponent={self._exponent})"

	def __repr__(self) -> str:
		return f"vqubit.PhaseGradientGate(num_qubits={self._num_qubits}, exponent={self._exponent!r})"


@value_equality
class PhasedFSimGate(Gate):
	"""The general excitation-preserving two-qubit gate ``PhasedFSim(theta, zeta, chi, gamma, phi)``."""

	def __init__(self, theta: Any, zeta: Any = 0.0, chi: Any = 0.0, gamma: Any = 0.0, phi: Any = 0.0) -> None:
		self._theta = _canonicalize_angle(theta)
		self._zeta = _canonicalize_angle(zeta)
		self._chi = _canonicalize_angle(chi)
		self._gamma = _canonicalize_angle(gamma)
		self._phi = _canonicalize_angle(phi)

	@property
	def theta(self) -> Any:
		return self._theta

	@property
	def zeta(self) -> Any:
		return self._zeta

	@property
	def chi(self) -> Any:
		return self._chi

	@property
	def gamma(self) -> Any:
		return self._gamma

	@property
	def phi(self) -> Any:
		return self._phi

	@staticmethod
	def from_fsim_rz(theta: Any, phi: Any, rz_angles_before: Any, rz_angles_after: Any) -> PhasedFSimGate:
		"""Create a :class:`PhasedFSimGate` from FSim plus Rz angles."""
		b0, b1 = rz_angles_before
		a0, a1 = rz_angles_after
		gamma = (-b0 - b1 - a0 - a1) / 2.0
		zeta = (b0 - b1 + a0 - a1) / 2.0
		chi = (b0 - b1 - a0 + a1) / 2.0
		return PhasedFSimGate(theta, zeta, chi, gamma, phi)

	@property
	def rz_angles_before(self) -> tuple[Any, Any]:
		return ((-self.gamma + self.zeta + self.chi) / 2.0, (-self.gamma - self.zeta - self.chi) / 2.0)

	@property
	def rz_angles_after(self) -> tuple[Any, Any]:
		return ((-self.gamma + self.zeta - self.chi) / 2.0, (-self.gamma - self.zeta + self.chi) / 2.0)

	def _num_qubits_(self) -> int:
		return 2

	def _is_parameterized_(self) -> bool:
		return any(is_symbolic(a) for a in (self._theta, self._zeta, self._chi, self._gamma, self._phi))

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		a = math.cos(self._theta)
		b = -1j * math.sin(self._theta)
		c = cmath.exp(-1j * self._phi)
		f1 = cmath.exp(-1j * self._gamma - 1j * self._zeta)
		f2 = cmath.exp(-1j * self._gamma + 1j * self._chi)
		f3 = cmath.exp(-1j * self._gamma - 1j * self._chi)
		f4 = cmath.exp(-1j * self._gamma + 1j * self._zeta)
		f5 = cmath.exp(-2j * self._gamma)
		return [
			[1 + 0j, 0j, 0j, 0j],
			[0j, f1 * a, f2 * b, 0j],
			[0j, f3 * b, f4 * a, 0j],
			[0j, 0j, 0j, f5 * c],
		]

	def _resolve_parameters_(self, resolver: ParamResolver) -> PhasedFSimGate:
		if not self._is_parameterized_():
			return self
		return PhasedFSimGate(
			*(resolve_value(a, resolver) for a in (self._theta, self._zeta, self._chi, self._gamma, self._phi))
		)

	def _value_equality_values_(self) -> Any:
		return (self._theta, self._zeta, self._chi, self._gamma, self._phi)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("PhFSim", "PhFSim")

	def __str__(self) -> str:
		return f"PhasedFSimGate(theta={self._theta}, zeta={self._zeta}, chi={self._chi}, gamma={self._gamma}, phi={self._phi})"

	def __repr__(self) -> str:
		return (
			f"vqubit.PhasedFSimGate(theta={self._theta!r}, zeta={self._zeta!r}, chi={self._chi!r}, "
			f"gamma={self._gamma!r}, phi={self._phi!r})"
		)


@value_equality
class UniformSuperpositionGate(Gate):
	"""Prepares the uniform superposition over the states ``[0, M)`` on ``num_qubits``."""

	def __init__(self, m_value: int, num_qubits: int) -> None:
		if not (isinstance(m_value, int) and m_value > 0):
			raise ValueError("m_value must be a positive integer.")
		bits = m_value.bit_length()
		if (m_value & (m_value - 1)) == 0:
			bits -= 1
		if not (isinstance(num_qubits, int) and num_qubits >= bits):
			raise ValueError("num_qubits must be an integer greater than or equal to log2(m_value).")
		self._m_value = m_value
		self._num_qubits = num_qubits

	@property
	def m_value(self) -> int:
		return self._m_value

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _decompose_(self, qubits: Sequence[Qid]) -> Any:
		from ._gates import H, Ry, X

		qreg = list(qubits)[::-1]
		m = self._m_value
		if m == 1:
			return
		if (m & (m - 1)) == 0:
			for q in qreg[: m.bit_length() - 1]:
				yield H(q)
			return
		k = m.bit_length()
		l_value = [i for i in range(m.bit_length()) if (m >> i) & 1]
		for q_bit in l_value[1:k]:
			yield X(qreg[q_bit])
		m_current = 2 ** l_value[0]
		theta = -2 * math.acos(math.sqrt(m_current / m))
		if l_value[0] > 0:
			for q in qreg[: l_value[0]]:
				yield H(q)
		yield Ry(theta)(qreg[l_value[1]])
		for i in range(l_value[0], l_value[1]):
			yield H.controlled_by(qreg[l_value[1]], control_values=[False]).on(qreg[l_value[1]], qreg[i])
		for index in range(1, len(l_value) - 1):
			theta = -2 * math.acos(math.sqrt(2 ** l_value[index] / (m - m_current)))
			yield (
				Ry(theta)
				.controlled_by(qreg[l_value[index]], control_values=[0])
				.on(qreg[l_value[index]], qreg[l_value[index + 1]])
			)
			for i in range(l_value[index], l_value[index + 1]):
				yield H.controlled_by(qreg[l_value[index + 1]], control_values=[0]).on(
					qreg[l_value[index + 1]], qreg[i]
				)
			m_current += 2 ** l_value[index]

	def _value_equality_values_(self) -> Any:
		return (self._m_value, self._num_qubits)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("USup",) * self._num_qubits

	def __str__(self) -> str:
		return f"UniformSuperpositionGate(m_value={self._m_value}, num_qubits={self._num_qubits})"

	def __repr__(self) -> str:
		return f"vqubit.UniformSuperpositionGate(m_value={self._m_value}, num_qubits={self._num_qubits})"


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
	"PhaseGradientGate",
	"PhasedFSimGate",
	"PhasedISwapPowGate",
	"QuantumFourierTransformGate",
	"QubitPermutationGate",
	"SwapPowGate",
	"ThreeQubitDiagonalGate",
	"TwoQubitDiagonalGate",
	"UniformSuperpositionGate",
	"cphase",
	"givens",
	"global_phase_operation",
	"identity_each",
	"ms",
	"q",
	"qft",
	"reset_each",
	"riswap",
	"rx",
	"ry",
	"rz",
	"wait",
]
