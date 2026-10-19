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
import contextlib
import math
from collections.abc import Iterable, Sequence
from typing import Any

from ._linalg import Matrix, identity, kron, to_flat
from ._ops import GateOperation
from ._parameter import ParamResolver, is_symbolic, resolve_value
from ._qubits import Qid
from ._value import value_equality

_EPS = 1e-12
_MISSING = object()


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

	def _cached_unitary_(self) -> Matrix | None:
		"""Like ``_unitary_`` but memoised on the gate (gates are immutable values)."""
		cached: Any = getattr(self, "_unitary_cache", _MISSING)
		if cached is _MISSING:
			cached = self._unitary_()
			with contextlib.suppress(AttributeError):  # pragma: no cover - slotted gate
				object.__setattr__(self, "_unitary_cache", cached)
		return cached

	def _cached_flat_unitary_(self) -> list[float] | None:
		"""The native interleaved form of :meth:`_cached_unitary_`, memoised too."""
		flat: Any = getattr(self, "_flat_cache", _MISSING)
		if flat is _MISSING:
			matrix = self._cached_unitary_()
			flat = None if matrix is None else to_flat(matrix)
			with contextlib.suppress(AttributeError):  # pragma: no cover - slotted gate
				object.__setattr__(self, "_flat_cache", flat)
		return flat

	def _inverse_(self) -> Gate | None:
		"""Return the inverse gate, or ``None`` when there is none."""
		return None

	def _is_measurement_(self) -> bool:
		return False

	def _measurement_key_name_(self) -> str | None:
		return None

	def _measurement_key_obj_(self) -> Any:
		return None

	def _is_parameterized_(self) -> bool:
		return False

	def _resolve_parameters_(self, resolver: ParamResolver) -> Gate:
		return self

	def _decompose_(self, qubits: Sequence[Qid]) -> list[GateOperation] | None:
		"""Return a decomposition into simpler operations, or ``None``."""
		return None

	def _kraus_(self) -> list[Matrix] | None:
		"""Return the Kraus operators of a channel, or ``None`` for a gate."""
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

	def controlled_by(
		self,
		*control_qubits: Qid | Iterable[Qid],
		control_values: Sequence[bool | int] | None = None,
	) -> ControlledGate:
		"""Return this gate with additional control qubits.

		The number of controls comes from the qubits passed here; the actual
		qubits are supplied when the returned gate is applied with ``on``, the
		controls first.
		"""
		controls: list[Qid] = []
		for qubit in control_qubits:
			if isinstance(qubit, Qid):
				controls.append(qubit)
			else:
				controls.extend(qubit)
		values = None if control_values is None else tuple(bool(v) for v in control_values)
		num_controls = len(controls) if controls else (len(values) if values is not None else 0)
		return ControlledGate(self, num_controls=num_controls, control_values=values)

	def __pow__(self, power: float) -> Gate:
		raise NotImplementedError(f"{type(self).__name__} does not support powers")

	def __call__(self, *qubits: Qid | Iterable[Qid]) -> GateOperation:
		return self.on(*qubits)

	def __str__(self) -> str:
		return type(self).__name__


class InterchangeableQubitsGate:
	"""Marker mixin for gates that are symmetric under permuting their qubits.

	Subclassing this (alongside :class:`Gate`) tells optimizers that the order of
	the gate's qubits does not matter.
	"""


@value_equality
class EigenGate(Gate):
	"""A single-qubit gate defined by its eigencomponents and an exponent."""

	def __init__(self, *, exponent: Any = 1.0, global_shift: Any = 0.0) -> None:
		self._exponent = exponent if is_symbolic(exponent) else float(exponent)
		self._global_shift = global_shift if is_symbolic(global_shift) else float(global_shift)

	@property
	def exponent(self) -> float:
		return self._exponent

	@property
	def global_shift(self) -> float:
		return self._global_shift

	def _eigen_components(self) -> list[tuple[float, Matrix]]:
		raise NotImplementedError

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
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

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._exponent) or is_symbolic(self._global_shift)

	def _resolve_parameters_(self, resolver: ParamResolver) -> EigenGate:
		if not self._is_parameterized_():
			return self
		return type(self)(
			exponent=resolve_value(self._exponent, resolver),
			global_shift=resolve_value(self._global_shift, resolver),
		)

	def _value_equality_values_(self) -> Any:
		return (self._exponent, self._global_shift)

	def __pow__(self, power: float) -> EigenGate:
		return type(self)(exponent=self._exponent * power, global_shift=self._global_shift)

	def _inverse_(self) -> EigenGate:
		return type(self)(exponent=-self._exponent, global_shift=self._global_shift)

	def __str__(self) -> str:
		return f"{type(self).__name__}(exponent={_describe(self._exponent)})"

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}(exponent={self._exponent!r}, global_shift={self._global_shift!r})"


def _pow_coefficients(exponent: float, global_shift: float) -> tuple[complex, complex]:
	"""Return the ``A`` (eigenvalue ``+1``) and ``B`` (eigenvalue ``-1``) coefficients.

	``complex(-1.0) ** exponent`` gives exactly ``+-1`` for integer exponents, so
	this matches the eigencomponent expansion bit for bit.
	"""
	return cmath.exp(1j * math.pi * global_shift), complex(-1.0) ** exponent * cmath.exp(-1j * math.pi * global_shift)


def _xpow_entries(exponent: float) -> tuple[complex, complex]:
	"""Return the ``(diagonal, off_diagonal)`` entries of ``XPowGate(exponent, 0)``."""
	a, b = _pow_coefficients(exponent, 0.0)
	return (a + b) / 2, (a - b) / 2


def _clean(value: float) -> float | int:
	"""Return ``value`` as an int when it is a whole number, else as a float."""
	return int(value) if abs(value - round(value)) < _EPS else float(value)


def _describe(value: Any) -> str:
	"""Render a numeric or symbolic parameter for a gate's ``str``."""
	return str(value) if is_symbolic(value) else str(_clean(value))


def _exponent_suffix(exponent: Any) -> str:
	"""The ``^0.5`` a Cirq diagram appends to a power gate's symbol.

	An exponent of ``1`` adds nothing, so the plain gates keep the bare symbol
	they have always had: ``vq.CZ`` is still ``@`` and ``vq.CZ ** 0.5`` is
	``@^0.5``, exactly as in Cirq.
	"""
	if is_symbolic(exponent):
		return f"^{exponent}"
	if abs(float(exponent) - 1.0) < _EPS:
		return ""
	return f"^{_clean(exponent)}"


def is_plain_gate(gate: Any) -> bool:
	"""Whether ``gate`` is the unexponentiated form of a power gate.

	The plain gates *are* power gates here (``vq.CZ`` is a ``CZPowGate``), so a
	table that maps a gate onto a circuit instruction has to ask this before
	mapping: only the unexponentiated gate *is* the ``CZ``/``CX``/``SWAP``
	instruction, not ``vq.CZ ** 0.5``. A gate that is not a power gate at all is
	not this function's business, so it counts as plain and the caller's own
	``isinstance`` decides.
	"""
	if not isinstance(gate, _PowGate):
		return True
	return (
		not is_symbolic(gate._exponent)
		and not is_symbolic(gate._global_shift)
		and abs(float(gate._exponent) - 1.0) < _EPS
		and abs(float(gate._global_shift)) < _EPS
	)


class XPowGate(EigenGate):
	"""The power of the Pauli ``X`` gate."""

	def _num_qubits_(self) -> int:
		return 1

	def _eigen_components(self) -> list[tuple[float, Matrix]]:
		return [
			(1.0, [[0.5, 0.5], [0.5, 0.5]]),
			(-1.0, [[0.5, -0.5], [-0.5, 0.5]]),
		]

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		a, b = _pow_coefficients(self._exponent, self._global_shift)
		diagonal = (a + b) / 2
		off = (a - b) / 2
		return [[diagonal, off], [off, diagonal]]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("X" + _exponent_suffix(self._exponent),)


class YPowGate(EigenGate):
	"""The power of the Pauli ``Y`` gate."""

	def _num_qubits_(self) -> int:
		return 1

	def _eigen_components(self) -> list[tuple[float, Matrix]]:
		return [
			(1.0, [[0.5, -0.5j], [0.5j, 0.5]]),
			(-1.0, [[0.5, 0.5j], [-0.5j, 0.5]]),
		]

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		a, b = _pow_coefficients(self._exponent, self._global_shift)
		diagonal = (a + b) / 2
		off = 0.5j * (a - b)
		return [[diagonal, -off], [off, diagonal]]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Y" + _exponent_suffix(self._exponent),)


class ZPowGate(EigenGate):
	"""The power of the Pauli ``Z`` gate."""

	def _num_qubits_(self) -> int:
		return 1

	def _eigen_components(self) -> list[tuple[float, Matrix]]:
		return [
			(1.0, [[1.0, 0.0], [0.0, 0.0]]),
			(-1.0, [[0.0, 0.0], [0.0, 1.0]]),
		]

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		a, b = _pow_coefficients(self._exponent, self._global_shift)
		return [[a, 0j], [0j, b]]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		if not is_symbolic(self._exponent):
			if abs(self._exponent - 0.5) < _EPS:
				return ("S",)
			if abs(self._exponent - 0.25) < _EPS:
				return ("T",)
		return ("Z" + _exponent_suffix(self._exponent),)


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

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		a, b = _pow_coefficients(self._exponent, self._global_shift)
		s = math.sqrt(0.5)
		top_left = a * (1 + s) / 2 + b * (1 - s) / 2
		bottom_right = a * (1 - s) / 2 + b * (1 + s) / 2
		off = s / 2 * (a - b)
		return [[top_left, off], [off, bottom_right]]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("H" + _exponent_suffix(self._exponent),)


# --- multi-qubit power gates ------------------------------------------------
#
# In Cirq the plain multi-qubit gates *are* these classes at exponent 1:
# ``cirq.CZ`` is a ``CZPowGate``, so ``isinstance(cirq.CZ, cirq.CZPowGate)``
# holds and ``cirq.CZ == cirq.CZPowGate()``. That is true here too, which is why
# the plain gate classes below are aliases of these, declared next to the
# single-qubit power gates rather than in ``_gates_more``.
#
# They are *not* ``EigenGate`` subclasses: vqubit's ``EigenGate._unitary_`` is
# 2x2 by construction, so it cannot hold a multi-qubit gate. That is the
# documented, permanent ``*-EigenGate`` divergence.


@value_equality
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
		return ("@",) * self._NumControls + (self._Symbol + _exponent_suffix(self._exponent),)


class CXPowGate(_ControlledPowGate):
	"""A power of the controlled-``X`` gate."""

	_SubGate = XPowGate
	_Symbol = "X"


class CYPowGate(_ControlledPowGate):
	"""A power of the controlled-``Y`` gate."""

	_SubGate = YPowGate
	_Symbol = "Y"


class CZPowGate(_ControlledPowGate):
	"""A power of the controlled-``Z`` gate, which draws as two ``@``."""

	_SubGate = ZPowGate
	_Symbol = "@"


class CCXPowGate(_ControlledPowGate):
	"""A power of the doubly-controlled ``X`` gate."""

	_SubGate = XPowGate
	_NumControls = 2
	_Symbol = "X"

	def _decompose_(self, qubits: Sequence[Qid]) -> list[GateOperation]:
		"""The eleven-gate H/T/CNOT circuit, which only exists at exponent 1."""
		if not is_plain_gate(self):
			raise ValueError(f"cannot decompose {self!r}: the Toffoli circuit only applies at exponent 1")
		control0, control1, target = qubits
		return [
			H.on(target),
			CNOT.on(control1, target),
			(T**-1).on(target),
			CNOT.on(control0, target),
			T.on(target),
			CNOT.on(control1, target),
			(T**-1).on(target),
			CNOT.on(control0, target),
			T.on(control1),
			T.on(target),
			H.on(target),
			CNOT.on(control0, control1),
			T.on(control0),
			(T**-1).on(control1),
			CNOT.on(control0, control1),
		]


class CCYPowGate(_ControlledPowGate):
	"""A power of the doubly-controlled ``Y`` gate."""

	_SubGate = YPowGate
	_NumControls = 2
	_Symbol = "Y"


class CCZPowGate(_ControlledPowGate):
	"""A power of the doubly-controlled ``Z`` gate, which draws as three ``@``."""

	_SubGate = ZPowGate
	_NumControls = 2
	_Symbol = "@"


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
		return ("×", "×" + _exponent_suffix(self._exponent))


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
		return ("iSwap", "iSwap" + _exponent_suffix(self._exponent))


class SqrtISwapGate(ISwapPowGate):
	"""The square root of the ``iSWAP`` gate, at exponent ``0.5`` by default.

	Cirq has no such class -- ``cirq.SQRT_ISWAP`` is a plain ``ISwapPowGate``,
	and so is ``vq.SQRT_ISWAP`` -- so this name is vqubit's own, kept because
	``vq.SqrtISwapGate()`` used to be the class of ``vq.SQRT_ISWAP``. Being a
	subclass, an instance of it still passes ``isinstance(gate, vq.ISwapPowGate)``.
	"""

	def __init__(self, *, exponent: Any = 0.5, global_shift: Any = 0.0) -> None:
		super().__init__(exponent=exponent, global_shift=global_shift)


@value_equality
class PhasedXPowGate(Gate):
	"""A power of ``X`` conjugated by a phase: ``Z^p X^e Z^-p``."""

	def __init__(self, *, phase_exponent: Any = 0.0, exponent: Any = 1.0, global_shift: Any = 0.0) -> None:
		self._phase_exponent = phase_exponent if is_symbolic(phase_exponent) else float(phase_exponent)
		self._exponent = exponent if is_symbolic(exponent) else float(exponent)
		self._global_shift = global_shift if is_symbolic(global_shift) else float(global_shift)

	@property
	def phase_exponent(self) -> Any:
		return self._phase_exponent

	@property
	def exponent(self) -> Any:
		return self._exponent

	@property
	def global_shift(self) -> Any:
		return self._global_shift

	def _num_qubits_(self) -> int:
		return 1

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		diagonal, off = _xpow_entries(self._exponent)
		z = cmath.exp(1j * math.pi * self._phase_exponent)
		z_inv = cmath.exp(-1j * math.pi * self._phase_exponent)
		shift = cmath.exp(1j * math.pi * self._global_shift * self._exponent) if self._global_shift else 1 + 0j
		return [[shift * diagonal, shift * off * z_inv], [shift * off * z, shift * diagonal]]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._phase_exponent) or is_symbolic(self._exponent) or is_symbolic(self._global_shift)

	def _resolve_parameters_(self, resolver: ParamResolver) -> PhasedXPowGate:
		if not self._is_parameterized_():
			return self
		return PhasedXPowGate(
			phase_exponent=resolve_value(self._phase_exponent, resolver),
			exponent=resolve_value(self._exponent, resolver),
			global_shift=resolve_value(self._global_shift, resolver),
		)

	def _value_equality_values_(self) -> Any:
		return (self._phase_exponent, self._exponent, self._global_shift)

	def __pow__(self, power: float) -> PhasedXPowGate:
		return PhasedXPowGate(
			phase_exponent=self._phase_exponent,
			exponent=self._exponent * power,
			global_shift=self._global_shift,
		)

	def _inverse_(self) -> PhasedXPowGate:
		return PhasedXPowGate(
			phase_exponent=self._phase_exponent,
			exponent=-self._exponent,
			global_shift=self._global_shift,
		)

	def _decompose_(self, qubits: Sequence[Qid]) -> list[GateOperation]:
		qubit = qubits[0]
		return [
			ZPowGate(exponent=self._phase_exponent).on(qubit),
			XPowGate(exponent=self._exponent).on(qubit),
			ZPowGate(exponent=-self._phase_exponent).on(qubit),
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("PhX",)

	def __str__(self) -> str:
		return f"PhasedXPowGate(phase_exponent={_describe(self._phase_exponent)}, exponent={_describe(self._exponent)})"

	def __repr__(self) -> str:
		return f"vqubit.PhasedXPowGate(phase_exponent={self._phase_exponent!r}, exponent={self._exponent!r}, global_shift={self._global_shift!r})"


@value_equality
class PhasedXZGate(Gate):
	"""The standard single-qubit gate: ``Z^z`` after a phased ``X^x``."""

	def __init__(self, *, x_exponent: Any = 0.0, z_exponent: Any = 0.0, axis_phase_exponent: Any = 0.0) -> None:
		self._x_exponent = x_exponent if is_symbolic(x_exponent) else float(x_exponent)
		self._z_exponent = z_exponent if is_symbolic(z_exponent) else float(z_exponent)
		self._axis_phase_exponent = (
			axis_phase_exponent if is_symbolic(axis_phase_exponent) else float(axis_phase_exponent)
		)

	@property
	def x_exponent(self) -> Any:
		return self._x_exponent

	@property
	def z_exponent(self) -> Any:
		return self._z_exponent

	@property
	def axis_phase_exponent(self) -> Any:
		return self._axis_phase_exponent

	def _num_qubits_(self) -> int:
		return 1

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		diagonal, off = _xpow_entries(self._x_exponent)
		axis = cmath.exp(1j * math.pi * self._axis_phase_exponent)
		axis_inv = cmath.exp(-1j * math.pi * self._axis_phase_exponent)
		z = cmath.exp(1j * math.pi * self._z_exponent)
		return [[diagonal, off * axis_inv], [z * off * axis, z * diagonal]]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._x_exponent) or is_symbolic(self._z_exponent) or is_symbolic(self._axis_phase_exponent)

	def _resolve_parameters_(self, resolver: ParamResolver) -> PhasedXZGate:
		if not self._is_parameterized_():
			return self
		return PhasedXZGate(
			x_exponent=resolve_value(self._x_exponent, resolver),
			z_exponent=resolve_value(self._z_exponent, resolver),
			axis_phase_exponent=resolve_value(self._axis_phase_exponent, resolver),
		)

	def _value_equality_values_(self) -> Any:
		return (self._x_exponent, self._z_exponent, self._axis_phase_exponent)

	def _inverse_(self) -> PhasedXZGate:
		# (Z^z PhX(x, a))^-1 = Z^-z PhX(-x, a + z).
		return PhasedXZGate(
			x_exponent=-self._x_exponent,
			z_exponent=-self._z_exponent,
			axis_phase_exponent=self._axis_phase_exponent + self._z_exponent,
		)

	def _decompose_(self, qubits: Sequence[Qid]) -> list[GateOperation]:
		qubit = qubits[0]
		return [
			PhasedXPowGate(phase_exponent=self._axis_phase_exponent, exponent=self._x_exponent).on(qubit),
			ZPowGate(exponent=self._z_exponent).on(qubit),
		]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("PhXZ",)

	def __str__(self) -> str:
		return (
			f"PhasedXZGate(x_exponent={_describe(self._x_exponent)}, "
			f"z_exponent={_describe(self._z_exponent)}, "
			f"axis_phase_exponent={_describe(self._axis_phase_exponent)})"
		)

	def __repr__(self) -> str:
		return f"vqubit.PhasedXZGate(x_exponent={self._x_exponent!r}, z_exponent={self._z_exponent!r}, axis_phase_exponent={self._axis_phase_exponent!r})"


@value_equality
class _RotationGate(Gate):
	"""Shared machinery for the single-axis rotation gates."""

	_axis_symbol = "R"

	def __init__(self, radians: Any) -> None:
		self._radians = radians if is_symbolic(radians) else float(radians)

	@property
	def radians(self) -> Any:
		return self._radians

	def _num_qubits_(self) -> int:
		return 1

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		return self._rotation_matrix(self._radians)

	def _rotation_matrix(self, radians: float) -> Matrix:
		raise NotImplementedError

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._radians)

	def _resolve_parameters_(self, resolver: ParamResolver) -> _RotationGate:
		if not self._is_parameterized_():
			return self
		return type(self)(resolve_value(self._radians, resolver))

	def _value_equality_values_(self) -> Any:
		return (self._radians,)

	def _inverse_(self) -> _RotationGate:
		return type(self)(-self._radians)

	def __pow__(self, power: float) -> _RotationGate:
		"""Scale the angle, as ``cirq.Rz(rads) ** x`` does (``0.3 -> 0.15``)."""
		return type(self)(self._radians * power)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return (self._axis_symbol,)

	def __str__(self) -> str:
		return f"{self._axis_symbol}({_describe(self._radians)})"

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}({self._radians!r})"


class RxGate(_RotationGate):
	"""A rotation of ``radians`` around the X axis."""

	_axis_symbol = "Rx"

	def _rotation_matrix(self, radians: float) -> Matrix:
		c = math.cos(radians / 2)
		s = math.sin(radians / 2)
		return [[c + 0j, -1j * s], [-1j * s, c + 0j]]


class RyGate(_RotationGate):
	"""A rotation of ``radians`` around the Y axis."""

	_axis_symbol = "Ry"

	def _rotation_matrix(self, radians: float) -> Matrix:
		c = math.cos(radians / 2)
		s = math.sin(radians / 2)
		return [[c + 0j, -s + 0j], [s + 0j, c + 0j]]


class RzGate(_RotationGate):
	"""A rotation of ``radians`` around the Z axis."""

	_axis_symbol = "Rz"

	def _rotation_matrix(self, radians: float) -> Matrix:
		return [[cmath.exp(-0.5j * radians), 0j], [0j, cmath.exp(0.5j * radians)]]


@value_equality
class IdentityGate(Gate):
	"""The identity on ``num_qubits`` qubits."""

	def __init__(self, num_qubits: int = 1) -> None:
		self._num_qubits = int(num_qubits)

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _unitary_(self) -> Matrix | None:
		return identity(1 << self._num_qubits)

	def _inverse_(self) -> IdentityGate:
		return IdentityGate(self._num_qubits)

	def __pow__(self, power: float) -> IdentityGate:
		"""The identity is its own power, as ``cirq.I ** x`` is ``cirq.I``."""
		return IdentityGate(self._num_qubits)

	def _value_equality_values_(self) -> Any:
		return (self._num_qubits,)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("I",) * self._num_qubits

	def __repr__(self) -> str:
		return f"vqubit.IdentityGate({self._num_qubits})"


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

	def _is_measurement_(self) -> bool:
		return True

	def _measurement_key_name_(self) -> str | None:
		return self._key

	def _measurement_key_obj_(self) -> Any:
		from ._measurement_key import MeasurementKey

		return MeasurementKey(self._key)

	def __str__(self) -> str:
		return f"measure({self._key!r})"

	def __repr__(self) -> str:
		return (
			f"vqubit.MeasurementGate({self._key!r}, num_qubits={self._num_qubits}, invert_mask={self._invert_mask!r})"
		)


@value_equality
class _FixedGate(Gate):
	"""Shared machinery for the fixed multi-qubit gates.

	``CSWAP`` is the only plain gate left that is not a power gate (``vq.CZ``
	and friends are ``CZPowGate``/``SwapPowGate``/… at exponent 1, as in Cirq),
	so this is what it builds on. The matrix and the diagram symbols are produced
	by static methods rather than stored as class attributes, so a gate instance
	never hands out a shared, mutable matrix.
	"""

	@staticmethod
	def _base_matrix() -> Matrix:
		raise NotImplementedError

	@staticmethod
	def _base_symbols() -> tuple[str, ...]:
		raise NotImplementedError

	@staticmethod
	def _base_num_qubits() -> int:
		return 2

	def _num_qubits_(self) -> int:
		return self._base_num_qubits()

	def _unitary_(self) -> Matrix:
		return [list(row) for row in self._base_matrix()]

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return self._base_symbols()

	def _value_equality_values_(self) -> Any:
		return ()

	def _inverse_(self) -> Gate:
		return type(self)()

	def __str__(self) -> str:
		return type(self).__name__

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}()"


# The plain gates *are* the power gates at exponent 1, as in Cirq: ``vq.CZ`` is
# a ``CZPowGate``, so ``isinstance(vq.CZ, vq.CZPowGate)`` and
# ``vq.CZ == vq.CZPowGate()`` both hold. The historical class names stay as
# aliases, so ``vq.CZGate`` and ``vq.OpIdentifier(vq.CNOTGate, 2)`` keep working.
#
# Tables that map a gate to a circuit instruction (QASM, Clifford, JSON, ...)
# must therefore use ``is_plain_gate`` below: after this fold
# ``isinstance(vq.CZ ** 0.5, vq.CZGate)`` is ``True``, and a fractional power is
# not the ``CZ`` instruction.
CNOTGate = CXPowGate
CZGate = CZPowGate
SWAPGate = SwapPowGate
CCXGate = CCXPowGate
ISWAPGate = ISwapPowGate


def _swap_basis_states(num_qubits: int, swaps: Sequence[tuple[int, int]]) -> Matrix:
	"""Return the identity on ``num_qubits`` with the listed basis states swapped."""
	matrix = identity(1 << num_qubits)
	for a, b in swaps:
		matrix[a][a] = 0j
		matrix[b][b] = 0j
		matrix[a][b] = 1.0 + 0j
		matrix[b][a] = 1.0 + 0j
	return matrix


class CSWAPGate(_FixedGate):
	"""The Fredkin gate: swap two qubits when the control is ``1``."""

	@staticmethod
	def _base_num_qubits() -> int:
		return 3

	@staticmethod
	def _base_symbols() -> tuple[str, ...]:
		return ("@", "×", "×")

	@staticmethod
	def _base_matrix() -> Matrix:
		# Swap |101> and |110>, which share the control at 1.
		return _swap_basis_states(3, [(0b101, 0b110)])

	def _decompose_(self, qubits: Sequence[Qid]) -> list[GateOperation]:
		control, first, second = qubits
		return [
			CNOT.on(second, first),
			CCX.on(control, first, second),
			CNOT.on(second, first),
		]


_PAULI_2X2: dict[str, Matrix] = {
	"X": [[0j, 1 + 0j], [1 + 0j, 0j]],
	"Y": [[0j, -1j], [1j, 0j]],
	"Z": [[1 + 0j, 0j], [0j, -1 + 0j]],
}


@value_equality
class _PauliProductPowGate(Gate):
	"""A two-qubit gate that is a power of a product of two Paulis."""

	_letter = "X"
	_symbols: tuple[str, ...] = ()

	def __init__(self, *, exponent: Any = 1.0) -> None:
		self._exponent = exponent if is_symbolic(exponent) else float(exponent)

	@property
	def exponent(self) -> Any:
		return self._exponent

	def _num_qubits_(self) -> int:
		return 2

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		pauli = kron(_PAULI_2X2[self._letter], _PAULI_2X2[self._letter])
		# Eigenvalues +1 and -1, so the power is a projection combination with no
		# global phase: at exponent 1 this is the Pauli product itself.
		phase = cmath.exp(1j * math.pi * self._exponent)
		a = (1 + phase) / 2
		b = (1 - phase) / 2
		return [[(a if i == j else 0j) + b * pauli[i][j] for j in range(4)] for i in range(4)]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._exponent)

	def _resolve_parameters_(self, resolver: ParamResolver) -> _PauliProductPowGate:
		if not self._is_parameterized_():
			return self
		return type(self)(exponent=resolve_value(self._exponent, resolver))

	def _value_equality_values_(self) -> Any:
		return (self._exponent,)

	def __pow__(self, power: float) -> _PauliProductPowGate:
		return type(self)(exponent=self._exponent * power)

	def _inverse_(self) -> _PauliProductPowGate:
		return type(self)(exponent=-self._exponent)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		# The exponent rides on the last symbol only, as in Cirq: ``XX**0.5``
		# draws as ``XX`` over ``XX^0.5``.
		return (*self._symbols[:-1], self._symbols[-1] + _exponent_suffix(self._exponent))

	def __str__(self) -> str:
		return f"{type(self).__name__}(exponent={_describe(self._exponent)})"


class XXPowGate(_PauliProductPowGate):
	"""A power of ``X⊗X``, an Ising interaction."""

	_letter = "X"
	_symbols = ("XX", "XX")


class YYPowGate(_PauliProductPowGate):
	"""A power of ``Y⊗Y``."""

	_letter = "Y"
	_symbols = ("YY", "YY")


class ZZPowGate(_PauliProductPowGate):
	"""A power of ``Z⊗Z``."""

	_letter = "Z"
	_symbols = ("ZZ", "ZZ")


@value_equality
class MatrixGate(Gate):
	"""A gate defined directly by a unitary matrix.

	The matrix dimension must be a power of two; the gate acts on
	``log2(dim)`` qubits.
	"""

	def __init__(self, matrix: Matrix) -> None:
		rows = [[complex(entry) for entry in row] for row in matrix]
		dim = len(rows)
		if dim == 0 or dim & (dim - 1):
			raise ValueError("MatrixGate needs a matrix whose dimension is a power of two")
		if any(len(row) != dim for row in rows):
			raise ValueError("MatrixGate needs a square matrix")
		self._matrix = rows

	@property
	def matrix(self) -> Matrix:
		return [list(row) for row in self._matrix]

	def _num_qubits_(self) -> int:
		return len(self._matrix).bit_length() - 1

	def _unitary_(self) -> Matrix:
		return [list(row) for row in self._matrix]

	def _value_equality_values_(self) -> Any:
		return (tuple(tuple(row) for row in self._matrix),)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Matrix",) * self.num_qubits

	def __str__(self) -> str:
		return f"MatrixGate({self.num_qubits}q)"

	def __repr__(self) -> str:
		return f"vqubit.MatrixGate({self._matrix!r})"


@value_equality
class ControlledGate(Gate):
	"""A gate applied only when the control qubits have the given values."""

	def __init__(
		self,
		sub_gate: Gate,
		*,
		num_controls: int,
		control_values: Sequence[bool | int] | None = None,
	) -> None:
		if num_controls < 1:
			raise ValueError("ControlledGate needs at least one control qubit")
		if not isinstance(sub_gate, Gate):
			raise TypeError("ControlledGate needs a Gate to control")
		values = (True,) * num_controls if control_values is None else tuple(bool(v) for v in control_values)
		if len(values) != num_controls:
			raise ValueError(f"control_values has length {len(values)}, expected {num_controls}")
		self._sub_gate = sub_gate
		self._num_controls = num_controls
		self._control_values = values

	@property
	def sub_gate(self) -> Gate:
		return self._sub_gate

	@property
	def num_controls(self) -> int:
		return self._num_controls

	@property
	def control_values(self) -> tuple[bool, ...]:
		return self._control_values

	def _num_qubits_(self) -> int:
		return self._num_controls + self._sub_gate.num_qubits

	def _unitary_(self) -> Matrix | None:
		sub = self._sub_gate._unitary_()
		if sub is None:
			return None
		sub = [list(row) for row in sub]
		block = len(sub)
		blocks = 1 << self._num_controls
		total = blocks * block
		target = 0
		for value in self._control_values:
			target = (target << 1) | int(value)
		offset = target * block
		matrix = identity(total)
		for i in range(block):
			for j in range(block):
				matrix[offset + i][offset + j] = sub[i][j]
		return matrix

	def _inverse_(self) -> ControlledGate | None:
		sub_inverse = self._sub_gate._inverse_()
		if sub_inverse is None:
			return None
		return ControlledGate(sub_inverse, num_controls=self._num_controls, control_values=self._control_values)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		controls = tuple("@" if value else "○" for value in self._control_values)
		return controls + self._sub_gate._diagram_symbols_(qubits[self._num_controls :])

	def _value_equality_values_(self) -> Any:
		return (self._sub_gate, self._num_controls, self._control_values)

	def __str__(self) -> str:
		return f"ControlledGate({self._sub_gate})"

	def __repr__(self) -> str:
		return f"vqubit.ControlledGate({self._sub_gate!r}, num_controls={self._num_controls}, control_values={self._control_values!r})"


# Instances, as in Cirq. The plain gates are the power gates at exponent 1, so
# ``CZ`` and ``CZPowGate()`` are the same value, as they are in Cirq.
X: XPowGate = XPowGate()
Y: YPowGate = YPowGate()
Z: ZPowGate = ZPowGate()
H: HPowGate = HPowGate()
S: ZPowGate = ZPowGate(exponent=0.5)
T: ZPowGate = ZPowGate(exponent=0.25)
I: IdentityGate = IdentityGate()  # noqa: E741 - the name matches Cirq's identity gate
CNOT: CXPowGate = CXPowGate()
CX: CXPowGate = CXPowGate()
CZ: CZPowGate = CZPowGate()
SWAP: SwapPowGate = SwapPowGate()
CCX: CCXPowGate = CCXPowGate()
TOFFOLI: CCXPowGate = CCXPowGate()
CSWAP: CSWAPGate = CSWAPGate()
FREDKIN: CSWAPGate = CSWAPGate()
CSwapGate = CSWAPGate
XX: XXPowGate = XXPowGate()
YY: YYPowGate = YYPowGate()
ZZ: ZZPowGate = ZZPowGate()
ISWAP: ISwapPowGate = ISwapPowGate()
CY: CYPowGate = CYPowGate()
CCNOT: CCXPowGate = CCXPowGate()
CCZ: CCZPowGate = CCZPowGate()
CCY: CCYPowGate = CCYPowGate()
ISWAP_INV: ISwapPowGate = ISwapPowGate(exponent=-1)
SQRT_ISWAP_INV: ISwapPowGate = ISwapPowGate(exponent=-0.5)
SWAP_INV: SwapPowGate = SwapPowGate(exponent=-1)


def Rx(radians: Any) -> RxGate:
	"""Return a rotation of ``radians`` around the X axis."""
	return RxGate(radians)


def Ry(radians: Any) -> RyGate:
	"""Return a rotation of ``radians`` around the Y axis."""
	return RyGate(radians)


def Rz(radians: Any) -> RzGate:
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


M = measure


@value_equality
class GlobalPhaseGate(Gate):
	"""A phase applied to the whole state, acting on no qubits."""

	def __init__(self, exponent: Any = 0.0) -> None:
		self._exponent = exponent if is_symbolic(exponent) else float(exponent)

	@property
	def exponent(self) -> Any:
		return self._exponent

	def _num_qubits_(self) -> int:
		return 0

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		return [[cmath.exp(1j * math.pi * self._exponent)]]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._exponent)

	def _resolve_parameters_(self, resolver: ParamResolver) -> GlobalPhaseGate:
		if not self._is_parameterized_():
			return self
		return GlobalPhaseGate(resolve_value(self._exponent, resolver))

	def _value_equality_values_(self) -> Any:
		return (self._exponent,)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ()

	def __str__(self) -> str:
		return f"GlobalPhaseGate(exponent={_describe(self._exponent)})"

	def __repr__(self) -> str:
		return f"vqubit.GlobalPhaseGate({self._exponent!r})"


@value_equality
class WaitGate(Gate):
	"""An idle wait of ``duration_nanos`` on ``num_qubits`` qubits."""

	def __init__(self, duration_nanos: float, num_qubits: int = 1) -> None:
		if duration_nanos < 0:
			raise ValueError("WaitGate duration must be non-negative")
		self._duration = float(duration_nanos)
		self._num_qubits = int(num_qubits)

	@property
	def duration_nanos(self) -> float:
		return self._duration

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _unitary_(self) -> Matrix:
		return identity(1 << self._num_qubits)

	def _value_equality_values_(self) -> Any:
		return (self._duration, self._num_qubits)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("Wait",) * self._num_qubits

	def __str__(self) -> str:
		return f"WaitGate({self._duration}, num_qubits={self._num_qubits})"

	def __repr__(self) -> str:
		return f"vqubit.WaitGate({self._duration!r}, num_qubits={self._num_qubits})"


def global_phase(exponent: Any) -> GlobalPhaseGate:
	"""Return a :class:`GlobalPhaseGate` with the given exponent."""
	return GlobalPhaseGate(exponent)


@value_equality
class FSimGate(Gate):
	"""The fSim gate: an iSWAP with a conditional phase, ``fSim(theta, phi)``."""

	def __init__(self, theta: Any, phi: Any) -> None:
		self._theta = theta if is_symbolic(theta) else float(theta)
		self._phi = phi if is_symbolic(phi) else float(phi)

	@property
	def theta(self) -> Any:
		return self._theta

	@property
	def phi(self) -> Any:
		return self._phi

	def _num_qubits_(self) -> int:
		return 2

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		c = math.cos(self._theta)
		s = math.sin(self._theta)
		phase = cmath.exp(-1j * self._phi)
		return [
			[1 + 0j, 0j, 0j, 0j],
			[0j, c + 0j, -1j * s, 0j],
			[0j, -1j * s, c + 0j, 0j],
			[0j, 0j, 0j, phase],
		]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._theta) or is_symbolic(self._phi)

	def _resolve_parameters_(self, resolver: ParamResolver) -> FSimGate:
		if not self._is_parameterized_():
			return self
		return FSimGate(resolve_value(self._theta, resolver), resolve_value(self._phi, resolver))

	def _inverse_(self) -> FSimGate:
		return FSimGate(-self._theta, -self._phi)

	def _value_equality_values_(self) -> Any:
		return (self._theta, self._phi)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("fSim", "fSim")

	def __str__(self) -> str:
		return f"FSimGate({_describe(self._theta)}, {_describe(self._phi)})"

	def __repr__(self) -> str:
		return f"vqubit.FSimGate({self._theta!r}, {self._phi!r})"


@value_equality
class MSGate(Gate):
	"""The Mølmer-Sørensen gate on two qubits: ``exp(-i rads X X)``."""

	def __init__(self, *, rads: Any) -> None:
		self._rads = rads if is_symbolic(rads) else float(rads)

	@property
	def rads(self) -> Any:
		return self._rads

	def _num_qubits_(self) -> int:
		return 2

	def _unitary_(self) -> Matrix | None:
		if self._is_parameterized_():
			return None
		x = _PAULI_2X2["X"]
		xx = kron(x, x)
		c = math.cos(self._rads)
		s = math.sin(self._rads)
		return [[(c if i == j else 0j) - 1j * s * xx[i][j] for j in range(4)] for i in range(4)]

	def _is_parameterized_(self) -> bool:
		return is_symbolic(self._rads)

	def _resolve_parameters_(self, resolver: ParamResolver) -> MSGate:
		if not self._is_parameterized_():
			return self
		return MSGate(rads=resolve_value(self._rads, resolver))

	def _value_equality_values_(self) -> Any:
		return (self._rads,)

	def _diagram_symbols_(self, qubits: Sequence[Qid]) -> tuple[str, ...]:
		return ("MS", "MS")

	def __str__(self) -> str:
		return f"MSGate(rads={_describe(self._rads)})"

	def __repr__(self) -> str:
		return f"vqubit.MSGate(rads={self._rads!r})"


SQRT_ISWAP: ISwapPowGate = ISwapPowGate(exponent=0.5)
