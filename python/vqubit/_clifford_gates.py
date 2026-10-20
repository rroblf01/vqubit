"""Clifford gates.

A Clifford gate maps the Pauli group to itself under conjugation. The
single-qubit Clifford group has 24 elements; :class:`SingleQubitCliffordGate`
represents one of them by the images of ``X`` and ``Z``. The multi-qubit
:class:`CliffordGate` stores its full unitary and can be built from a list of
operations or from a Clifford tableau.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from ._gates import CNOT, CZ, SWAP, Gate, H, I, S, X, Y, Z, identity
from ._linalg import Matrix, matmul
from ._linalg_public import to_nested_list
from ._pauli import Pauli, _product
from ._qubits import LineQubit, Qid
from ._value import value_equality

_PAULI_MATRIX: dict[Pauli, Matrix] = {
	Pauli.I: [[1 + 0j, 0j], [0j, 1 + 0j]],
	Pauli.X: [[0j, 1 + 0j], [1 + 0j, 0j]],
	Pauli.Y: [[0j, -1j], [1j, 0j]],
	Pauli.Z: [[1 + 0j, 0j], [0j, -1 + 0j]],
}
_PAULI_GATE: dict[Pauli, Any] = {Pauli.I: I, Pauli.X: X, Pauli.Y: Y, Pauli.Z: Z}


def _conj_t(m: Matrix) -> Matrix:
	return [[m[j][i].conjugate() for j in range(len(m))] for i in range(len(m[0]))]


def _matrix_close(a: Matrix, b: Matrix, atol: float = 1e-8) -> bool:
	return all(abs(a[i][j] - b[i][j]) <= atol for i in range(len(a)) for j in range(len(a[0])))


def _pauli_image(unitary: Matrix, pauli: Pauli) -> tuple[Pauli, bool]:
	"""Return ``(P, flip)`` with ``U pauli U† = (-1)**flip * P``."""
	image = matmul(matmul(unitary, _PAULI_MATRIX[pauli]), _conj_t(unitary))
	for candidate in Pauli:
		if _matrix_close(image, _PAULI_MATRIX[candidate]):
			return candidate, False
		if _matrix_close(image, [[-z for z in row] for row in _PAULI_MATRIX[candidate]]):
			return candidate, True
	raise ValueError("matrix is not a single-qubit Clifford gate")


def _canonical_key(unitary: Matrix) -> tuple[tuple[complex, ...], ...]:
	for row in unitary:
		for entry in row:
			if abs(entry) > 1e-9:
				phase = entry / abs(entry)
				return tuple(
					tuple(
						round((z / phase).real, 6) + 0j
						if abs((z / phase).imag) < 1e-6
						else round((z / phase).real, 6) + round((z / phase).imag, 6) * 1j
						for z in r
					)
					for r in unitary
				)
	return tuple(tuple(complex(z) for z in r) for r in unitary)


def _generate_single_qubit_cliffords() -> list[tuple[Pauli, bool, Pauli, bool, Matrix]]:
	generators = [to_nested_list(H._unitary_()), to_nested_list(S._unitary_())]
	seen: dict[Any, Matrix] = {_canonical_key(identity(2)): identity(2)}
	frontier = [identity(2)]
	while frontier:
		next_frontier: list[Matrix] = []
		for unitary in frontier:
			for generator in generators:
				candidate = matmul(unitary, generator)
				key = _canonical_key(candidate)
				if key not in seen:
					seen[key] = candidate
					next_frontier.append(candidate)
		frontier = next_frontier
	result = []
	for unitary in seen.values():
		x_to = _pauli_image(unitary, Pauli.X)
		z_to = _pauli_image(unitary, Pauli.Z)
		result.append((x_to[0], x_to[1], z_to[0], z_to[1], unitary))
	return result


_ALL_SINGLE_QUBIT_CLIFFORDS = _generate_single_qubit_cliffords()


@value_equality
class SingleQubitCliffordGate(Gate):
	"""A single-qubit Clifford gate, stored by the images of ``X`` and ``Z``."""

	I: SingleQubitCliffordGate  # noqa: E741 - the name matches Cirq's identity gate
	X: SingleQubitCliffordGate
	Y: SingleQubitCliffordGate
	Z: SingleQubitCliffordGate
	H: SingleQubitCliffordGate
	S: SingleQubitCliffordGate
	X_sqrt: SingleQubitCliffordGate
	Y_sqrt: SingleQubitCliffordGate
	Z_sqrt: SingleQubitCliffordGate
	X_nsqrt: SingleQubitCliffordGate
	Y_nsqrt: SingleQubitCliffordGate
	Z_nsqrt: SingleQubitCliffordGate
	all_single_qubit_cliffords: tuple[SingleQubitCliffordGate, ...]

	def __init__(self, *, x_to: tuple[Pauli, bool], z_to: tuple[Pauli, bool], unitary: Matrix) -> None:
		self._x_to = x_to
		self._z_to = z_to
		self._unitary_matrix = [list(row) for row in unitary]

	def _num_qubits_(self) -> int:
		return 1

	@staticmethod
	def from_xz_map(x_to: tuple[Pauli, bool], z_to: tuple[Pauli, bool]) -> SingleQubitCliffordGate:
		"""Return the gate mapping ``X`` to ``x_to`` and ``Z`` to ``z_to``."""
		for entry in _ALL_SINGLE_QUBIT_CLIFFORDS:
			if (entry[0], entry[1]) == (x_to[0], bool(x_to[1])) and (entry[2], entry[3]) == (z_to[0], bool(z_to[1])):
				return SingleQubitCliffordGate(x_to=entry[0:2], z_to=entry[2:4], unitary=entry[4])
		raise ValueError(f"No single-qubit Clifford maps X to {x_to} and Z to {z_to}")

	@staticmethod
	def from_single_map(
		pauli_map_to: dict[Pauli, tuple[Pauli, bool]] | None = None,
		*,
		x_to: tuple[Pauli, bool] | None = None,
		y_to: tuple[Pauli, bool] | None = None,
		z_to: tuple[Pauli, bool] | None = None,
	) -> SingleQubitCliffordGate:
		"""Return the gate from a single Pauli image, inferring the rest."""
		mapping = _merge_maps(pauli_map_to, x_to, y_to, z_to)
		if len(mapping) != 1:
			raise ValueError("from_single_map takes exactly one transform")
		((trans_from, (trans_to, flip)),) = tuple(mapping.items())
		if trans_from == trans_to:
			trans_from2 = _relative_pauli(trans_to, 1)
			trans_to2 = _relative_pauli(trans_from, 1)
			flip2 = False
		else:
			trans_from2 = trans_to
			trans_to2 = trans_from
			flip2 = not flip
		mapping[trans_from2] = (trans_to2, flip2)
		return SingleQubitCliffordGate.from_double_map(mapping)

	@staticmethod
	def from_double_map(
		pauli_map_to: dict[Pauli, tuple[Pauli, bool]] | None = None,
		*,
		x_to: tuple[Pauli, bool] | None = None,
		y_to: tuple[Pauli, bool] | None = None,
		z_to: tuple[Pauli, bool] | None = None,
	) -> SingleQubitCliffordGate:
		"""Return the gate from two Pauli images, inferring the third."""
		mapping = _merge_maps(pauli_map_to, x_to, y_to, z_to)
		if len(mapping) != 2:
			raise ValueError("from_double_map takes exactly two transforms")
		(from1, trans1), (from2, trans2) = tuple(mapping.items())
		from3 = _third_pauli(from1, from2)
		to3 = _third_pauli(trans1[0], trans2[0])
		flip3 = trans1[1] ^ trans2[1] ^ ((from1.value < from2.value) != (trans1[0].value < trans2[0].value))
		mapping[from3] = (to3, flip3)
		x_image = mapping[Pauli.X]
		z_image = mapping[Pauli.Z]
		return SingleQubitCliffordGate.from_xz_map(x_image, z_image)

	@staticmethod
	def from_pauli(pauli: Pauli, sqrt: bool = False) -> SingleQubitCliffordGate:
		"""Return the Pauli gate, or its square root."""
		prev_pauli = _relative_pauli(pauli, -1)
		next_pauli = _relative_pauli(pauli, 1)
		if sqrt:
			mapping = {prev_pauli: (next_pauli, True), pauli: (pauli, False), next_pauli: (prev_pauli, False)}
		else:
			mapping = {prev_pauli: (prev_pauli, True), pauli: (pauli, False), next_pauli: (next_pauli, True)}
		# The map is already complete, so read X and Z straight off it.
		return SingleQubitCliffordGate.from_xz_map(mapping[Pauli.X], mapping[Pauli.Z])

	@staticmethod
	def from_quarter_turns(pauli: Pauli, quarter_turns: int) -> SingleQubitCliffordGate:
		"""Return a quarter turn about a Pauli axis."""
		quarter_turns %= 4
		if quarter_turns == 0:
			return SingleQubitCliffordGate.I
		if quarter_turns == 1:
			return SingleQubitCliffordGate.from_pauli(pauli, True)
		if quarter_turns == 2:
			return SingleQubitCliffordGate.from_pauli(pauli)
		return SingleQubitCliffordGate.from_pauli(pauli, True) ** -1

	@staticmethod
	def from_unitary(u: Sequence[Sequence[complex]]) -> SingleQubitCliffordGate | None:
		"""Return the Clifford gate with the given unitary, or ``None``."""
		matrix = to_nested_list(u)
		if len(matrix) != 2 or len(matrix[0]) != 2:
			return None
		try:
			x_to = _pauli_image(matrix, Pauli.X)
			z_to = _pauli_image(matrix, Pauli.Z)
		except ValueError:
			return None
		return SingleQubitCliffordGate.from_xz_map(x_to, z_to)

	@staticmethod
	def from_unitary_with_global_phase(
		u: Sequence[Sequence[complex]],
	) -> tuple[SingleQubitCliffordGate, complex] | None:
		"""Return the Clifford gate and the global phase relating it to ``u``."""
		matrix = to_nested_list(u)
		gate = SingleQubitCliffordGate.from_unitary(matrix)
		if gate is None:
			return None
		best = max(((i, j) for i in range(2) for j in range(2)), key=lambda t: abs(matrix[t[0]][t[1]]))
		gate_unitary = gate._unitary_()
		return gate, matrix[best[0]][best[1]] / gate_unitary[best[0]][best[1]]

	def pauli_tuple(self, pauli: Pauli) -> tuple[Pauli, bool]:
		"""Return the image ``(Pauli, flip)`` of ``pauli`` under this gate."""
		if pauli is Pauli.X:
			return self._x_to
		if pauli is Pauli.Z:
			return self._z_to
		to_x, flip_x = self._x_to
		to_z, flip_z = self._z_to
		product, phase = _product(to_x.value, to_z.value)
		sign = 1j * (-1 if flip_x else 1) * (-1 if flip_z else 1) * phase
		return Pauli(product), bool(sign.real < 0)

	def decompose_rotation(self) -> list[tuple[Pauli, int]]:
		"""Decompose into a sequence of Pauli rotations ``(axis, quarter_turns)``."""
		x_rot = self.pauli_tuple(Pauli.X)
		y_rot = self.pauli_tuple(Pauli.Y)
		z_rot = self.pauli_tuple(Pauli.Z)
		whole = (x_rot[0] is Pauli.X, y_rot[0] is Pauli.Y, z_rot[0] is Pauli.Z)
		num_whole = sum(whole)
		flips = (x_rot[1], y_rot[1], z_rot[1])
		if num_whole == 3:
			if not any(flips):
				return []
			index = flips.index(False)
			return [(_index_to_pauli(index), 2)]
		if num_whole == 1:
			index = whole.index(True)
			pauli = _index_to_pauli(index)
			next_pauli = _index_to_pauli(index + 1)
			output: list[tuple[Pauli, int]] = []
			if flips[index]:
				output.append((next_pauli, 2))
			if self.pauli_tuple(next_pauli)[1]:
				output.append((pauli, -1))
			else:
				output.append((pauli, 1))
			return output
		if x_rot[0] is Pauli.Y:
			return [(Pauli.X, -1 if y_rot[1] else 1), (Pauli.Z, -1 if x_rot[1] else 1)]
		return [(Pauli.Z, 1 if y_rot[1] else -1), (Pauli.X, 1 if z_rot[1] else -1)]

	def decompose_gate(self) -> list[Gate]:
		"""Decompose into a sequence of ``H`` and Pauli-rotation gates."""
		if self == SingleQubitCliffordGate.H:
			return [H]
		return [_PAULI_GATE[pauli] ** (quarter_turns / 2) for pauli, quarter_turns in self.decompose_rotation()]

	def merged_with(self, second: SingleQubitCliffordGate) -> SingleQubitCliffordGate:
		"""Return the gate equal to ``self`` followed by ``second`` (up to phase)."""
		product = matmul(second._unitary_(), self._unitary_())
		return SingleQubitCliffordGate.from_xz_map(_pauli_image(product, Pauli.X), _pauli_image(product, Pauli.Z))

	def equivalent_gate_before(self, after: SingleQubitCliffordGate) -> SingleQubitCliffordGate:
		"""Return the gate making ``--output--self--`` equal to ``--self--gate--``."""
		return self.merged_with(after).merged_with(self**-1)

	def commutes_with_pauli(self, pauli: Pauli) -> bool:
		"""Whether the gate leaves the Pauli axis unchanged."""
		to, flip = self.pauli_tuple(pauli)
		return to is pauli and not flip

	def commutes_with_single_qubit_gate(self, gate: SingleQubitCliffordGate) -> bool:
		"""Whether the two Cliffords commute."""
		return self.merged_with(gate) == gate.merged_with(self)

	def _unitary_(self) -> Matrix:
		return [list(row) for row in self._unitary_matrix]

	def _has_unitary_(self) -> bool:
		return True

	def __pow__(self, power: float) -> SingleQubitCliffordGate:
		if int(power) == power:
			# The inverse is the adjoint unitary, which is a Clifford gate as
			# well; taking ``self ** -1`` here would recurse forever.
			base = self
			if power < 0:
				inverse = SingleQubitCliffordGate.from_unitary(_conj_t(self._unitary_()))
				assert inverse is not None  # the adjoint of a Clifford is one
				base = inverse
			# Every element of the group has order 1, 2, 3 or 4, so twelve
			# turns always come back to the identity; reducing modulo four
			# would be wrong for the elements of order three.
			result = SingleQubitCliffordGate.I
			for _ in range(abs(int(power)) % 12):
				result = result.merged_with(base)
			return result
		if int(2 * power) == 2 * power:
			half = self._square_root()
			if half is not None:
				# ``half`` is one of the group's elements, so the integral power
				# applies; the cast keeps the type checker happy about it.
				return half ** int(2 * power)
		raise NotImplementedError(f"cannot raise {self!r} to the power {power}")

	def _square_root(self) -> SingleQubitCliffordGate | None:
		for entry in _ALL_SINGLE_QUBIT_CLIFFORDS:
			candidate = SingleQubitCliffordGate(x_to=entry[0:2], z_to=entry[2:4], unitary=entry[4])
			if candidate.merged_with(candidate) == self:
				return candidate
		return None

	def _value_equality_values_(self) -> Any:
		return (self._x_to, self._z_to)

	def __str__(self) -> str:
		return "SingleQubitCliffordGate(...)"

	def __repr__(self) -> str:
		return f"vqubit.SingleQubitCliffordGate.from_xz_map({self._x_to!r}, {self._z_to!r})"


def _index_to_pauli(index: int) -> Pauli:
	return (Pauli.X, Pauli.Y, Pauli.Z)[index % 3]


def _relative_pauli(pauli: Pauli, offset: int) -> Pauli:
	order = [Pauli.X, Pauli.Y, Pauli.Z]
	return order[(order.index(pauli) + offset) % 3]


def _third_pauli(first: Pauli, second: Pauli) -> Pauli:
	order = [Pauli.X, Pauli.Y, Pauli.Z]
	return next(p for p in order if p is not first and p is not second)


def _merge_maps(
	pauli_map_to: dict[Pauli, tuple[Pauli, bool]] | None,
	x_to: tuple[Pauli, bool] | None,
	y_to: tuple[Pauli, bool] | None,
	z_to: tuple[Pauli, bool] | None,
) -> dict[Pauli, tuple[Pauli, bool]]:
	if pauli_map_to is not None:
		return dict(pauli_map_to)
	mapping: dict[Pauli, tuple[Pauli, bool]] = {}
	for pauli, transform in ((Pauli.X, x_to), (Pauli.Y, y_to), (Pauli.Z, z_to)):
		if transform is not None:
			mapping[pauli] = transform
	return mapping


def _initialize_class_attributes() -> None:
	by_map = {(entry[0], entry[1], entry[2], entry[3]): entry[4] for entry in _ALL_SINGLE_QUBIT_CLIFFORDS}

	def make(x_to: tuple[Pauli, bool], z_to: tuple[Pauli, bool]) -> SingleQubitCliffordGate:
		return SingleQubitCliffordGate.from_xz_map(x_to, z_to)

	SingleQubitCliffordGate.I = make((Pauli.X, False), (Pauli.Z, False))
	SingleQubitCliffordGate.X = make((Pauli.X, False), (Pauli.Z, True))
	SingleQubitCliffordGate.Y = make((Pauli.X, True), (Pauli.Z, True))
	SingleQubitCliffordGate.Z = make((Pauli.X, True), (Pauli.Z, False))
	SingleQubitCliffordGate.H = make((Pauli.Z, False), (Pauli.X, False))
	SingleQubitCliffordGate.S = make((Pauli.Y, False), (Pauli.Z, False))
	SingleQubitCliffordGate.X_sqrt = make((Pauli.X, False), (Pauli.Y, True))
	SingleQubitCliffordGate.Y_sqrt = make((Pauli.Z, True), (Pauli.X, False))
	SingleQubitCliffordGate.Z_sqrt = make((Pauli.Y, False), (Pauli.Z, False))
	SingleQubitCliffordGate.X_nsqrt = make((Pauli.X, False), (Pauli.Y, False))
	SingleQubitCliffordGate.Y_nsqrt = make((Pauli.Z, False), (Pauli.X, True))
	SingleQubitCliffordGate.Z_nsqrt = make((Pauli.Y, True), (Pauli.Z, False))
	SingleQubitCliffordGate.all_single_qubit_cliffords = tuple(
		make((entry[0], entry[1]), (entry[2], entry[3])) for entry in _ALL_SINGLE_QUBIT_CLIFFORDS
	)
	_ = by_map


_initialize_class_attributes()


@value_equality
class CliffordGate(Gate):
	"""A multi-qubit Clifford gate, stored by its unitary."""

	def __init__(self, *, unitary: Sequence[Sequence[complex]], num_qubits: int) -> None:
		self._unitary_matrix = to_nested_list(unitary)
		self._num_qubits = num_qubits

	@classmethod
	def from_unitary(cls, u: Sequence[Sequence[complex]]) -> CliffordGate:
		"""Return a Clifford gate from its unitary (assumed Clifford)."""
		matrix = to_nested_list(u)
		dim = len(matrix)
		return cls(unitary=matrix, num_qubits=dim.bit_length() - 1)

	@classmethod
	def from_op_list(cls, operations: Sequence[Any], qubit_order: Sequence[Qid]) -> CliffordGate:
		"""Compose the given operations into a Clifford gate."""
		from ._circuits import Circuit

		circuit = Circuit(list(operations))
		matrix = to_nested_list(circuit.unitary())
		return cls(unitary=matrix, num_qubits=len(qubit_order))

	def _num_qubits_(self) -> int:
		return self._num_qubits

	def _unitary_(self) -> Matrix:
		return [list(row) for row in self._unitary_matrix]

	def _has_stabilizer_effect_(self) -> bool:
		return True

	def __pow__(self, power: float) -> CliffordGate:
		if int(power) != power:
			raise NotImplementedError("CliffordGate only supports integer powers")
		dim = 1 << self._num_qubits
		result = identity(dim)
		base = self._unitary_matrix if power >= 0 else _invert(self._unitary_matrix)
		for _ in range(abs(int(power))):
			result = matmul(base, result)
		return CliffordGate(unitary=result, num_qubits=self._num_qubits)

	def _value_equality_values_(self) -> Any:
		# A Clifford gate is defined up to a global phase, and a matrix computed
		# by multiplying two others carries floating-point noise (``CZSWAP ** 3``
		# against ``CZSWAP``). Normalise the phase away and round, so equal gates
		# compare equal, which is what comparing Cirq's tableaux amounts to.
		flat = [complex(entry) for row in self._unitary_matrix for entry in row]
		pivot = next((entry for entry in flat if abs(entry) > 1e-9), 1 + 0j)
		phase = pivot / abs(pivot)
		return (
			tuple((round((entry / phase).real, 9), round((entry / phase).imag, 9)) for entry in flat),
			self._num_qubits,
		)

	def __str__(self) -> str:
		return f"CliffordGate({self._num_qubits}q)"

	def __repr__(self) -> str:
		return f"vqubit.CliffordGate.from_unitary({self._unitary_matrix!r})"


def _invert(matrix: Matrix) -> Matrix:
	return _conj_t(matrix)


_Q0, _Q1 = LineQubit.range(2)
CXSWAP = CliffordGate.from_op_list([CNOT(_Q0, _Q1), SWAP(_Q0, _Q1)], [_Q0, _Q1])
CZSWAP = CliffordGate.from_op_list([CZ(_Q0, _Q1), SWAP(_Q0, _Q1)], [_Q0, _Q1])
