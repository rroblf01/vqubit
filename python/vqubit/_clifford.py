"""A Clifford (stabilizer) simulator.

Clifford circuits -- those built from Hadamard, phase and CNOT-like gates --
never leave the stabilizer states, and a stabilizer state is described by a
small binary tableau (Aaronson and Gottesman, *Improved Simulation of
Stabilizer Circuits*). Tracking that tableau instead of the state vector is why
this simulator can reach many more qubits than :class:`~vqubit.Simulator`.

The tableau and the simulation loop live in V: Python builds a flat instruction
stream and reshapes the results, so the whole run happens in one native call.
"""

from __future__ import annotations

import random
from collections.abc import Sequence
from typing import Any

import vqubit_native as _native

from ._circuits import Circuit
from ._compat import HAS_NUMPY, _np
from ._gates import (
	CNOTGate,
	CZGate,
	HPowGate,
	IdentityGate,
	MeasurementGate,
	SWAPGate,
	XPowGate,
	YPowGate,
	ZPowGate,
)
from ._linalg import vector_from_flat
from ._pauli import Pauli
from ._qubits import LineQubit, Qid

_INSTR_H = 0
_INSTR_S = 1
_INSTR_SDG = 2
_INSTR_X = 3
_INSTR_Y = 4
_INSTR_Z = 5
_INSTR_CNOT = 6
_INSTR_CZ = 7
_INSTR_SWAP = 8
_INSTR_MEASURE = 9


def _instruction_code(gate: Any) -> int | None:
	"""Return the instruction code for a Clifford gate, or ``None`` for identity."""
	if isinstance(gate, IdentityGate):
		return None
	if isinstance(gate, HPowGate) and gate.exponent == 1:
		return _INSTR_H
	if isinstance(gate, ZPowGate):
		if gate.exponent == 1:
			return _INSTR_Z
		if gate.exponent == 0.5:
			return _INSTR_S
		if gate.exponent == -0.5:
			return _INSTR_SDG
	if isinstance(gate, XPowGate) and gate.exponent in (1, -1):
		return _INSTR_X
	if isinstance(gate, YPowGate) and gate.exponent in (1, -1):
		return _INSTR_Y
	if isinstance(gate, CNOTGate):
		return _INSTR_CNOT
	if isinstance(gate, CZGate):
		return _INSTR_CZ
	if isinstance(gate, SWAPGate):
		return _INSTR_SWAP
	raise TypeError(f"{gate!r} is not a Clifford gate")


def _build_instructions(
	program: Circuit,
	index: dict[Qid, int],
	allow_measurements: bool,
) -> tuple[list[int], list[tuple[str, int, int]], dict[str, int]]:
	"""Return ``(instructions, measure_meta, key_sizes)`` for a circuit."""
	instructions: list[int] = []
	measure_meta: list[tuple[str, int, int]] = []
	key_sizes: dict[str, int] = {}
	for operation in program.all_operations():
		gate = operation.gate
		if isinstance(gate, MeasurementGate):
			if not allow_measurements:
				raise ValueError("simulate() cannot run a circuit with measurements; use run()")
			invert_mask = gate.full_invert_mask
			for position, qubit in enumerate(operation.qubits):
				instructions.extend((_INSTR_MEASURE, index[qubit], 0))
				measure_meta.append((gate.key, position, int(invert_mask[position])))
				key_sizes[gate.key] = max(key_sizes.get(gate.key, 0), position + 1)
			continue
		code = _instruction_code(gate)
		if code is None:
			continue
		targets = [index[q] for q in operation.qubits]
		second = targets[1] if len(targets) > 1 else 0
		instructions.extend((code, targets[0], second))
	return instructions, measure_meta, key_sizes


def _bool_matrix(rows: Sequence[Sequence[Any]]) -> Any:
	if HAS_NUMPY:
		return _np.array(rows, dtype=bool)
	return [[bool(value) for value in row] for row in rows]


def _bool_vector(values: Sequence[Any]) -> Any:
	if HAS_NUMPY:
		return _np.array(values, dtype=bool)
	return [bool(value) for value in values]


def _to_int_matrix(value: Any) -> list[list[int]]:
	if HAS_NUMPY and isinstance(value, _np.ndarray):
		return [[int(bool(v)) for v in row] for row in value.tolist()]
	return [[int(bool(v)) for v in row] for row in value]


def _to_int_vector(value: Any) -> list[int]:
	if HAS_NUMPY and isinstance(value, _np.ndarray):
		return [int(bool(v)) for v in value.tolist()]
	return [int(bool(v)) for v in value]


def _matmul(a: Sequence[Sequence[int]], b: Sequence[Sequence[int]]) -> list[list[int]]:
	rows = len(a)
	inner = len(b)
	cols = len(b[0]) if inner else 0
	return [[sum(a[i][k] * b[k][j] for k in range(inner)) for j in range(cols)] for i in range(rows)]


def _matvec(a: Sequence[Sequence[int]], v: Sequence[int]) -> list[int]:
	return [sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a))]


def _transpose(a: Sequence[Sequence[int]]) -> list[list[int]]:
	return [list(col) for col in zip(*a, strict=False)]


def _mod2(matrix: Sequence[Sequence[int]]) -> list[list[int]]:
	return [[value % 2 for value in row] for row in matrix]


def _mod2_vector(values: Sequence[int]) -> list[int]:
	return [value % 2 for value in values]


def _vec_add(a: Sequence[int], b: Sequence[int]) -> list[int]:
	return [a[i] + b[i] for i in range(len(a))]


def _tril_add(a: Sequence[Sequence[int]], b: Sequence[Sequence[int]]) -> list[list[int]]:
	size = len(a)
	return [[0 if i <= j else a[i][j] + b[i][j] for j in range(size)] for i in range(size)]


def _resolve_random(seed: Any) -> Any:
	if seed is None or isinstance(seed, int):
		return random.Random(seed)
	return seed


def _draw_bit(prng: Any) -> int:
	try:
		return int(prng.randint(2))
	except TypeError:
		return int(prng.randint(0, 1))


class CliffordTableau:
	"""Tableau representation of a stabilizer state.

	The tableau stores the stabilizer generators of the state using three binary
	arrays: ``xs``, ``zs`` and ``rs``. Each row is a Pauli string ``P`` that is an
	eigenoperator of the state vector with eigenvalue one: ``P|psi> = |psi>``.

	References:
		Aaronson and Gottesman, https://arxiv.org/abs/quant-ph/0406196
	"""

	def __init__(
		self,
		num_qubits: int,
		initial_state: int = 0,
		rs: Any | None = None,
		xs: Any | None = None,
		zs: Any | None = None,
	) -> None:
		self.n = num_qubits
		self.initial_state = initial_state
		self._rs = self._reconstruct_rs(rs)
		self._xs = self._reconstruct_xs(xs)
		self._zs = self._reconstruct_zs(zs)

	def _reconstruct_rs(self, rs: Any | None) -> list[int]:
		size = 2 * self.n
		if rs is None:
			new_rs = [0] * (size + 1)
			for i in range(self.n):
				new_rs[self.n + i] = (self.initial_state >> (self.n - 1 - i)) & 1
			return new_rs
		values = _to_int_vector(rs)
		if len(values) == size:
			return [*values, 0]
		raise ValueError("The value you passed for rs is not the correct shape.")

	def _reconstruct_xs(self, xs: Any | None) -> list[list[int]]:
		if xs is None:
			new_xs = [[0] * self.n for _ in range(2 * self.n + 1)]
			for i in range(self.n):
				new_xs[i][i] = 1
			return new_xs
		rows = _to_int_matrix(xs)
		if len(rows) == 2 * self.n and all(len(row) == self.n for row in rows):
			return [*rows, [0] * self.n]
		raise ValueError("The value you passed for xs is not the correct shape.")

	def _reconstruct_zs(self, zs: Any | None) -> list[list[int]]:
		if zs is None:
			new_zs = [[0] * self.n for _ in range(2 * self.n + 1)]
			for i in range(self.n):
				new_zs[self.n + i][i] = 1
			return new_zs
		rows = _to_int_matrix(zs)
		if len(rows) == 2 * self.n and all(len(row) == self.n for row in rows):
			return [*rows, [0] * self.n]
		raise ValueError("The value you passed for zs is not the correct shape.")

	@classmethod
	def from_flat(cls, flat: Sequence[int], num_qubits: int) -> CliffordTableau:
		"""Rebuild a tableau from the flattened native form."""
		n = num_qubits
		side = 2 * n + 1
		xs = [[int(flat[i * side + j]) for j in range(n)] for i in range(2 * n)]
		zs = [[int(flat[i * side + n + j]) for j in range(n)] for i in range(2 * n)]
		rs = [int(flat[i * side + 2 * n]) for i in range(2 * n)]
		return cls(n, rs=rs, xs=xs, zs=zs)

	@property
	def xs(self) -> Any:
		return _bool_matrix(self._xs[:-1])

	@xs.setter
	def xs(self, new_xs: Any) -> None:
		rows = _to_int_matrix(new_xs)
		assert len(rows) == 2 * self.n and all(len(row) == self.n for row in rows)
		self._xs[:-1] = rows

	@property
	def zs(self) -> Any:
		return _bool_matrix(self._zs[:-1])

	@zs.setter
	def zs(self, new_zs: Any) -> None:
		rows = _to_int_matrix(new_zs)
		assert len(rows) == 2 * self.n and all(len(row) == self.n for row in rows)
		self._zs[:-1] = rows

	@property
	def rs(self) -> Any:
		return _bool_vector(self._rs[:-1])

	@rs.setter
	def rs(self, new_rs: Any) -> None:
		values = _to_int_vector(new_rs)
		assert len(values) == 2 * self.n
		self._rs[:-1] = values

	def matrix(self) -> Any:
		"""Return the ``2n`` by ``2n`` matrix representation of the tableau."""
		rows = [self._xs[i] + self._zs[i] for i in range(2 * self.n)]
		return _bool_matrix(rows)

	def _validate(self) -> bool:
		"""Check whether the tableau satisfies the symplectic property."""
		size = 2 * self.n
		table = [self._xs[i] + self._zs[i] for i in range(size)]
		perm = list(range(self.n, size)) + list(range(self.n))
		skew_eye = [[1 if j == perm[i] else 0 for j in range(size)] for i in range(size)]
		product = _mod2(_matmul(_matmul(_transpose(table), skew_eye), table))
		return product == skew_eye

	def __eq__(self, other: Any) -> bool:
		if not isinstance(other, type(self)):
			return NotImplemented
		return self.n == other.n and self._rs == other._rs and self._xs == other._xs and self._zs == other._zs

	def __hash__(self) -> int:
		return hash(
			(
				self.n,
				tuple(tuple(row) for row in self._xs),
				tuple(tuple(row) for row in self._zs),
				tuple(self._rs),
			)
		)

	def __copy__(self) -> CliffordTableau:
		return self.copy()

	def copy(self, deep_copy_buffers: bool = True) -> CliffordTableau:
		"""Return a deep copy of this tableau."""
		state = CliffordTableau(self.n)
		state._xs = [list(row) for row in self._xs]
		state._zs = [list(row) for row in self._zs]
		state._rs = list(self._rs)
		return state

	def __repr__(self) -> str:
		return (
			f"vqubit.CliffordTableau({self.n}, "
			f"rs={self._rs[:-1]!r}, xs={self._xs[:-1]!r}, zs={self._zs[:-1]!r}, "
			f"initial_state={self.initial_state})"
		)

	def __str__(self) -> str:
		words = []
		for i in range(self.n, 2 * self.n):
			words.append("- " if self._rs[i] else "+ ")
			for k in range(self.n):
				if self._xs[i][k] and not self._zs[i][k]:
					words.append("X ")
				elif not self._xs[i][k] and self._zs[i][k]:
					words.append("Z ")
				elif self._xs[i][k] and self._zs[i][k]:
					words.append("Y ")
				else:
					words.append("I ")
			words.append("\n")
		return "".join(words[:-1])

	def then(self, second: CliffordTableau) -> CliffordTableau:
		"""Return the composed tableau ``this then second``, up to global phase."""
		if not isinstance(second, CliffordTableau):
			raise TypeError("The type for second tableau must be the CliffordTableau type")
		if self.n != second.n:
			raise ValueError(f"Mismatched number of qubits of two tableaux: {self.n} vs {second.n}.")

		n = self.n
		size = 2 * n
		m1 = [self._xs[i] + self._zs[i] for i in range(size)]
		m2 = [second._xs[i] + second._zs[i] for i in range(size)]

		num_ys1 = [sum(m1[i][j] * m1[i][n + j] for j in range(n)) for i in range(size)]
		num_ys2 = [sum(m2[i][j] * m2[i][n + j] for j in range(n)) for i in range(size)]
		p1 = [value % 2 for value in num_ys1]
		p2 = [value % 2 for value in num_ys2]

		s1 = [self._rs[i] + (num_ys1[i] % 4) // 2 for i in range(size)]
		s2 = [second._rs[i] + (num_ys2[i] % 4) // 2 for i in range(size)]

		lmbda = [[0] * size for _ in range(size)]
		for i in range(n):
			lmbda[i][n + i] = 1

		m_12 = _mod2(_matmul(m1, m2))
		p_12 = _mod2_vector(_vec_add(p1, _matvec(m1, p2)))

		outer = [[p2[i] * p2[j] for j in range(size)] for i in range(size)]
		m2_lam_m2t = _matmul(_matmul(m2, lmbda), _transpose(m2))
		inner = _tril_add(outer, m2_lam_m2t)
		m1_inner = _matmul(_matmul(m1, inner), _transpose(m1))
		diagonal = [m1_inner[i][i] for i in range(size)]

		m1_s2 = _matvec(m1, s2)
		m1_p2 = _matvec(m1, p2)
		s_12 = [s1[i] + m1_s2[i] + p1[i] * m1_p2[i] + diagonal[i] for i in range(size)]

		num_ys12 = [sum(m_12[i][j] * m_12[i][n + j] for j in range(n)) for i in range(size)]
		merged_sign = [((p_12[i] + 2 * s_12[i] - num_ys12[i]) % 4) // 2 for i in range(size)]

		merged = CliffordTableau(n)
		merged.xs = [row[:n] for row in m_12]
		merged.zs = [row[n:] for row in m_12]
		merged.rs = merged_sign
		return merged

	def inverse(self) -> CliffordTableau:
		"""Return the inverse tableau of this tableau."""
		n = self.n
		size = 2 * n
		xs = self._xs[:size]
		zs = self._zs[:size]
		ret = CliffordTableau(n)
		new_xs = [None] * size
		new_zs = [None] * size
		top_zs = _transpose(zs[n:])
		bottom_zs = _transpose(zs[:n])
		bottom_xs = _transpose(xs[n:])
		top_xs = _transpose(xs[:n])
		for i in range(n):
			new_xs[i] = top_zs[i]
			new_zs[i] = bottom_zs[i]
			new_xs[n + i] = bottom_xs[i]
			new_zs[n + i] = top_xs[i]
		ret.xs = new_xs
		ret.zs = new_zs
		ret.rs = ret.then(self).rs
		return ret

	def __matmul__(self, second: CliffordTableau) -> CliffordTableau:
		if not isinstance(second, CliffordTableau):
			return NotImplemented
		return second.then(self)

	def _rowsum(self, q1: int, q2: int) -> None:
		"""Multiply the stabilizer in row ``q1`` by the one in row ``q2``."""

		def g(x1: int, z1: int, x2: int, z2: int) -> int:
			if not x1 and not z1:
				return 0
			if x1 and z1:
				return int(z2) - int(x2)
			if x1 and not z1:
				return int(z2) * (2 * int(x2) - 1)
			return int(x2) * (1 - 2 * int(z2))

		total = 2 * int(self._rs[q1]) + 2 * int(self._rs[q2])
		for j in range(self.n):
			total += g(self._xs[q2][j], self._zs[q2][j], self._xs[q1][j], self._zs[q1][j])
		self._rs[q1] = bool(total % 4)
		for j in range(self.n):
			self._xs[q1][j] ^= self._xs[q2][j]
			self._zs[q1][j] ^= self._zs[q2][j]

	def _row_to_dense_pauli(self, i: int) -> Any:
		from ._pauli import DensePauliString

		coefficient = -1 if self._rs[i] else 1
		mask = []
		for k in range(self.n):
			if self._xs[i][k] and not self._zs[i][k]:
				mask.append("X")
			elif not self._xs[i][k] and self._zs[i][k]:
				mask.append("Z")
			elif self._xs[i][k] and self._zs[i][k]:
				mask.append("Y")
			else:
				mask.append("I")
		return DensePauliString("".join(mask), coefficient=coefficient)

	def stabilizers(self) -> list[Any]:
		"""Return the stabilizer generators ``{S_i}`` with ``S_i|psi> = |psi>``."""
		return [self._row_to_dense_pauli(i) for i in range(self.n, 2 * self.n)]

	def destabilizers(self) -> list[Any]:
		"""Return the destabilizer generators, which complete the Pauli group."""
		return [self._row_to_dense_pauli(i) for i in range(self.n)]

	def stabilizer_pauli_strings(self) -> list[Any]:
		"""Return the stabilizer generators as sparse Pauli strings."""
		from ._pauli import PauliString

		result: list[Any] = []
		for k in range(self.n):
			row = self.n + k
			qubit_pauli_map: dict[Qid, Pauli] = {}
			for j in range(self.n):
				pauli = _PAULI_BITS[(self._xs[row][j], self._zs[row][j])]
				if pauli is not Pauli.I:
					qubit_pauli_map[LineQubit(j)] = pauli
			sign = -1 if self._rs[row] else 1
			result.append(PauliString(qubit_pauli_map, coefficient=sign))
		return result

	def _measure(self, q: int, prng: Any) -> int:
		"""Measure the ``q``-th qubit in the Z basis, updating the tableau."""
		is_commuting = True
		p = -1
		for i in range(self.n, 2 * self.n):
			if self._xs[i][q]:
				p = i
				is_commuting = False
				break

		if is_commuting:
			scratch = 2 * self.n
			for j in range(self.n):
				self._xs[scratch][j] = 0
				self._zs[scratch][j] = 0
			self._rs[scratch] = False
			for i in range(self.n):
				if self._xs[i][q]:
					self._rowsum(scratch, self.n + i)
			return int(self._rs[scratch])

		for i in range(2 * self.n):
			if i != p and self._xs[i][q]:
				self._rowsum(i, p)

		for j in range(self.n):
			self._xs[p - self.n][j] = self._xs[p][j]
			self._zs[p - self.n][j] = self._zs[p][j]
		self._rs[p - self.n] = self._rs[p]

		for j in range(self.n):
			self._xs[p][j] = 0
			self._zs[p][j] = 0

		self._zs[p][q] = 1
		self._rs[p] = bool(_draw_bit(prng))
		return int(self._rs[p])

	def apply_x(self, axis: int, exponent: float = 1, global_shift: float = 0) -> None:
		if exponent % 2 == 0:
			return
		if exponent % 0.5 != 0.0:
			raise ValueError("X exponent must be half integer")
		effective_exponent = exponent % 2
		if effective_exponent == 0.5:
			for i in range(2 * self.n):
				self._xs[i][axis] ^= self._zs[i][axis]
			for i in range(2 * self.n):
				self._rs[i] ^= self._xs[i][axis] & self._zs[i][axis]
		elif effective_exponent == 1:
			for i in range(2 * self.n):
				self._rs[i] ^= self._zs[i][axis]
		elif effective_exponent == 1.5:
			for i in range(2 * self.n):
				self._rs[i] ^= self._xs[i][axis] & self._zs[i][axis]
			for i in range(2 * self.n):
				self._xs[i][axis] ^= self._zs[i][axis]

	def apply_y(self, axis: int, exponent: float = 1, global_shift: float = 0) -> None:
		if exponent % 2 == 0:
			return
		if exponent % 0.5 != 0.0:
			raise ValueError("Y exponent must be half integer")
		effective_exponent = exponent % 2
		if effective_exponent == 0.5:
			for i in range(2 * self.n):
				self._rs[i] ^= self._xs[i][axis] & (1 - self._zs[i][axis])
			for i in range(2 * self.n):
				self._xs[i][axis], self._zs[i][axis] = self._zs[i][axis], self._xs[i][axis]
		elif effective_exponent == 1:
			for i in range(2 * self.n):
				self._rs[i] ^= self._xs[i][axis] ^ self._zs[i][axis]
		elif effective_exponent == 1.5:
			for i in range(2 * self.n):
				self._rs[i] ^= (1 - self._xs[i][axis]) & self._zs[i][axis]
			for i in range(2 * self.n):
				self._xs[i][axis], self._zs[i][axis] = self._zs[i][axis], self._xs[i][axis]

	def apply_z(self, axis: int, exponent: float = 1, global_shift: float = 0) -> None:
		if exponent % 2 == 0:
			return
		if exponent % 0.5 != 0.0:
			raise ValueError("Z exponent must be half integer")
		effective_exponent = exponent % 2
		if effective_exponent == 0.5:
			for i in range(2 * self.n):
				self._rs[i] ^= self._xs[i][axis] & self._zs[i][axis]
			for i in range(2 * self.n):
				self._zs[i][axis] ^= self._xs[i][axis]
		elif effective_exponent == 1:
			for i in range(2 * self.n):
				self._rs[i] ^= self._xs[i][axis]
		elif effective_exponent == 1.5:
			for i in range(2 * self.n):
				self._rs[i] ^= self._xs[i][axis] & (1 - self._zs[i][axis])
			for i in range(2 * self.n):
				self._zs[i][axis] ^= self._xs[i][axis]

	def apply_h(self, axis: int, exponent: float = 1, global_shift: float = 0) -> None:
		if exponent % 2 == 0:
			return
		if exponent % 1 != 0:
			raise ValueError("H exponent must be integer")
		self.apply_y(axis, 0.5)
		self.apply_x(axis)

	def apply_cz(self, control_axis: int, target_axis: int, exponent: float = 1, global_shift: float = 0) -> None:
		if exponent % 2 == 0:
			return
		if exponent % 1 != 0:
			raise ValueError("CZ exponent must be integer")
		for i in range(2 * self.n):
			self._xs[i][target_axis], self._zs[i][target_axis] = (
				self._zs[i][target_axis],
				self._xs[i][target_axis],
			)
		for i in range(2 * self.n):
			self._rs[i] ^= self._xs[i][target_axis] & self._zs[i][target_axis]
		for i in range(2 * self.n):
			self._rs[i] ^= (
				self._xs[i][control_axis]
				& self._zs[i][target_axis]
				& (1 - (self._xs[i][target_axis] ^ self._zs[i][control_axis]))
			)
		for i in range(2 * self.n):
			self._xs[i][target_axis] ^= self._xs[i][control_axis]
		for i in range(2 * self.n):
			self._zs[i][control_axis] ^= self._zs[i][target_axis]
		for i in range(2 * self.n):
			self._xs[i][target_axis], self._zs[i][target_axis] = (
				self._zs[i][target_axis],
				self._xs[i][target_axis],
			)
		for i in range(2 * self.n):
			self._rs[i] ^= self._xs[i][target_axis] & self._zs[i][target_axis]

	def apply_cx(self, control_axis: int, target_axis: int, exponent: float = 1, global_shift: float = 0) -> None:
		if exponent % 2 == 0:
			return
		if exponent % 1 != 0:
			raise ValueError("CX exponent must be integer")
		for i in range(2 * self.n):
			self._rs[i] ^= (
				self._xs[i][control_axis]
				& self._zs[i][target_axis]
				& (1 - (self._xs[i][target_axis] ^ self._zs[i][control_axis]))
			)
		for i in range(2 * self.n):
			self._xs[i][target_axis] ^= self._xs[i][control_axis]
		for i in range(2 * self.n):
			self._zs[i][control_axis] ^= self._zs[i][target_axis]

	def apply_global_phase(self, coefficient: Any) -> None:
		pass

	def measure(self, axes: Sequence[int], seed: Any = None) -> list[int]:
		"""Measure each axis in ``axes``, returning the outcome in order."""
		prng = _resolve_random(seed)
		return [self._measure(axis, prng) for axis in axes]

	def sample(self, axes: Sequence[int], repetitions: int = 1, seed: Any = None) -> Any:
		"""Sample ``axes`` ``repetitions`` times, without mutating this tableau."""
		prng = _resolve_random(seed)
		measurements = []
		for _ in range(repetitions):
			state = self.copy()
			measurements.append(state.measure(axes, prng))
		if HAS_NUMPY:
			return _np.array(measurements, dtype=_np.uint8)
		return measurements

	def kron(self, other: Any) -> CliffordTableau:
		raise NotImplementedError()

	def factor(self, axes: Sequence[int], *, validate: bool = True, atol: float = 1e-7) -> Any:
		raise NotImplementedError()

	def reindex(self, axes: Sequence[int]) -> CliffordTableau:
		raise NotImplementedError()

	@property
	def supports_factor(self) -> bool:
		return False

	@property
	def can_represent_mixed_states(self) -> bool:
		return False


_PAULI_BITS = {
	(0, 0): Pauli.I,
	(1, 0): Pauli.X,
	(1, 1): Pauli.Y,
	(0, 1): Pauli.Z,
}


class CliffordState:
	"""A stabilizer state, with the qubits it lives on."""

	def __init__(self, tableau: Any, qubits: Sequence[Qid]) -> None:
		self._tableau = tableau
		self._qubits = tuple(qubits)

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits

	def state_vector(self) -> list[complex]:
		return vector_from_flat(self._tableau.state_vector())

	def stabilizers(self) -> list[Any]:
		flat = self._tableau.stabilizer_tableau()
		tableau = CliffordTableau.from_flat(flat, len(self._qubits))
		return tableau.stabilizer_pauli_strings()

	def dirac_notation(self, decimals: int = 3) -> str:
		from ._sim import StateVector

		return StateVector(self.state_vector(), self._qubits).dirac_notation(decimals)

	def __repr__(self) -> str:
		return f"vqubit.CliffordState({len(self._qubits)} qubits)"


class CliffordSimulationResult:
	"""The result of :meth:`CliffordSimulator.simulate`."""

	def __init__(self, state: CliffordState) -> None:
		self._final_state = state

	@property
	def final_state(self) -> CliffordState:
		return self._final_state

	def state_vector(self) -> list[complex]:
		return self._final_state.state_vector()

	def __repr__(self) -> str:
		return f"vqubit.CliffordSimulationResult({self._final_state!r})"


class CliffordSimulator:
	"""A stabilizer simulator for Clifford circuits."""

	def __init__(self, *, seed: int | None = None) -> None:
		self._seed = seed

	def simulate(
		self,
		program: Circuit,
		*,
		qubit_order: Sequence[Qid] | None = None,
	) -> CliffordSimulationResult:
		qubits = self._resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		tableau = _native.NativeTableau(len(qubits))
		for operation in program.all_operations():
			gate = operation.gate
			if isinstance(gate, MeasurementGate):
				raise ValueError("simulate() cannot run a circuit with measurements; use run()")
			code = _instruction_code(gate)
			if code is None:
				continue
			targets = [index[q] for q in operation.qubits]
			tableau.apply(code, targets[0], targets[1] if len(targets) > 1 else 0)
		return CliffordSimulationResult(CliffordState(tableau, qubits))

	def run(
		self,
		program: Circuit,
		repetitions: int = 1,
		*,
		qubit_order: Sequence[Qid] | None = None,
		seed: int | None = None,
	) -> Any:
		from ._sim import Result

		if repetitions < 1:
			raise ValueError("repetitions must be at least 1")
		qubits = self._resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		instructions, measure_meta, key_sizes = _build_instructions(program, index, allow_measurements=True)
		effective_seed = seed if seed is not None else self._seed
		if effective_seed is None:
			effective_seed = random.randrange(1 << 31)
		invert_mask = [meta[2] for meta in measure_meta]
		bits = _native.clifford_run(instructions, len(qubits), invert_mask, repetitions, int(effective_seed))

		records: dict[str, list[list[int]]] = {
			key: [[0] * size for _ in range(repetitions)] for key, size in key_sizes.items()
		}
		num_measurements = len(measure_meta)
		for rep in range(repetitions):
			for m, (key, position, _) in enumerate(measure_meta):
				records[key][rep][position] = bits[rep * num_measurements + m]
		return Result(records)

	sample = run

	@staticmethod
	def _resolve_qubits(program: Circuit, qubit_order: Sequence[Qid] | None) -> tuple[Qid, ...]:
		if qubit_order is None:
			return program.all_qubits()
		qubits = tuple(qubit_order)
		missing = set(program.all_qubits()) - set(qubits)
		if missing:
			raise ValueError(f"qubit_order is missing qubits used by the circuit: {sorted(missing)!r}")
		return qubits
