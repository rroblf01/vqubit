"""A Clifford (stabilizer) simulator.

Clifford circuits -- those built from Hadamard, phase and CNOT-like gates --
never leave the stabilizer states, and a stabilizer state is described by a
small binary tableau (Aaronson and Gottesman, *Improved Simulation of
Stabilizer Circuits*). Tracking that tableau instead of the state vector is why
this simulator can reach many more qubits than :class:`~vqubit.Simulator`.

The tableau has two blocks of rows: the first ``n`` are *destabilizers*, the
last ``n`` are the *stabilizer* generators of the state. Each row is a Pauli
operator stored as a pair of bit-vectors plus a phase bit.
"""

from __future__ import annotations

import random
from collections.abc import Sequence
from typing import Any

from ._circuits import Circuit
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
from ._linalg import identity, kron, matmul
from ._ops import Operation
from ._pauli import Pauli
from ._qubits import Qid

_PAULI_BITS = {
	(0, 0): Pauli.I,
	(1, 0): Pauli.X,
	(1, 1): Pauli.Y,
	(0, 1): Pauli.Z,
}


class CliffordTableau:
	"""The binary tableau of a stabilizer state."""

	def __init__(self, num_qubits: int) -> None:
		self.num_qubits = num_qubits
		self.rows = 2 * num_qubits + 1  # one scratch row
		self.x = [[0] * num_qubits for _ in range(self.rows)]
		self.z = [[0] * num_qubits for _ in range(self.rows)]
		self.r = [0] * self.rows
		for i in range(num_qubits):
			self.x[i][i] = 1
			self.z[num_qubits + i][i] = 1

	@staticmethod
	def _g(x1: int, z1: int, x2: int, z2: int) -> int:
		"""Return the exponent of ``i`` from multiplying two Paulis."""
		if x1 == 0 and z1 == 0:
			return 0
		if x1 == 1 and z1 == 1:
			return z2 - x2
		if x1 == 1 and z1 == 0:
			return z2 * (2 * x2 - 1)
		return x2 * (1 - 2 * z2)

	def rowsum(self, h: int, i: int) -> None:
		"""Set row ``h`` to the product of rows ``h`` and ``i``."""
		total = 2 * self.r[h] + 2 * self.r[i]
		for j in range(self.num_qubits):
			total += self._g(self.x[i][j], self.z[i][j], self.x[h][j], self.z[h][j])
		self.r[h] = 0 if total % 4 == 0 else 1
		for j in range(self.num_qubits):
			self.x[h][j] ^= self.x[i][j]
			self.z[h][j] ^= self.z[i][j]

	def hadamard(self, a: int) -> None:
		for i in range(2 * self.num_qubits):
			self.r[i] ^= self.x[i][a] * self.z[i][a]
			self.x[i][a], self.z[i][a] = self.z[i][a], self.x[i][a]

	def phase(self, a: int) -> None:
		for i in range(2 * self.num_qubits):
			self.r[i] ^= self.x[i][a] * self.z[i][a]
			self.z[i][a] ^= self.x[i][a]

	def cnot(self, a: int, b: int) -> None:
		for i in range(2 * self.num_qubits):
			self.r[i] ^= self.x[i][a] * self.z[i][b] * (self.x[i][b] ^ self.z[i][a] ^ 1)
			self.x[i][b] ^= self.x[i][a]
			self.z[i][a] ^= self.z[i][b]

	def pauli_x(self, a: int) -> None:
		for i in range(2 * self.num_qubits):
			self.r[i] ^= self.z[i][a]

	def pauli_z(self, a: int) -> None:
		for i in range(2 * self.num_qubits):
			self.r[i] ^= self.x[i][a]

	def pauli_y(self, a: int) -> None:
		self.pauli_x(a)
		self.pauli_z(a)

	def cz(self, a: int, b: int) -> None:
		self.hadamard(b)
		self.cnot(a, b)
		self.hadamard(b)

	def swap(self, a: int, b: int) -> None:
		self.cnot(a, b)
		self.cnot(b, a)
		self.cnot(a, b)

	def measure(self, a: int, rng: random.Random) -> int:
		"""Measure qubit ``a`` in the Z basis, updating the tableau."""
		n = self.num_qubits
		p = -1
		for i in range(n, 2 * n):
			if self.x[i][a]:
				p = i
				break
		if p != -1:
			for i in range(2 * n):
				if i != p and self.x[i][a]:
					self.rowsum(i, p)
			for j in range(n):
				self.x[p - n][j] = self.x[p][j]
				self.z[p - n][j] = self.z[p][j]
			self.r[p - n] = self.r[p]
			for j in range(n):
				self.x[p][j] = 0
				self.z[p][j] = 0
			self.z[p][a] = 1
			self.r[p] = rng.getrandbits(1)
			return self.r[p]

		scratch = 2 * n
		for j in range(n):
			self.x[scratch][j] = 0
			self.z[scratch][j] = 0
		self.r[scratch] = 0
		for i in range(n):
			if self.x[i][a]:
				self.rowsum(scratch, i + n)
		return self.r[scratch]

	def _stabilizer_matrix(self, row: int) -> list[list[complex]]:
		matrix: list[list[complex]] = [[1 + 0j]]
		for j in range(self.num_qubits):
			matrix = kron(matrix, _PAULI_BITS[(self.x[row][j], self.z[row][j])].matrix())
		sign = -1 if self.r[row] else 1
		return [[sign * entry for entry in row_values] for row_values in matrix]

	def stabilizer_pauli_strings(self) -> list[Any]:
		"""Return the stabilizer generators as Pauli strings."""
		from ._pauli import PauliString
		from ._qubits import LineQubit

		result: list[Any] = []
		for k in range(self.num_qubits):
			row = self.num_qubits + k
			qubit_pauli_map: dict[Qid, Pauli] = {}
			for j in range(self.num_qubits):
				pauli = _PAULI_BITS[(self.x[row][j], self.z[row][j])]
				if pauli is not Pauli.I:
					qubit_pauli_map[LineQubit(j)] = pauli
			result.append(PauliString(qubit_pauli_map, coefficient=-1 if self.r[row] else 1))
		return result

	def state_vector(self) -> list[complex]:
		"""Return the state vector, defined up to a global phase."""
		dim = 1 << self.num_qubits
		accumulator = identity(dim)
		for k in range(self.num_qubits):
			stabilizer = self._stabilizer_matrix(self.num_qubits + k)
			projector = [[stabilizer[i][j] + (1 if i == j else 0) for j in range(dim)] for i in range(dim)]
			accumulator = matmul(accumulator, projector)
		for column in range(dim):
			vector = [accumulator[i][column] for i in range(dim)]
			norm = sum(abs(entry) ** 2 for entry in vector) ** 0.5
			if norm > 1e-9:
				return [entry / norm for entry in vector]
		raise ValueError("the tableau does not describe a pure state")


class CliffordState:
	"""A stabilizer state, with the qubits it lives on."""

	def __init__(self, tableau: CliffordTableau, qubits: Sequence[Qid]) -> None:
		self._tableau = tableau
		self._qubits = tuple(qubits)

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits

	def state_vector(self) -> list[complex]:
		return self._tableau.state_vector()

	def stabilizers(self) -> list[Any]:
		return self._tableau.stabilizer_pauli_strings()

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
		tableau = CliffordTableau(len(qubits))
		for operation in program.all_operations():
			if isinstance(operation.gate, MeasurementGate):
				raise ValueError("simulate() cannot run a circuit with measurements; use run()")
			self._apply(tableau, operation, index)
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
		rng = random.Random(self._seed if seed is None else seed)
		records: dict[str, list[list[int]]] = {}
		for _ in range(repetitions):
			tableau = CliffordTableau(len(qubits))
			for operation in program.all_operations():
				if isinstance(operation.gate, MeasurementGate):
					self._measure(tableau, operation, index, rng, records)
				else:
					self._apply(tableau, operation, index)
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

	@staticmethod
	def _apply(tableau: CliffordTableau, operation: Operation, index: dict[Qid, int]) -> None:
		gate = operation.gate
		targets = [index[q] for q in operation.qubits]
		if isinstance(gate, IdentityGate):
			return
		if isinstance(gate, HPowGate) and gate.exponent == 1:
			tableau.hadamard(targets[0])
			return
		if isinstance(gate, ZPowGate):
			if gate.exponent == 1:
				tableau.pauli_z(targets[0])
			elif gate.exponent == 0.5:
				tableau.phase(targets[0])
			elif gate.exponent == -0.5:
				for _ in range(3):
					tableau.phase(targets[0])
			else:
				raise TypeError(f"{gate!r} is not a Clifford gate")
			return
		if isinstance(gate, XPowGate) and gate.exponent in (1, -1):
			tableau.pauli_x(targets[0])
			return
		if isinstance(gate, YPowGate) and gate.exponent in (1, -1):
			tableau.pauli_y(targets[0])
			return
		if isinstance(gate, CNOTGate):
			tableau.cnot(targets[0], targets[1])
			return
		if isinstance(gate, CZGate):
			tableau.cz(targets[0], targets[1])
			return
		if isinstance(gate, SWAPGate):
			tableau.swap(targets[0], targets[1])
			return
		raise TypeError(f"{gate!r} is not a Clifford gate")

	@staticmethod
	def _measure(
		tableau: CliffordTableau,
		operation: Operation,
		index: dict[Qid, int],
		rng: random.Random,
		records: dict[str, list[list[int]]],
	) -> None:
		gate = operation.gate
		assert isinstance(gate, MeasurementGate)
		invert_mask = gate.full_invert_mask
		bits = [tableau.measure(index[q], rng) ^ int(invert_mask[k]) for k, q in enumerate(operation.qubits)]
		records.setdefault(gate.key, []).append(bits)
