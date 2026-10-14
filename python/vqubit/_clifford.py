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


class CliffordTableau:
	"""The binary tableau of a stabilizer state."""

	def __init__(self, num_qubits: int, flat: Sequence[int] | None = None) -> None:
		self.num_qubits = num_qubits
		self.rows = 2 * num_qubits + 1
		self._side = 2 * num_qubits + 1
		self._flat = list(flat) if flat is not None else []

	@classmethod
	def from_flat(cls, flat: Sequence[int], num_qubits: int) -> CliffordTableau:
		"""Rebuild a tableau from the flattened native form."""
		return cls(num_qubits, flat)

	def _bit(self, row: int, col: int) -> int:
		return self._flat[row * self._side + col]

	def stabilizer_pauli_strings(self) -> list[Any]:
		"""Return the stabilizer generators as Pauli strings."""
		from ._pauli import PauliString

		result: list[Any] = []
		for k in range(self.num_qubits):
			row = self.num_qubits + k
			qubit_pauli_map: dict[Qid, Pauli] = {}
			for j in range(self.num_qubits):
				xbit = self._bit(row, j)
				zbit = self._bit(row, self.num_qubits + j)
				pauli = _PAULI_BITS[(xbit, zbit)]
				if pauli is not Pauli.I:
					qubit_pauli_map[LineQubit(j)] = pauli
			sign = -1 if self._bit(row, 2 * self.num_qubits) else 1
			result.append(PauliString(qubit_pauli_map, coefficient=sign))
		return result


_PAULI_BITS = {
	(0, 0): Pauli.I,
	(1, 0): Pauli.X,
	(1, 1): Pauli.Y,
	(0, 1): Pauli.Z,
}


class CliffordState:
	"""A stabilizer state, with the qubits it lives on."""

	def __init__(self, tableau: CliffordTableau, qubits: Sequence[Qid]) -> None:
		self._tableau = tableau
		self._qubits = tuple(qubits)

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits

	def state_vector(self) -> list[complex]:
		vector = _native.clifford_state_vector(self._tableau._flat, self._tableau.num_qubits)
		return vector_from_flat(vector)

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
		instructions, _, _ = _build_instructions(program, index, allow_measurements=False)
		flat = _native.clifford_tableau(instructions, len(qubits))
		tableau = CliffordTableau.from_flat(flat, len(qubits))
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
