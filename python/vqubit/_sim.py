"""Local simulators.

The :class:`Simulator` evolves the state vector, the
:class:`DensityMatrixSimulator` a density matrix, and both do the heavy linear
algebra in the native kernels. States are kept as flat, interleaved buffers of
``float`` (``[re0, im0, re1, im1, ...]``) end to end, so a simulation does not
pay for converting Python objects on every gate and uses far less memory.

No third-party package is involved: a circuit is built, simulated and sampled
entirely on your machine.
"""

from __future__ import annotations

import random
from collections import Counter
from collections.abc import Callable, Sequence
from typing import Any

import vqubit_native as _native

from ._circuits import Circuit
from ._compat import HAS_NUMPY, is_numpy_array, native_complex_view
from ._gates import MeasurementGate
from ._linalg import to_flat, vector_from_flat, vector_to_flat
from ._ops import Moment
from ._qubits import Qid
from ._sim_abc import (
	SimulatesAmplitudes,
	SimulatesFinalState,
	SimulatesIntermediateState,
	SimulatesSamples,
	SimulationTrialResultBase,
	SimulatorBase,
)


def _is_basis_zero(initial_state: Any) -> bool:
	"""Report whether ``initial_state`` is the computational basis state ``|0...0>``."""
	return isinstance(initial_state, int) and not isinstance(initial_state, bool) and initial_state == 0


def _resolve_qubits(circuit: Circuit, qubit_order: Sequence[Qid] | None) -> tuple[Qid, ...]:
	if qubit_order is None:
		return circuit.all_qubits()
	qubits = tuple(qubit_order)
	missing = set(circuit.all_qubits()) - set(qubits)
	if missing:
		raise ValueError(f"qubit_order is missing qubits used by the circuit: {sorted(missing)!r}")
	return qubits


def _initial_state_vector(initial_state: Any, dim: int) -> list[complex]:
	if isinstance(initial_state, StateVector):
		return initial_state.state_vector()
	if isinstance(initial_state, int) and not isinstance(initial_state, bool):
		if not 0 <= initial_state < dim:
			raise ValueError(f"initial_state {initial_state} out of range for {dim} amplitudes")
		vector = [0j] * dim
		vector[initial_state] = 1.0 + 0j
		return vector
	sequence = list(initial_state)
	if len(sequence) != dim:
		raise ValueError(f"initial_state has length {len(sequence)}, expected {dim}")
	return [complex(x) for x in sequence]


class StateVector:
	"""A state vector together with the qubits it is expressed in."""

	def __init__(self, state: Sequence[complex], qubits: Sequence[Qid], *, owner: Any = None) -> None:
		qubits = tuple(qubits)
		if is_numpy_array(state):
			self._state = state
		else:
			self._state = list(state)
		if len(self._state) != 1 << len(qubits):
			raise ValueError(f"state of length {len(self._state)} does not match {len(qubits)} qubits")
		self._qubits = qubits
		self._owner = owner

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits

	def state_vector(self) -> Any:
		"""Return the amplitudes (a zero-copy view when numpy is installed)."""
		if is_numpy_array(self._state):
			return self._state
		return list(self._state)

	def dirac_notation(self, decimals: int = 3) -> str:
		"""Return a human-readable form such as ``0.707|00⟩ + 0.707|11⟩``."""
		terms: list[str] = []
		for i, amplitude in enumerate(self._state):
			if abs(amplitude) < 10 ** (-decimals):
				continue
			bits = format(i, f"0{len(self._qubits)}b")
			terms.append(f"{_format_complex(amplitude, decimals)}|{bits}⟩")
		return " + ".join(terms) if terms else "0"

	def __str__(self) -> str:
		return self.dirac_notation()

	def __repr__(self) -> str:
		return f"vqubit.StateVector({self.dirac_notation()})"


def _format_complex(value: complex, decimals: int) -> str:
	if abs(value.imag) < 10 ** (-decimals):
		return f"{value.real:.{decimals}f}"
	return f"{value:.{decimals}f}"


class SimulationResult(SimulationTrialResultBase):
	"""The result of :meth:`Simulator.simulate`."""

	def __init__(self, qubits: Sequence[Qid], state: StateVector) -> None:
		self._qubits = tuple(qubits)
		self._final_state = state

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits

	@property
	def final_state(self) -> StateVector:
		return self._final_state

	def state_vector(self) -> Any:
		return self._final_state.state_vector()

	def dirac_notation(self, decimals: int = 3) -> str:
		return self._final_state.dirac_notation(decimals)

	def __str__(self) -> str:
		return self._final_state.dirac_notation()

	def __repr__(self) -> str:
		return f"vqubit.SimulationResult({self._final_state!r})"


class Result:
	"""Measurement results: one row of bits per repetition, per key."""

	def __init__(self, measurements: dict[str, list[list[int]]]) -> None:
		self._measurements = {key: [list(row) for row in rows] for key, rows in measurements.items()}

	@property
	def measurements(self) -> dict[str, list[list[int]]]:
		return self._measurements

	def __len__(self) -> int:
		return len(next(iter(self._measurements.values()), []))

	def histogram(self, key: str = "", fold_func: Callable[[list[int]], Any] | None = None) -> Counter:
		"""Count how often each outcome appears for ``key``."""
		if not key:
			if len(self._measurements) != 1:
				raise ValueError("histogram needs a key when there is not exactly one measurement")
			key = next(iter(self._measurements))
		rows = self._measurements[key]
		if fold_func is None:

			def fold_func(bits: list[int]) -> int:
				return int("".join(str(b) for b in bits), 2) if bits else 0

		return Counter(fold_func(row) for row in rows)

	def __repr__(self) -> str:
		return f"vqubit.Result(keys={list(self._measurements)!r})"


def _record_bits(
	gate: MeasurementGate,
	outcome: int,
	num_targets: int,
	measurements: dict[str, list[list[int]]],
) -> None:
	bits = [(outcome >> (num_targets - 1 - k)) & 1 for k in range(num_targets)]
	invert_mask = gate.full_invert_mask
	bits = [bit ^ int(invert_mask[k]) for k, bit in enumerate(bits)]
	measurements.setdefault(gate.key, []).append(bits)


def _measurements_are_terminal(program: Circuit) -> bool:
	"""Report whether every measurement happens after all non-measurements.

	When that holds, the final state can be computed once and sampled directly
	instead of re-simulating the circuit for every repetition.
	"""
	operations = list(program.all_operations())
	for index, operation in enumerate(operations):
		if not isinstance(operation.gate, MeasurementGate):
			continue
		qubits = set(operation.qubits)
		for later in operations[index + 1 :]:
			if not isinstance(later.gate, MeasurementGate) and qubits.intersection(later.qubits):
				return False
	return True


def _draw_outcome(probabilities: Sequence[float], rng: random.Random) -> int:
	threshold = rng.random()
	accumulated = 0.0
	for index, probability in enumerate(probabilities):
		accumulated += probability
		if threshold < accumulated:
			return index
	return len(probabilities) - 1


class Simulator(SimulatorBase, SimulatesSamples, SimulatesAmplitudes, SimulatesFinalState, SimulatesIntermediateState):
	"""A state-vector simulator.

	``simulate`` evolves a circuit with no measurements and returns the final
	state. ``run`` samples measurements and returns a :class:`Result`.
	"""

	def __init__(self, *, seed: int | None = None) -> None:
		self._seed = seed

	def simulate(
		self,
		program: Circuit,
		*,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
	) -> SimulationResult:
		qubits = _resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		num_qubits = len(qubits)
		state = _native.NativeState(num_qubits)
		if _is_basis_zero(initial_state):
			state.set_zero_state()
		else:
			state.set_amplitudes(vector_to_flat(_initial_state_vector(initial_state, 1 << num_qubits)))
		for moment in program.moments:
			self._apply_moment(state, moment, index, rng=None, measurements=None)
		view = native_complex_view(state)
		if view is not None:
			return SimulationResult(qubits, StateVector(view, qubits, owner=state))
		return SimulationResult(qubits, StateVector(vector_from_flat(state.amplitudes()), qubits))

	def run(
		self,
		program: Circuit,
		repetitions: int = 1,
		*,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
		seed: int | None = None,
	) -> Result:
		if repetitions < 1:
			raise ValueError("repetitions must be at least 1")
		qubits = _resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		num_qubits = len(qubits)
		dim = 1 << num_qubits
		amplitudes = vector_to_flat(_initial_state_vector(initial_state, dim))
		rng = random.Random(self._seed if seed is None else seed)
		records: dict[str, list[list[int]]] = {}
		if _measurements_are_terminal(program):
			state = _native.NativeState(num_qubits)
			if _is_basis_zero(initial_state):
				state.set_zero_state()
			else:
				state.set_amplitudes(amplitudes)
			measurements: list[tuple[MeasurementGate, list[int]]] = []
			for moment in program.moments:
				for operation in moment.operations:
					if isinstance(operation.gate, MeasurementGate):
						measurements.append((operation.gate, [index[q] for q in operation.qubits]))
						continue
					matrix = operation.gate._cached_unitary_()
					if matrix is None:
						raise TypeError(f"{operation.gate!r} has no unitary and cannot be simulated")
					state.apply_gate(to_flat([list(row) for row in matrix]), [index[q] for q in operation.qubits])
			probabilities = [list(state.probabilities(targets)) for _, targets in measurements]
			if HAS_NUMPY:
				import numpy as _np

				np_rng = _np.random.default_rng(self._seed if seed is None else seed)
				for (gate, targets), probs in zip(measurements, probabilities, strict=True):
					cumulative = _np.cumsum(_np.asarray(probs, dtype=_np.float64))
					outcomes = _np.searchsorted(cumulative, np_rng.random(repetitions))
					for outcome in outcomes.tolist():
						_record_bits(gate, int(outcome), len(targets), records)
				return Result(records)
			for _ in range(repetitions):
				for (gate, targets), probs in zip(measurements, probabilities, strict=True):
					_record_bits(gate, _draw_outcome(probs, rng), len(targets), records)
			return Result(records)
		for _ in range(repetitions):
			state = _native.NativeState(num_qubits)
			state.set_amplitudes(amplitudes)
			for moment in program.moments:
				self._apply_moment(state, moment, index, rng=rng, measurements=records)
		return Result(records)

	sample = run

	def _apply_moment(
		self,
		state: Any,
		moment: Moment,
		index: dict[Qid, int],
		*,
		rng: random.Random | None,
		measurements: dict[str, list[list[int]]] | None,
	) -> None:
		for operation in moment.operations:
			gate = operation.gate
			targets = [index[q] for q in operation.qubits]
			if isinstance(gate, MeasurementGate):
				if rng is None or measurements is None:
					raise ValueError("simulate() cannot run a circuit with measurements; use run()")
				outcome = state.measure(targets, rng.random())
				_record_bits(gate, outcome, len(targets), measurements)
				continue
			matrix = gate._cached_unitary_()
			if matrix is None:
				raise TypeError(f"{gate!r} has no unitary and cannot be simulated")
			state.apply_gate(to_flat([list(row) for row in matrix]), targets)


class DensityMatrix:
	"""A density matrix together with the qubits it is expressed in."""

	def __init__(self, matrix: Any, qubits: Sequence[Qid], *, owner: Any = None) -> None:
		self._qubits = tuple(qubits)
		if is_numpy_array(matrix):
			self._matrix = matrix
		else:
			self._matrix = [[complex(entry) for entry in row] for row in matrix]
		self._owner = owner

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits

	def density_matrix(self) -> Any:
		"""Return the density matrix (a zero-copy view when numpy is installed)."""
		if is_numpy_array(self._matrix):
			return self._matrix
		return [list(row) for row in self._matrix]

	def dirac_notation(self, decimals: int = 3) -> str:
		return StateVector(
			[self._matrix[i][i].real ** 0.5 + 0j for i in range(len(self._matrix))], self._qubits
		).dirac_notation(decimals)

	def __repr__(self) -> str:
		return f"vqubit.DensityMatrix({len(self._qubits)} qubits)"


class DensityMatrixTrialResult(SimulationTrialResultBase):
	"""The result of :meth:`DensityMatrixSimulator.simulate`."""

	def __init__(self, qubits: Sequence[Qid], density: DensityMatrix) -> None:
		self._qubits = tuple(qubits)
		self._final_density_matrix = density

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits

	@property
	def final_density_matrix(self) -> DensityMatrix:
		return self._final_density_matrix

	def density_matrix(self) -> Any:
		return self._final_density_matrix.density_matrix()

	def __repr__(self) -> str:
		return f"vqubit.DensityMatrixTrialResult({self._final_density_matrix!r})"


def _density_flat_to_nested(flat: Sequence[float], dim: int) -> list[list[complex]]:
	return [[complex(flat[2 * (i * dim + j)], flat[2 * (i * dim + j) + 1]) for j in range(dim)] for i in range(dim)]


class DensityMatrixSimulator(SimulatorBase, SimulatesSamples, SimulatesFinalState, SimulatesIntermediateState):
	"""A simulator that evolves a density matrix, so it can run noisy circuits."""

	def __init__(self, *, seed: int | None = None) -> None:
		self._seed = seed

	def simulate(
		self,
		program: Circuit,
		*,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
	) -> DensityMatrixTrialResult:
		qubits = _resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		num_qubits = len(qubits)
		density = _native.NativeDensity(num_qubits)
		if _is_basis_zero(initial_state):
			density.set_zero_state()
		else:
			density.set_from_amplitudes(vector_to_flat(_initial_state_vector(initial_state, 1 << num_qubits)))
		for moment in program.moments:
			self._apply_moment(density, moment, index, rng=None, measurements=None)
		view = native_complex_view(density)
		if view is not None:
			matrix = DensityMatrix(view.reshape(1 << num_qubits, 1 << num_qubits), qubits, owner=density)
		else:
			matrix = DensityMatrix(_density_flat_to_nested(density.density_matrix(), 1 << num_qubits), qubits)
		return DensityMatrixTrialResult(qubits, matrix)

	def run(
		self,
		program: Circuit,
		repetitions: int = 1,
		*,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
		seed: int | None = None,
	) -> Result:
		if repetitions < 1:
			raise ValueError("repetitions must be at least 1")
		qubits = _resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		num_qubits = len(qubits)
		amplitudes = vector_to_flat(_initial_state_vector(initial_state, 1 << num_qubits))
		rng = random.Random(self._seed if seed is None else seed)
		records: dict[str, list[list[int]]] = {}
		for _ in range(repetitions):
			density = _native.NativeDensity(num_qubits)
			density.set_from_amplitudes(amplitudes)
			for moment in program.moments:
				self._apply_moment(density, moment, index, rng=rng, measurements=records)
		return Result(records)

	sample = run

	def _apply_moment(
		self,
		density: Any,
		moment: Moment,
		index: dict[Qid, int],
		*,
		rng: random.Random | None,
		measurements: dict[str, list[list[int]]] | None,
	) -> None:
		for operation in moment.operations:
			gate = operation.gate
			targets = [index[q] for q in operation.qubits]
			if isinstance(gate, MeasurementGate):
				if rng is None or measurements is None:
					raise ValueError("simulate() cannot run a circuit with measurements; use run()")
				outcome = density.measure(targets, rng.random())
				_record_bits(gate, outcome, len(targets), measurements)
				continue
			operators = gate._kraus_()
			if operators is None:
				matrix = gate._cached_unitary_()
				if matrix is None:
					raise TypeError(f"{gate!r} has neither a unitary nor Kraus operators")
				operators = [matrix]
			flat: list[float] = []
			for operator in operators:
				flat.extend(to_flat([list(row) for row in operator]))
			density.apply_channel(flat, len(operators), targets)
