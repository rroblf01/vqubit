"""Step-by-step simulation: states, step results, trial results and helpers.

This mirrors the ``cirq.sim`` step API on top of the native buffers used by
:mod:`vqubit._sim`.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence
from typing import Any

import vqubit_native as _native

from ._compat import asarray, asvector
from ._gates import MeasurementGate
from ._linalg import to_flat, vector_from_flat, vector_to_flat
from ._qis_extra import density_matrix_from_state_vector
from ._qubits import Qid
from ._sim import StateVector, _initial_state_vector, _qubit_pauli, _reduced_density_matrix
from ._sim_abc import (
	SimulationStateBase,
	SimulationTrialResultBase,
	StepResultBase,
)
from ._value import value_equality


def _probabilities_from_state(state: Sequence[complex], indices: Sequence[int]) -> list[float]:
	flat = vector_to_flat([complex(z) for z in state])
	num_qubits = len(state).bit_length() - 1
	return list(_native.state_probabilities(flat, list(indices), num_qubits))


def _collapse_from_state(state: Sequence[complex], indices: Sequence[int], outcome: int) -> list[complex]:
	flat = vector_to_flat([complex(z) for z in state])
	num_qubits = len(state).bit_length() - 1
	# The kernel returns an all-zero vector for an impossible outcome rather
	# than an empty one, so the probability has to be checked up front.
	probabilities = list(_native.state_probabilities(flat, list(indices), num_qubits))
	if not 0 <= outcome < len(probabilities) or probabilities[outcome] <= 0:
		raise ValueError(f"outcome {outcome} has probability zero")
	return vector_from_flat(_native.collapse_state(flat, list(indices), outcome, num_qubits))


def _checked_num_qubits(size: int, qid_shape: Sequence[int] | None, what: str) -> int:
	"""Validate the shape of a state and return its number of qubits."""
	if qid_shape is None:
		if size == 0 or size & (size - 1):
			raise ValueError(f"{what} has {size} amplitudes, which is not a power of two")
		return size.bit_length() - 1
	dims = tuple(int(dimension) for dimension in qid_shape)
	total = math.prod(dims)
	if total != size:
		raise ValueError(f"{what} has {size} amplitudes, but qid_shape {dims} needs {total}")
	if any(dimension != 2 for dimension in dims):
		raise ValueError("the native kernels measure qubits only; qid_shape must be all twos")
	return len(dims)


def _draw(probabilities: Sequence[float], rng: random.Random) -> int:
	threshold = rng.random()
	accumulated = 0.0
	for index, probability in enumerate(probabilities):
		accumulated += probability
		if threshold < accumulated:
			return index
	return len(probabilities) - 1


def measure_state_vector(
	state: Sequence[complex],
	indices: Sequence[int],
	*,
	qid_shape: Sequence[int] | None = None,
	out: Any = None,
	seed: int | None = None,
) -> tuple[list[int], list[complex]]:
	"""Measure ``indices`` of a state vector; return the outcomes and the collapsed state."""
	_checked_num_qubits(len(state), qid_shape, "state_vector")
	probabilities = _probabilities_from_state(state, indices)
	outcome = _draw(probabilities, random.Random(seed))
	collapsed = _collapse_from_state(state, indices, outcome)
	if out is not None:
		out[:] = collapsed
	return [outcome], collapsed


def measure_density_matrix(
	density_matrix: Sequence[Sequence[complex]],
	indices: Sequence[int],
	*,
	qid_shape: Sequence[int] | None = None,
	out: Any = None,
	seed: int | None = None,
) -> tuple[list[int], list[list[complex]]]:
	"""Measure ``indices`` of a density matrix; return the outcomes and the collapsed matrix."""
	_checked_num_qubits(len(density_matrix), qid_shape, "density_matrix")
	dim = len(density_matrix)
	num_qubits = dim.bit_length() - 1
	flat = to_flat([[complex(z) for z in row] for row in density_matrix])
	probabilities = list(_native.density_probabilities(flat, list(indices), num_qubits))
	outcome = _draw(probabilities, random.Random(seed))
	if not 0 <= outcome < len(probabilities) or probabilities[outcome] <= 0:
		raise ValueError(f"outcome {outcome} has probability zero")
	collapsed_flat = _native.density_collapse(flat, list(indices), outcome, num_qubits)
	collapsed = [
		[complex(collapsed_flat[2 * (i * dim + j)], collapsed_flat[2 * (i * dim + j) + 1]) for j in range(dim)]
		for i in range(dim)
	]
	if out is not None:
		out[:] = collapsed
	return [outcome], collapsed


def sample_state_vector(
	state: Sequence[complex],
	indices: Sequence[int],
	repetitions: int = 1,
	*,
	seed: int | None = None,
	qid_shape: Sequence[int] | None = None,
) -> list[int]:
	"""Sample measurement outcomes from a state vector."""
	_checked_num_qubits(len(state), qid_shape, "state_vector")
	probabilities = _probabilities_from_state(state, indices)
	rng = random.Random(seed)
	return [_draw(probabilities, rng) for _ in range(repetitions)]


def sample_density_matrix(
	density_matrix: Sequence[Sequence[complex]],
	indices: Sequence[int],
	repetitions: int = 1,
	*,
	seed: int | None = None,
	qid_shape: Sequence[int] | None = None,
) -> list[int]:
	"""Sample measurement outcomes from a density matrix."""
	_checked_num_qubits(len(density_matrix), qid_shape, "density_matrix")
	dim = len(density_matrix)
	num_qubits = dim.bit_length() - 1
	flat = to_flat([[complex(z) for z in row] for row in density_matrix])
	probabilities = list(_native.density_probabilities(flat, list(indices), num_qubits))
	rng = random.Random(seed)
	return [_draw(probabilities, rng) for _ in range(repetitions)]


def final_state_vector(
	program: Any,
	*,
	initial_state: Any = 0,
	param_resolver: Any = None,
	qubit_order: Any = None,
	ignore_terminal_measurements: bool = False,
	dtype: Any = None,
	seed: int | None = None,
) -> list[complex]:
	"""Return the final state vector of a circuit, as Cirq's helper does.

	A simulation *result* is accepted too, which is what this used to take: it
	returns its state without simulating anything.
	"""
	del dtype
	from ._circuits import Circuit

	if hasattr(program, "state_vector") and not isinstance(program, Circuit):
		return list(program.state_vector())
	from ._sim import Simulator

	circuit = program
	if ignore_terminal_measurements:
		circuit = _without_terminal_measurements(circuit)
	elif any(isinstance(operation.gate, MeasurementGate) for operation in circuit.all_operations()):
		# Cirq refuses here too: a measurement leaves a mixed state behind, so
		# there is no single state vector to return unless they are ignored.
		raise ValueError(
			"the circuit has measurements, and a state vector cannot hold them; "
			"pass ignore_terminal_measurements=True for measurements at the end, or simulate with run()"
		)
	return list(
		Simulator(seed=seed)
		.simulate(circuit, param_resolver=param_resolver, qubit_order=qubit_order, initial_state=initial_state)
		.state_vector()
	)


def final_density_matrix(
	program: Any,
	*,
	noise: Any = None,
	initial_state: Any = 0,
	param_resolver: Any = None,
	qubit_order: Any = None,
	dtype: Any = None,
	seed: int | None = None,
	ignore_measurement_results: bool = True,
) -> list[list[complex]]:
	"""Return the final density matrix of a circuit, as Cirq's helper does.

	A simulation *result* is accepted too, which is what this used to take: it
	returns its matrix without simulating anything. ``noise`` is accepted and
	ignored: add the channels to the circuit instead, which is what a noise model
	expands to.
	"""
	del dtype, ignore_measurement_results
	from ._circuits import Circuit

	if hasattr(program, "density_matrix") and not isinstance(program, Circuit):
		return [list(row) for row in program.density_matrix()]
	circuit = program
	if noise is not None:
		# A noise model is exactly a map from moments to moments with channels.
		circuit = Circuit(noise.noisy_moments(list(circuit.moments), sorted(circuit.all_qubits())))
	from ._sim import DensityMatrixSimulator

	return [
		list(row)
		for row in DensityMatrixSimulator(seed=seed)
		.simulate(circuit, param_resolver=param_resolver, qubit_order=qubit_order, initial_state=initial_state)
		.density_matrix()
	]


def _without_terminal_measurements(circuit: Any) -> Any:
	"""Drop the measurements that nothing follows, which a state vector cannot carry."""
	from ._circuits import Circuit
	from ._sim import _measurements_are_terminal

	if not _measurements_are_terminal(circuit):
		raise ValueError(
			"the circuit has measurements that are not at the end; a state vector cannot hold them "
			"(pass ignore_terminal_measurements=True only for terminal ones, or simulate with run())"
		)
	operations = [
		operation for operation in circuit.all_operations() if not isinstance(operation.gate, MeasurementGate)
	]
	return Circuit(operations)


def sample_sweep(program: Any, params: Any, repetitions: int = 1, **kwargs: Any) -> Any:
	"""Sample a parameterized circuit over a sweep."""
	from ._protocols import resolve_parameters
	from ._study import sample, to_resolvers, to_sweep

	return [
		sample(resolve_parameters(program, resolver), repetitions=repetitions, **kwargs)
		for resolver in to_resolvers(to_sweep(params))
	]


class StateVectorMixin:
	"""Mixin providing state-vector conveniences."""

	def dirac_notation(self, decimals: int = 3) -> str:
		"""The state written in Dirac notation."""
		subject: Any = self
		return StateVector(subject.state_vector(), subject.qubits).dirac_notation(decimals)

	def __str__(self) -> str:
		return self.dirac_notation()


class StepResult(StepResultBase):
	"""The result of simulating a single step."""

	def __init__(
		self,
		*,
		qubits: Sequence[Qid],
		measurements: dict[str, Any] | None = None,
		state_vector: Sequence[complex] | None = None,
		density_matrix: Sequence[Sequence[complex]] | None = None,
	) -> None:
		self._qubits = tuple(qubits)
		self.measurements = measurements if measurements is not None else {}
		self._state_vector = None if state_vector is None else [complex(z) for z in state_vector]
		self._density_matrix = None if density_matrix is None else [[complex(z) for z in row] for row in density_matrix]

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits of the state, in order."""
		return self._qubits

	@property
	def qubit_map(self) -> dict[Qid, int]:
		"""Which index each qubit holds in the state, as Cirq reports it."""
		return {qubit: index for index, qubit in enumerate(self._qubits)}

	def state_vector(self) -> list[complex]:
		"""The state as a flat vector of amplitudes."""
		if self._state_vector is None:
			raise ValueError("this step result has no state vector")
		return asvector(self._state_vector)

	def density_matrix(self) -> list[list[complex]]:
		"""The state as a density matrix."""
		if self._density_matrix is not None:
			return asarray(self._density_matrix)
		return density_matrix_from_state_vector(self.state_vector())

	def set_state_vector(self, state: Sequence[complex]) -> None:
		"""Replace the state with ``state``."""
		self._state_vector = [complex(z) for z in state]

	def dirac_notation(self, decimals: int = 3) -> str:
		"""The state written in Dirac notation."""
		return StateVector(self.state_vector(), self._qubits).dirac_notation(decimals)

	def sample(self, qubits: Sequence[Qid], repetitions: int = 1, *, seed: int | None = None) -> list[int]:
		"""Sample a measurement of the given qubits."""
		indices = [self._qubits.index(q) for q in qubits]
		return sample_state_vector(self.state_vector(), indices, repetitions, seed=seed)

	def sample_measurement_ops(
		self, measurement_ops: Sequence[Any], repetitions: int = 1, seed: Any = None
	) -> dict[str, list[list[int]]]:
		"""Sample the given measurement operations on this step's state.

		Each repetition starts from this step's state and applies the measurements
		in order (so they collapse it, as a real run would), which is what Cirq's
		step results do.
		"""
		prng = random.Random(seed) if not isinstance(seed, random.Random) else seed
		records: dict[str, list[list[int]]] = {}
		for _ in range(repetitions):
			state = _native.NativeState(len(self._qubits))
			state.set_amplitudes(vector_to_flat(self.state_vector()))
			for operation in measurement_ops:
				gate = operation.gate
				targets = [self._qubits.index(qubit) for qubit in operation.qubits]
				outcome = state.measure(targets, prng.random())
				bits = [(outcome >> (len(targets) - 1 - k)) & 1 for k in range(len(targets))]
				mask = gate.full_invert_mask
				rows = records.setdefault(gate.key, [])
				rows.append([bit ^ int(mask[k]) for k, bit in enumerate(bits)])
		return records

	def __repr__(self) -> str:
		return f"vqubit.StepResult({len(self._qubits)} qubits)"


class StateVectorStepResult(StepResult, StateVectorMixin):
	"""A step result holding a state vector."""

	def __init__(self, *, qubits: Sequence[Qid], measurements: dict[str, Any], state_vector: Sequence[complex]) -> None:
		super().__init__(qubits=qubits, measurements=measurements, state_vector=state_vector)

	def bloch_vector_of(self, qubit: Qid) -> list[float]:
		"""The Bloch vector of ``qubit`` in this step's state: ``[<X>, <Y>, <Z>]``."""
		return [
			_qubit_pauli(qubit, letter)
			.expectation_from_state_vector(self.state_vector(), self.qubit_map, check_preconditions=False)
			.real
			for letter in "XYZ"
		]

	def density_matrix_of(self, qubits: Sequence[Qid] | None = None) -> Any:
		"""The reduced density matrix of ``qubits`` in this step (all by default)."""
		return _reduced_density_matrix(self.state_vector(), self.qubit_map, qubits)


class DensityMatrixStepResult(StepResult):
	"""A step result holding a density matrix."""

	def __init__(
		self, *, qubits: Sequence[Qid], measurements: dict[str, Any], density_matrix: Sequence[Sequence[complex]]
	) -> None:
		super().__init__(qubits=qubits, measurements=measurements, density_matrix=density_matrix)


class CliffordSimulatorStepResult(StepResult):
	"""A step result from a Clifford simulation."""

	def __init__(self, *, qubits: Sequence[Qid], measurements: dict[str, Any], state: Any) -> None:
		self._clifford_state = state
		super().__init__(qubits=qubits, measurements=measurements, state_vector=state.state_vector())

	@property
	def state(self) -> Any:
		"""The Clifford state of this step, under the name Cirq gives it."""
		return self._clifford_state

	@property
	def clifford_state(self) -> Any:
		"""The stabilizer state after this step."""
		return self._clifford_state


class SparseSimulatorStep(StateVectorStepResult):
	"""A step result from a sparse state-vector simulation."""


class SimulationState(SimulationStateBase):
	"""The mutable state of a step-by-step simulation."""

	def __init__(self, qubits: Sequence[Qid], prng: random.Random | None = None) -> None:
		self._qubits = tuple(qubits)
		self._index = {q: i for i, q in enumerate(self._qubits)}
		self.prng = prng if prng is not None else random.Random()
		self.log_of_measurement_results: dict[str, Any] = {}

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits the state lives on."""
		return self._qubits

	def apply_operation(self, operation: Any) -> None:
		"""Apply one operation to the state."""
		raise NotImplementedError

	def measure(self, qubits: Sequence[Qid], key: str) -> int:
		"""Measure ``qubits`` and record the outcome under ``key``."""
		raise NotImplementedError

	def reset(self, qubits: Sequence[Qid]) -> None:
		"""Reset ``qubits`` to ``|0>``."""
		raise NotImplementedError

	def create_merged_state(self) -> SimulationState:
		"""Return a state that can hold a tensor product with this one."""
		return self

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}({len(self._qubits)} qubits)"


class StateVectorSimulationState(SimulationState):
	"""A mutable state-vector simulation state."""

	def __init__(self, *, qubits: Sequence[Qid], prng: random.Random | None = None, initial_state: Any = 0) -> None:
		super().__init__(qubits, prng)
		self._state = _native.NativeState(len(self._qubits))
		self._state.set_amplitudes(vector_to_flat(_initial_state_vector(initial_state, 1 << len(self._qubits))))

	def state_vector(self) -> list[complex]:
		"""The state as a flat vector of amplitudes."""
		return vector_from_flat(self._state.amplitudes())

	def apply_operation(self, operation: Any) -> None:
		"""Apply one operation to the state."""
		gate = operation.gate
		targets = [self._index[q] for q in operation.qubits]
		if isinstance(gate, MeasurementGate):
			outcome = self._state.measure(targets, self.prng.random())
			self.log_of_measurement_results[gate.key] = outcome
			return
		matrix = gate._unitary_()
		if matrix is None:
			raise TypeError(f"{gate!r} has no unitary")
		self._state.apply_gate(to_flat([list(row) for row in matrix]), targets)

	def measure(self, qubits: Sequence[Qid], key: str) -> int:
		"""Measure ``qubits`` and record the outcome under ``key``."""
		targets = [self._index[q] for q in qubits]
		outcome = self._state.measure(targets, self.prng.random())
		self.log_of_measurement_results[key] = outcome
		return outcome


class DensityMatrixSimulationState(SimulationState):
	"""A mutable density-matrix simulation state."""

	def __init__(self, *, qubits: Sequence[Qid], prng: random.Random | None = None, initial_state: Any = 0) -> None:
		super().__init__(qubits, prng)
		self._state = _native.NativeDensity(len(self._qubits))
		self._state.set_from_amplitudes(vector_to_flat(_initial_state_vector(initial_state, 1 << len(self._qubits))))

	def density_matrix(self) -> list[list[complex]]:
		"""The state as a density matrix."""
		dim = 1 << len(self._qubits)
		flat = self._state.density_matrix()
		return [[complex(flat[2 * (i * dim + j)], flat[2 * (i * dim + j) + 1]) for j in range(dim)] for i in range(dim)]

	def apply_operation(self, operation: Any) -> None:
		"""Apply one operation to the state."""
		gate = operation.gate
		targets = [self._index[q] for q in operation.qubits]
		if isinstance(gate, MeasurementGate):
			outcome = self._state.measure(targets, self.prng.random())
			self.log_of_measurement_results[gate.key] = outcome
			return
		operators = gate._kraus_()
		if operators is None:
			matrix = gate._unitary_()
			if matrix is None:
				raise TypeError(f"{gate!r} has neither a unitary nor Kraus operators")
			operators = [matrix]
		flat: list[float] = []
		for operator in operators:
			flat.extend(to_flat([list(row) for row in operator]))
		self._state.apply_channel(flat, len(operators), targets)

	def measure(self, qubits: Sequence[Qid], key: str) -> int:
		"""Measure ``qubits`` and record the outcome under ``key``."""
		targets = [self._index[q] for q in qubits]
		outcome = self._state.measure(targets, self.prng.random())
		self.log_of_measurement_results[key] = outcome
		return outcome


class CliffordTableauSimulationState(SimulationState):
	"""A mutable Clifford-tableau simulation state."""

	def __init__(self, *, tableau: Any, qubits: Sequence[Qid], prng: random.Random | None = None) -> None:
		super().__init__(qubits, prng)
		self._tableau = tableau

	@property
	def tableau(self) -> Any:
		"""The stabilizer tableau behind the state."""
		return self._tableau

	def state_vector(self) -> list[complex]:
		"""The state as a flat vector of amplitudes."""
		return vector_from_flat(self._tableau.state_vector())

	def apply_operation(self, operation: Any) -> None:
		"""Apply one Clifford operation to the tablea."""
		from ._clifford import _instruction_word

		gate = operation.gate
		word = _instruction_word(gate)
		if word is None:
			raise TypeError(f"{gate!r} is not a Clifford gate")
		targets = [self._index[q] for q in operation.qubits]
		for code, positions in word:
			first = targets[positions[0]]
			second = targets[positions[1]] if len(positions) > 1 else first
			self._tableau.apply(code, first, second)


class StabilizerSimulationState(CliffordTableauSimulationState):
	"""A stabilizer simulation state (alias of the tableau state)."""


@value_equality
class StabilizerStateChForm:
	"""A stabilizer state in CH form, stored by its state vector.

	vqubit keeps the amplitudes directly (its native kernels are state-vector
	first) where Cirq stores the C and H matrices; the state it represents is
	the same, which is what the simulator and the sampling methods use.
	"""

	def __init__(self, state_vector: Sequence[complex], num_qubits: int | None = None) -> None:
		self._state = [complex(z) for z in state_vector]
		self._num_qubits = num_qubits if num_qubits is not None else len(self._state).bit_length() - 1

	@property
	def num_qubits(self) -> int:
		"""The number of qubits the state spans."""
		return self._num_qubits

	def state_vector(self) -> list[complex]:
		"""The state as a flat vector of amplitudes."""
		return list(self._state)

	def density_matrix(self) -> list[list[complex]]:
		"""The state as a density matrix."""
		return density_matrix_from_state_vector(self._state)

	def _value_equality_values_(self) -> Any:
		return (tuple(self._state),)

	def __repr__(self) -> str:
		return f"vqubit.StabilizerStateChForm({self._num_qubits} qubits)"


class StabilizerChFormSimulationState(SimulationState):
	"""A simulation state backed by a CH-form stabilizer state.

	vqubit's CH form stores the amplitudes (see :class:`StabilizerStateChForm`),
	so an operation is applied the same way the state-vector state applies it.
	"""

	def __init__(self, *, qubits: Sequence[Qid], prng: random.Random | None = None, initial_state: Any = 0) -> None:
		super().__init__(qubits, prng)
		self._vector = _initial_state_vector(initial_state, 1 << len(self._qubits))
		self._state = StabilizerStateChForm(self._vector, len(self._qubits))

	def state_vector(self) -> list[complex]:
		"""The state as a flat vector of amplitudes."""
		return self._state.state_vector()

	def apply_operation(self, operation: Any) -> None:
		"""Apply an operation to the state, measuring if it is a measurement."""
		gate = operation.gate
		targets = [self._index[q] for q in operation.qubits]
		if isinstance(gate, MeasurementGate):
			outcome = self._measure_targets(targets)
			self.log_of_measurement_results[gate.key] = outcome
			return
		matrix = gate._unitary_()
		if matrix is None:
			raise TypeError(f"{gate!r} has no unitary")
		total = 1 << len(self._qubits)
		state = _native.NativeState(len(self._qubits))
		state.set_amplitudes(vector_to_flat(self._state.state_vector()))
		state.apply_gate(to_flat([list(row) for row in matrix]), targets)
		amplitudes = vector_from_flat(state.amplitudes())
		if len(amplitudes) == total:
			self._state = StabilizerStateChForm(amplitudes, len(self._qubits))

	def measure(self, qubits: Sequence[Qid], key: str) -> int:
		"""Measure ``qubits`` and record the outcome under ``key``."""
		outcome = self._measure_targets([self._index[q] for q in qubits])
		self.log_of_measurement_results[key] = outcome
		return outcome

	def _measure_targets(self, targets: Sequence[int]) -> int:
		"""Sample the ``targets`` and collapse the stored state onto the outcome."""
		state = _native.NativeState(len(self._qubits))
		state.set_amplitudes(vector_to_flat(self._state.state_vector()))
		outcome = state.measure(list(targets), self.prng.random())
		self._state = StabilizerStateChForm(vector_from_flat(state.amplitudes()), len(self._qubits))
		return outcome


class SimulationProductState:
	"""A tensor product of independent simulation states."""

	def __init__(self, states: Sequence[SimulationStateBase]) -> None:
		self._states = tuple(states)

	@property
	def states(self) -> tuple[SimulationStateBase, ...]:
		"""The one state per qubit that make up the product."""
		return self._states

	def __repr__(self) -> str:
		return f"vqubit.SimulationProductState({len(self._states)} states)"


class SimulationTrialResult(SimulationTrialResultBase, StepResult):
	"""The result of a full simulation trial."""

	def __init__(
		self,
		params: Any,
		measurements: dict[str, Any],
		*,
		qubits: Sequence[Qid],
		state_vector: Sequence[complex] | None = None,
		density_matrix: Sequence[Sequence[complex]] | None = None,
	) -> None:
		self._params = params
		super().__init__(
			qubits=qubits, measurements=measurements, state_vector=state_vector, density_matrix=density_matrix
		)

	@property
	def params(self) -> Any:
		"""The parameter resolver the simulation used."""
		return self._params


class StateVectorTrialResult(SimulationTrialResult, StateVectorMixin):
	"""A trial result holding a final state vector."""

	def __init__(
		self, params: Any, measurements: dict[str, Any], *, qubits: Sequence[Qid], state_vector: Sequence[complex]
	) -> None:
		super().__init__(params, measurements, qubits=qubits, state_vector=state_vector)


class CliffordTrialResult(SimulationTrialResult):
	"""A trial result from a Clifford simulation."""

	def __init__(self, params: Any, measurements: dict[str, Any], *, qubits: Sequence[Qid], final_state: Any) -> None:
		self._clifford_state = final_state
		super().__init__(params, measurements, qubits=qubits, state_vector=final_state.state_vector())

	@property
	def final_state(self) -> Any:
		"""The stabilizer state at the end of the simulation."""
		return self._clifford_state


class StabilizerSampler:
	"""A sampler backed by the Clifford simulator."""

	def __init__(self, *, seed: int | None = None) -> None:
		self._seed = seed

	def run(self, program: Any, repetitions: int = 1, **kwargs: Any) -> Any:
		"""Run a stabilizer circuit and sample its measurements."""
		from ._clifford import CliffordSimulator

		return CliffordSimulator(seed=self._seed).run(program, repetitions=repetitions, **kwargs)

	def sample(self, program: Any, repetitions: int = 1, **kwargs: Any) -> Any:
		"""Run a stabilizer circuit and sample its measurements."""
		return self.run(program, repetitions=repetitions, **kwargs)

	sample_expectation_values = None


__all__ = [
	"CliffordSimulatorStepResult",
	"CliffordTableauSimulationState",
	"CliffordTrialResult",
	"DensityMatrixSimulationState",
	"DensityMatrixStepResult",
	"SimulationProductState",
	"SimulationState",
	"SimulationTrialResult",
	"SparseSimulatorStep",
	"StabilizerChFormSimulationState",
	"StabilizerSampler",
	"StabilizerSimulationState",
	"StabilizerStateChForm",
	"StateVectorMixin",
	"StateVectorSimulationState",
	"StateVectorStepResult",
	"StateVectorTrialResult",
	"StepResult",
	"final_density_matrix",
	"final_state_vector",
	"measure_density_matrix",
	"measure_state_vector",
	"sample_density_matrix",
	"sample_state_vector",
	"sample_sweep",
]
