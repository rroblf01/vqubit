"""Step-by-step simulation: states, step results, trial results and helpers.

This mirrors the ``cirq.sim`` step API on top of the native buffers used by
:mod:`vqubit._sim`.
"""

from __future__ import annotations

import random
from collections.abc import Sequence
from typing import Any

import vqubit_native as _native

from ._compat import asarray, asvector
from ._gates import MeasurementGate
from ._linalg import to_flat, vector_from_flat, vector_to_flat
from ._qis_extra import density_matrix_from_state_vector
from ._qubits import Qid
from ._sim import StateVector, _initial_state_vector
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
	collapsed = _native.collapse_state(flat, list(indices), outcome, num_qubits)
	if not collapsed:
		raise ValueError(f"outcome {outcome} has probability zero")
	return vector_from_flat(collapsed)


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
) -> tuple[list[complex], int]:
	"""Measure ``indices`` of a state vector; return the collapsed state and outcome."""
	del qid_shape
	probabilities = _probabilities_from_state(state, indices)
	outcome = _draw(probabilities, random.Random(seed))
	collapsed = _collapse_from_state(state, indices, outcome)
	if out is not None:
		out[:] = collapsed
	return collapsed, outcome


def measure_density_matrix(
	density_matrix: Sequence[Sequence[complex]],
	indices: Sequence[int],
	*,
	qid_shape: Sequence[int] | None = None,
	out: Any = None,
	seed: int | None = None,
) -> tuple[list[list[complex]], int]:
	"""Measure ``indices`` of a density matrix; return the collapsed matrix and outcome."""
	del qid_shape
	dim = len(density_matrix)
	num_qubits = dim.bit_length() - 1
	flat = to_flat([[complex(z) for z in row] for row in density_matrix])
	probabilities = list(_native.density_probabilities(flat, list(indices), num_qubits))
	outcome = _draw(probabilities, random.Random(seed))
	collapsed_flat = _native.density_collapse(flat, list(indices), outcome, num_qubits)
	if not collapsed_flat:
		raise ValueError(f"outcome {outcome} has probability zero")
	collapsed = [
		[complex(collapsed_flat[2 * (i * dim + j)], collapsed_flat[2 * (i * dim + j) + 1]) for j in range(dim)]
		for i in range(dim)
	]
	if out is not None:
		out[:] = collapsed
	return collapsed, outcome


def sample_state_vector(
	state: Sequence[complex],
	indices: Sequence[int],
	repetitions: int = 1,
	*,
	seed: int | None = None,
	qid_shape: Sequence[int] | None = None,
) -> list[int]:
	"""Sample measurement outcomes from a state vector."""
	del qid_shape
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
	del qid_shape
	dim = len(density_matrix)
	num_qubits = dim.bit_length() - 1
	flat = to_flat([[complex(z) for z in row] for row in density_matrix])
	probabilities = list(_native.density_probabilities(flat, list(indices), num_qubits))
	rng = random.Random(seed)
	return [_draw(probabilities, rng) for _ in range(repetitions)]


def final_state_vector(result: Any) -> list[complex]:
	"""Return the final state vector of a simulation result."""
	return list(result.state_vector())


def final_density_matrix(result: Any) -> list[list[complex]]:
	"""Return the final density matrix of a simulation result."""
	return [list(row) for row in result.density_matrix()]


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
		return self._qubits

	def state_vector(self) -> list[complex]:
		if self._state_vector is None:
			raise ValueError("this step result has no state vector")
		return asvector(self._state_vector)

	def density_matrix(self) -> list[list[complex]]:
		if self._density_matrix is not None:
			return asarray(self._density_matrix)
		return density_matrix_from_state_vector(self.state_vector())

	def set_state_vector(self, state: Sequence[complex]) -> None:
		self._state_vector = [complex(z) for z in state]

	def dirac_notation(self, decimals: int = 3) -> str:
		return StateVector(self.state_vector(), self._qubits).dirac_notation(decimals)

	def sample(self, qubits: Sequence[Qid], repetitions: int = 1, *, seed: int | None = None) -> list[int]:
		indices = [self._qubits.index(q) for q in qubits]
		return sample_state_vector(self.state_vector(), indices, repetitions, seed=seed)

	def __repr__(self) -> str:
		return f"vqubit.StepResult({len(self._qubits)} qubits)"


class StateVectorStepResult(StepResult, StateVectorMixin):
	"""A step result holding a state vector."""

	def __init__(self, *, qubits: Sequence[Qid], measurements: dict[str, Any], state_vector: Sequence[complex]) -> None:
		super().__init__(qubits=qubits, measurements=measurements, state_vector=state_vector)


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
	def clifford_state(self) -> Any:
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
		return self._qubits

	def apply_operation(self, operation: Any) -> None:
		raise NotImplementedError

	def measure(self, qubits: Sequence[Qid], key: str) -> int:
		raise NotImplementedError

	def reset(self, qubits: Sequence[Qid]) -> None:
		raise NotImplementedError

	def create_merged_state(self) -> SimulationState:
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
		return vector_from_flat(self._state.amplitudes())

	def apply_operation(self, operation: Any) -> None:
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
		dim = 1 << len(self._qubits)
		flat = self._state.density_matrix()
		return [[complex(flat[2 * (i * dim + j)], flat[2 * (i * dim + j) + 1]) for j in range(dim)] for i in range(dim)]

	def apply_operation(self, operation: Any) -> None:
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
		return self._tableau

	def state_vector(self) -> list[complex]:
		return vector_from_flat(self._tableau.state_vector())

	def apply_operation(self, operation: Any) -> None:
		from ._clifford import _instruction_code

		gate = operation.gate
		code = _instruction_code(gate)
		if code is None:
			raise TypeError(f"{gate!r} is not a Clifford gate")
		targets = [self._index[q] for q in operation.qubits]
		if len(targets) == 1:
			self._tableau.apply(code, targets[0], targets[0])
		else:
			self._tableau.apply(code, targets[0], targets[1])


class StabilizerSimulationState(CliffordTableauSimulationState):
	"""A stabilizer simulation state (alias of the tableau state)."""


@value_equality
class StabilizerStateChForm:
	"""A stabilizer state in CH form, stored by its state vector."""

	def __init__(self, state_vector: Sequence[complex], num_qubits: int | None = None) -> None:
		self._state = [complex(z) for z in state_vector]
		self._num_qubits = num_qubits if num_qubits is not None else len(self._state).bit_length() - 1

	@property
	def num_qubits(self) -> int:
		return self._num_qubits

	def state_vector(self) -> list[complex]:
		return list(self._state)

	def density_matrix(self) -> list[list[complex]]:
		return density_matrix_from_state_vector(self._state)

	def _value_equality_values_(self) -> Any:
		return (tuple(self._state),)

	def __repr__(self) -> str:
		return f"vqubit.StabilizerStateChForm({self._num_qubits} qubits)"


class StabilizerChFormSimulationState(SimulationState):
	"""A simulation state backed by a CH-form stabilizer state."""

	def __init__(self, *, qubits: Sequence[Qid], prng: random.Random | None = None, initial_state: Any = 0) -> None:
		super().__init__(qubits, prng)
		vector = _initial_state_vector(initial_state, 1 << len(self._qubits))
		self._state = StabilizerStateChForm(vector, len(self._qubits))

	def state_vector(self) -> list[complex]:
		return self._state.state_vector()

	def apply_operation(self, operation: Any) -> None:
		raise NotImplementedError("CH-form simulation does not apply gates in this build")


class SimulationProductState:
	"""A tensor product of independent simulation states."""

	def __init__(self, states: Sequence[SimulationStateBase]) -> None:
		self._states = tuple(states)

	@property
	def states(self) -> tuple[SimulationStateBase, ...]:
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
		return self._clifford_state


class StabilizerSampler:
	"""A sampler backed by the Clifford simulator."""

	def __init__(self, *, seed: int | None = None) -> None:
		self._seed = seed

	def run(self, program: Any, repetitions: int = 1, **kwargs: Any) -> Any:
		from ._clifford import CliffordSimulator

		return CliffordSimulator(seed=self._seed).run(program, repetitions, **kwargs)

	def sample(self, program: Any, repetitions: int = 1, **kwargs: Any) -> Any:
		return self.run(program, repetitions, **kwargs)

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
