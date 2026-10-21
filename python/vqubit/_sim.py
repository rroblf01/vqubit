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

import math
import random
from bisect import bisect_right
from collections import Counter
from collections.abc import Callable, Iterator, Sequence
from typing import Any

import vqubit_native as _native

from ._circuits import Circuit
from ._compat import HAS_NUMPY, is_numpy_array, native_complex_view
from ._gates import H, MeasurementGate, S, measure
from ._linalg import to_flat, vector_from_flat, vector_to_flat
from ._ops import Moment
from ._protocols import is_parameterized, resolve_parameters
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


def _as_resolver(param_resolver: Any) -> Any:
	"""Return ``param_resolver`` as a :class:`ParamResolver` (an empty one for ``None``)."""
	from ._parameter import ParamResolver

	if param_resolver is None:
		return ParamResolver({})
	if isinstance(param_resolver, ParamResolver):
		return param_resolver
	return ParamResolver(param_resolver)


def _resolve_program(program: Circuit, param_resolver: Any) -> Circuit:
	"""Apply ``param_resolver`` to ``program`` and refuse to leave symbols behind.

	Cirq does the same: a simulator handed a parameterized circuit without a
	resolver raises, instead of failing later on some gate whose unitary is
	``None``.
	"""
	if param_resolver is not None:
		program = resolve_parameters(program, param_resolver)
	if is_parameterized(program):
		raise ValueError("the circuit still contains unsolved symbols; pass param_resolver= to fix them")
	return program


def _qubit_pauli(qubit: Qid, letter: str) -> Any:
	"""A one-qubit ``PauliString`` naming ``qubit``, for the Bloch vector."""
	from ._pauli import Pauli, PauliString

	return PauliString({qubit: Pauli(letter)})


def _reduced_density_matrix(state_vector: Any, qubit_map: Any, qubits: Sequence[Qid] | None) -> Any:
	"""The density matrix of ``qubits``, the rest traced out."""
	from ._qis import partial_trace
	from ._qis_extra import density_matrix_from_state_vector

	full = density_matrix_from_state_vector(list(state_vector))
	keep = [qubit_map[qubit] for qubit in qubits] if qubits is not None else list(range(len(qubit_map)))
	return partial_trace(full, sorted(keep))


def _copy_records(records: dict[str, list[list[int]]]) -> dict[str, list[list[int]]]:
	"""A defensive copy of measurement records, so a step does not keep growing."""
	return {key: [list(row) for row in rows] for key, rows in records.items()}


def _as_observables(observables: Any) -> list[Any]:
	"""The observables as a list of ``PauliSum``, accepting one or many."""
	from ._pauli import PauliSum

	if isinstance(observables, list):
		return [PauliSum.wrap(observable) for observable in observables]
	return [PauliSum.wrap(observables)]


class _ExpectationValues:
	"""Exact expectation values of Pauli observables, as Cirq's simulators give.

	Everything is written on top of ``simulate_sweep_iter`` and the per-simulator
	:meth:`_observable_expectation`, the way Cirq's ``SimulatesExpectationValues``
	is written on top of ``simulate_expectation_values_sweep_iter``.
	"""

	simulate_sweep_iter: Any

	def _observable_expectation(self, result: Any, observable: Any) -> complex:
		"""How this simulator reads an observable off a trial result."""
		raise NotImplementedError

	def simulate_expectation_values_sweep_iter(
		self,
		program: Circuit,
		observables: Any,
		params: Any,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
		permit_terminal_measurements: bool = False,
	) -> Iterator[list[complex]]:
		"""Exact expectation values at every point of the sweep ``params``."""
		if not permit_terminal_measurements and _has_terminal_measurements(program):
			raise ValueError(
				"the circuit has terminal measurements, which may skew expectation values; "
				"pass permit_terminal_measurements=True if that is what you want"
			)
		sums = _as_observables(observables)
		for result in self.simulate_sweep_iter(program, params, qubit_order, initial_state):
			yield [self._observable_expectation(result, observable) for observable in sums]

	def simulate_expectation_values_sweep(
		self,
		program: Circuit,
		observables: Any,
		params: Any,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
		permit_terminal_measurements: bool = False,
	) -> list[list[complex]]:
		"""Like :meth:`simulate_expectation_values_sweep_iter`, but a list."""
		return list(
			self.simulate_expectation_values_sweep_iter(
				program, observables, params, qubit_order, initial_state, permit_terminal_measurements
			)
		)

	def simulate_expectation_values(
		self,
		program: Circuit,
		observables: Any,
		param_resolver: Any = None,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
		permit_terminal_measurements: bool = False,
	) -> list[complex]:
		"""Exact expectation values of ``observables`` on the final state."""
		return self.simulate_expectation_values_sweep(
			program, observables, param_resolver, qubit_order, initial_state, permit_terminal_measurements
		)[0]


class _SampledExpectationValues:
	"""Expectation values estimated by sampling, as Cirq's samplers give.

	Each Pauli string is measured in its own basis (rotating it onto ``Z``) for
	``num_samples`` repetitions, and the parity of the measured bits averages to
	its expectation value -- which is what a real device (or Cirq's
	``Sampler.sample_expectation_values``) does.
	"""

	run: Any

	def sample_expectation_values(
		self,
		program: Circuit,
		observables: Any,
		*,
		num_samples: int,
		params: Any = None,
		permit_terminal_measurements: bool = False,
	) -> list[list[float]]:
		"""Estimate expectation values from ``num_samples`` repetitions."""
		if num_samples <= 0:
			raise ValueError(f"expectation values need at least one sample, got {num_samples}")
		sums = _as_observables(observables)
		if not sums:
			raise ValueError("at least one observable is needed")
		if not permit_terminal_measurements and _has_terminal_measurements(program):
			raise ValueError(
				"the circuit has terminal measurements, which may skew expectation values; "
				"pass permit_terminal_measurements=True if that is what you want"
			)
		from ._study import to_resolvers

		estimates: list[list[float]] = []
		for resolver in to_resolvers(params):
			row: list[float] = []
			for observable in sums:
				total = 0.0
				for term in observable:
					total += _sample_pauli_string(self, program, term, num_samples, resolver)
				row.append(total)
			estimates.append(row)
		return estimates


def _sample_pauli_string(
	simulator: Any, program: Circuit, pauli_string: Any, num_samples: int, param_resolver: Any
) -> float:
	"""Estimate ``<P>`` by measuring ``P`` in its own basis.

	The basis change is written with the stabilizer kernels' own gates (``X``
	needs a Hadamard, ``Y`` a Hadamard after an ``S†``; both verified against the
	unitaries), because those kernels can only apply the instructions they know,
	and ``PauliString.to_z_basis_ops`` hands out general Clifford gates.
	"""
	from ._parameter import ParamResolver
	from ._pauli import Pauli

	if not pauli_string:
		return pauli_string.coefficient.real
	stages: list[list[Any]] = []
	for qubit, pauli in pauli_string.items():
		if pauli is Pauli.Z:
			continue
		gates = (H,) if pauli is Pauli.X else (S**-1, H)
		for depth, gate in enumerate(gates):
			stages.extend([] for _ in range(depth + 1 - len(stages)))
			stages[depth].append(gate.on(qubit))
	measured = Circuit(
		*program.moments,
		*[Moment(stage) for stage in stages if stage],
		Moment([measure(*pauli_string.qubits, key="m")]),
	)
	result = simulator.run(measured, ParamResolver(param_resolver), num_samples)
	rows = result.measurements["m"]
	signs = sum(1 if sum(row) % 2 == 0 else -1 for row in rows)
	return pauli_string.coefficient.real * signs / num_samples


class _SweepsAndBatches:
	"""The sweep and batch entry points the simulators share.

	Each concrete simulator implements ``run``, ``simulate`` and (where it can)
	``simulate_moment_steps``; everything here is written in terms of those, the
	way Cirq's ``SimulatorBase`` is.
	"""

	#: Provided by the concrete simulators (annotated like ``_sim_abc`` does, so
	#: their own signatures stay free to be the Cirq ones).
	run: Any
	simulate: Any

	def run_sweep(self, program: Circuit, params: Any, repetitions: int = 1) -> list[Result]:
		"""Run ``program`` once per point of the parameter sweep ``params``."""
		return list(self.run_sweep_iter(program, params, repetitions))

	def run_sweep_iter(self, program: Circuit, params: Any, repetitions: int = 1) -> Iterator[Result]:
		"""Like :meth:`run_sweep`, but yielding each result as it is produced."""
		from ._study import to_resolvers

		for resolver in to_resolvers(params):
			yield self.run(program, resolver, repetitions)

	def simulate_sweep(
		self, program: Circuit, params: Any, qubit_order: Sequence[Qid] | None = None, initial_state: Any = 0
	) -> list[Any]:
		"""Like :meth:`simulate`, once per point of the sweep ``params``."""
		return list(self.simulate_sweep_iter(program, params, qubit_order, initial_state))

	def simulate_sweep_iter(
		self, program: Circuit, params: Any, qubit_order: Sequence[Qid] | None = None, initial_state: Any = 0
	) -> Iterator[Any]:
		"""Like :meth:`simulate_sweep`, but yielding each result as it is produced."""
		from ._study import to_resolvers

		for resolver in to_resolvers(params):
			yield self.simulate(program, resolver, qubit_order, initial_state)

	def run_batch(
		self,
		programs: Sequence[Circuit],
		params_list: Sequence[Any] | None = None,
		repetitions: int | Sequence[int] = 1,
	) -> list[list[Result]]:
		"""Run each circuit of ``programs``, optionally with its own sweep, as Cirq does."""
		if params_list is None:
			params_list = [None] * len(programs)
		if isinstance(repetitions, int):
			repetitions = [repetitions] * len(programs)
		if not (len(programs) == len(params_list) == len(repetitions)):
			raise ValueError("programs, params_list and repetitions must have the same length")
		return [
			self.run_sweep(program, params if params is not None else None, count)
			for program, params, count in zip(programs, params_list, repetitions, strict=True)
		]

	def run_async(self, program: Circuit, param_resolver: Any = None, repetitions: int = 1) -> Any:
		"""The coroutine form of :meth:`run`, as Cirq's simulators offer."""
		import asyncio

		return asyncio.sleep(0, result=self.run(program, param_resolver, repetitions))

	def run_sweep_async(self, program: Circuit, params: Any, repetitions: int = 1) -> Any:
		"""The coroutine form of :meth:`run_sweep`."""
		import asyncio

		return asyncio.sleep(0, result=self.run_sweep(program, params, repetitions))

	def run_batch_async(
		self,
		programs: Sequence[Circuit],
		params_list: Sequence[Any] | None = None,
		repetitions: int | Sequence[int] = 1,
	) -> Any:
		"""The coroutine form of :meth:`run_batch`."""
		import asyncio

		return asyncio.sleep(0, result=self.run_batch(programs, params_list, repetitions))


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
		"""The qubits of the state, in order."""
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

	def __init__(
		self,
		qubits: Sequence[Qid],
		state: StateVector,
		params: Any = None,
		measurements: dict[str, list[list[int]]] | None = None,
	) -> None:
		self._qubits = tuple(qubits)
		self._final_state = state
		self._params = _as_resolver(params)
		self._measurements = measurements if measurements is not None else {}

	@property
	def params(self) -> Any:
		"""The resolver used for the simulation (empty when there was none)."""
		return self._params

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits of the final state, in order."""
		return self._qubits

	@property
	def qubit_map(self) -> dict[Qid, int]:
		"""Which index each qubit holds in the state, as Cirq reports it."""
		return {qubit: index for index, qubit in enumerate(self._qubits)}

	def bloch_vector_of(self, qubit: Qid) -> list[float]:
		"""The Bloch vector of ``qubit`` in the final state: ``[<X>, <Y>, <Z>]``."""
		return [
			_qubit_pauli(qubit, letter)
			.expectation_from_state_vector(self.final_state_vector, self.qubit_map, check_preconditions=False)
			.real
			for letter in "XYZ"
		]

	def density_matrix_of(self, qubits: Sequence[Qid] | None = None) -> Any:
		"""The reduced density matrix of ``qubits`` (every qubit by default)."""
		return _reduced_density_matrix(self.final_state_vector, self.qubit_map, qubits)

	def get_state_containing_qubit(self, qubit: Qid) -> Any:
		"""A simulation state holding ``qubit``, as Cirq's non-split simulators return.

		vqubit keeps one state for the whole circuit, so this is that state (with
		the qubit in it), not a factor of it.
		"""
		from ._sim_step import StateVectorSimulationState

		del qubit
		return StateVectorSimulationState(qubits=self.qubits, initial_state=self.final_state_vector)

	@property
	def measurements(self) -> dict[str, list[list[int]]]:
		"""The measurement outcomes the simulation sampled (empty when none)."""
		return self._measurements

	@property
	def final_state(self) -> StateVector:
		"""The state at the end of the simulation."""
		return self._final_state

	@property
	def final_state_vector(self) -> Any:
		"""The state vector, under the name Cirq gives it."""
		return self._final_state.state_vector()

	def state_vector(self) -> Any:
		"""The final state as a flat vector of amplitudes."""
		return self._final_state.state_vector()

	def dirac_notation(self, decimals: int = 3) -> str:
		"""The final state written in Dirac notation."""
		return self._final_state.dirac_notation(decimals)

	def __str__(self) -> str:
		return self._final_state.dirac_notation()

	def __repr__(self) -> str:
		return f"vqubit.SimulationResult({self._final_state!r})"


class Result:
	"""Measurement results: one row of bits per repetition, per key.

	Carries the ``params`` of the run that produced it -- the resolver that fixed
	the circuit's symbols (empty when there were none) -- as Cirq's ``Result``
	does, which is how a caller tells the points of a sweep apart.
	"""

	def __init__(self, measurements: dict[str, list[list[int]]], params: Any = None) -> None:
		self._measurements = {key: [list(row) for row in rows] for key, rows in measurements.items()}
		self._params = _as_resolver(params)

	@property
	def params(self) -> Any:
		"""The resolver used to run the circuit."""
		return self._params

	@property
	def measurements(self) -> dict[str, list[list[int]]]:
		"""The measurement outcomes, keyed by measurement name."""
		return self._measurements

	def __len__(self) -> int:
		return len(next(iter(self._measurements.values()), []))

	def histogram(self, key: str = "", fold_func: Callable[[list[int]], Any] | None = None) -> Counter:
		"""Count how often each outcome appears for ``key``.

		>>> import vqubit
		>>> q0 = vqubit.LineQubit(0)
		>>> circuit = vqubit.Circuit(vqubit.X(q0), vqubit.measure(q0, key="m"))
		>>> vqubit.Simulator().run(circuit, repetitions=3).histogram(key="m")
		Counter({1: 3})
		"""
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


def _has_terminal_measurements(program: Circuit) -> bool:
	"""Whether any measurement is left un-followed on its qubits, as Cirq checks.

	Cirq's ``are_any_measurements_terminal``: a measurement at the end of the
	circuit usually means the caller wanted samples, not expectation values.
	"""
	operations = list(program.all_operations())
	for index, operation in enumerate(operations):
		if not isinstance(operation.gate, MeasurementGate):
			continue
		qubits = set(operation.qubits)
		if not any(qubits & set(later.qubits) for later in operations[index + 1 :]):
			return True
	return False


def _measurements_are_terminal(program: Circuit) -> bool:
	"""Report whether every measurement happens after all non-measurements.

	When that holds, the final state can be computed once and sampled directly
	instead of re-simulating the circuit for every repetition.
	"""
	operations = list(program.all_operations())
	for index, operation in enumerate(operations):
		# A *gate* measurement, not a sub-circuit that happens to contain one:
		# that one is not terminal, it has to run.
		if not isinstance(_operation_gate(operation), MeasurementGate):
			continue
		qubits = set(operation.qubits)
		for later in operations[index + 1 :]:
			later_gate = _operation_gate(later)
			if not isinstance(later_gate, MeasurementGate) and qubits.intersection(later.qubits):
				return False
	return True


def _operation_gate(operation: Any) -> Any | None:
	"""The gate of ``operation``, or ``None`` when it does not wrap one.

	A ``CircuitOperation`` or a ``TaggedOperation`` has no gate of its own (the
	property raises), and the simulator has to take those apart instead of asking
	them for a matrix.
	"""
	try:
		return operation.gate
	except NotImplementedError:
		return None


def _is_unitary_circuit(program: Circuit) -> bool:
	"""Report whether every operation has a unitary.

	A channel needs one trajectory per repetition, so the "simulate once and
	sample the final distribution" shortcut only applies when this holds.
	"""
	for operation in program.all_operations():
		gate = _operation_gate(operation)
		if isinstance(gate, MeasurementGate):
			continue
		if gate is None or gate._cached_flat_unitary_() is None:
			# A composite operation or a channel: its state depends on the shot.
			return False
	return True


def _kraus_matrices(gate: Any) -> list[list[list[complex]]] | None:
	"""The Kraus operators of ``gate``, or ``None`` when it is not a channel."""
	getter = getattr(gate, "_kraus_", None)
	if getter is None:
		return None
	operators = getter()
	if operators is None:
		return None
	return [[[complex(value) for value in row] for row in operator] for operator in operators]


def _diagonal_effects(effects: Sequence[Sequence[Sequence[complex]]], tolerance: float = 1e-12) -> bool:
	"""Whether every ``K†K`` is diagonal, so the probabilities need marginals only."""
	for effect in effects:
		size = len(effect)
		for row in range(size):
			for column in range(size):
				if row != column and abs(effect[row][column]) > tolerance:
					return False
	return True


def _scatter_bits(rest: int, target: int, targets: Sequence[int], num_qubits: int) -> int:
	"""The full state index for a "rest" index and a target-subsystem index, MSB first."""
	rest_qubits = [qubit for qubit in range(num_qubits) if qubit not in targets]
	index = 0
	for position, qubit in enumerate(rest_qubits):
		bit = (rest >> (len(rest_qubits) - 1 - position)) & 1
		index |= bit << (num_qubits - 1 - qubit)
	for position, qubit in enumerate(targets):
		bit = (target >> (len(targets) - 1 - position)) & 1
		index |= bit << (num_qubits - 1 - qubit)
	return index


def _target_density_matrix(amplitudes: Sequence[float], num_qubits: int, targets: Sequence[int]) -> list[list[complex]]:
	"""The reduced density matrix of ``targets``, from the interleaved amplitudes."""
	size = 1 << len(targets)
	rest_size = 1 << (num_qubits - len(targets))
	reduced = [[0j] * size for _ in range(size)]
	for rest in range(rest_size):
		row = [complex(0, 0)] * size
		indices = [_scatter_bits(rest, target, targets, num_qubits) for target in range(size)]
		for target in range(size):
			index = 2 * indices[target]
			row[target] = complex(amplitudes[index], amplitudes[index + 1])
		for first in range(size):
			for second in range(size):
				reduced[first][second] += row[first] * row[second].conjugate()
	return reduced


def _channel_probabilities(
	state: Any, effects: Sequence[Sequence[Sequence[complex]]], targets: Sequence[int]
) -> list[float]:
	"""The probability of each Kraus branch, given the current state."""
	size = len(effects[0])
	if _diagonal_effects(effects):
		# Every channel this library defines has diagonal effects, so the native
		# marginals are enough: no pass over the whole state in Python.
		marginal = list(state.probabilities(list(targets)))
		probabilities = []
		for effect in effects:
			total = 0.0
			for index in range(size):
				total += float(effect[index][index].real) * marginal[index]
			probabilities.append(max(total, 0.0))
		return probabilities
	# A general Kraus channel: its effects are not diagonal, so the coherences of
	# the target subsystem matter. One pass builds its reduced density matrix.
	amplitudes = list(state.amplitudes())
	num_qubits = len(amplitudes).bit_length() - 2
	reduced = _target_density_matrix(amplitudes, num_qubits, list(targets))
	probabilities = []
	for effect in effects:
		accumulated = 0j
		for row in range(size):
			for column in range(size):
				accumulated += effect[row][column] * reduced[column][row]
		probabilities.append(max(accumulated.real, 0.0))
	return probabilities


def _apply_channel(state: Any, gate: Any, targets: Sequence[int], rng: random.Random) -> None:
	"""Apply one trajectory of a channel to the state."""
	from ._channels import ResetChannel

	if isinstance(gate, ResetChannel):
		for target in targets:
			if state.measure([target], rng.random()):
				state.apply_gate([0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0], [target])
		return
	operators = _kraus_matrices(gate)
	if not operators:
		raise TypeError(f"{gate!r} has neither a unitary nor Kraus operators and cannot be simulated")
	effects = []
	for operator in operators:
		size = len(operator)
		conjugate_transpose = [[operator[column][row].conjugate() for column in range(size)] for row in range(size)]
		effects.append(
			[
				[sum(conjugate_transpose[row][k] * operator[k][column] for k in range(size)) for column in range(size)]
				for row in range(size)
			]
		)
	probabilities = _channel_probabilities(state, effects, targets)
	total = sum(probabilities)
	if total <= 0:
		raise ValueError(f"{gate!r} has no Kraus branch with positive probability on this state")
	threshold = rng.random() * total
	accumulated = 0.0
	chosen = len(probabilities) - 1
	for index, probability in enumerate(probabilities):
		accumulated += probability
		if threshold < accumulated:
			chosen = index
			break
	scale = 1.0 / math.sqrt(probabilities[chosen])
	operator = operators[chosen]
	size = len(operator)
	flat: list[float] = []
	for row in range(size):
		for column in range(size):
			value = operator[row][column] * scale
			flat.extend((value.real, value.imag))
	state.apply_gate(flat, list(targets))


def _draw_outcome(probabilities: Sequence[float], rng: random.Random) -> int:
	threshold = rng.random()
	accumulated = 0.0
	for index, probability in enumerate(probabilities):
		accumulated += probability
		if threshold < accumulated:
			return index
	return len(probabilities) - 1


def _cumulative_probabilities(probabilities: Sequence[float]) -> list[float]:
	"""The running totals of ``probabilities``, in the same summation order."""
	cumulative: list[float] = []
	accumulated = 0.0
	for probability in probabilities:
		accumulated += probability
		cumulative.append(accumulated)
	return cumulative


def _draw_cumulative(cumulative: Sequence[float], rng: random.Random) -> int:
	"""Draw an outcome from precomputed running totals.

	The same draw as :func:`_draw_outcome` —first total strictly above the
	threshold, or the last outcome— without walking the distribution on every
	shot: ``bisect_right`` skips the totals that are not above the threshold.
	"""
	return min(bisect_right(cumulative, rng.random()), len(cumulative) - 1)


def _outcome_of(targets: Sequence[int], positions: dict[int, int], union_outcome: int) -> int:
	"""Read one measurement's outcome out of a joint outcome over the union.

	The bits are most-significant-first, the same convention ``_draw_outcome``
	and ``NativeState.measure`` use for a single measurement.
	"""
	size = len(positions)
	outcome = 0
	for target in targets:
		outcome = (outcome << 1) | ((union_outcome >> (size - 1 - positions[target])) & 1)
	return outcome


class Simulator(
	_SweepsAndBatches,
	_ExpectationValues,
	_SampledExpectationValues,
	SimulatorBase,
	SimulatesSamples,
	SimulatesAmplitudes,
	SimulatesFinalState,
	SimulatesIntermediateState,
):
	"""A state-vector simulator.

	``simulate`` evolves a circuit with no measurements and returns the final
	state. ``run`` samples measurements and returns a :class:`Result`.
	"""

	def __init__(self, *, seed: int | None = None) -> None:
		self._seed = seed

	def simulate(
		self,
		program: Circuit,
		param_resolver: Any = None,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
	) -> SimulationResult:
		"""Simulate a circuit and return its final state.

		>>> import vqubit
		>>> q0 = vqubit.LineQubit(0)
		>>> result = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.H(q0)))
		>>> [round(float(abs(z)), 3) for z in result.state_vector()]
		[0.707, 0.707]
		"""
		program = _resolve_program(program, param_resolver)
		qubits = _resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		num_qubits = len(qubits)
		state = _native.NativeState(num_qubits)
		if _is_basis_zero(initial_state):
			state.set_zero_state()
		else:
			state.set_amplitudes(vector_to_flat(_initial_state_vector(initial_state, 1 << num_qubits)))
		# A measurement is sampled and collapses the state, as in Cirq's simulator.
		records: dict[str, list[list[int]]] = {}
		for moment in program.moments:
			self._apply_moment(state, moment, index, rng=random.Random(self._seed), measurements=records)
		view = native_complex_view(state)
		if view is not None:
			return SimulationResult(qubits, StateVector(view, qubits, owner=state), param_resolver, records)
		return SimulationResult(
			qubits, StateVector(vector_from_flat(state.amplitudes()), qubits), param_resolver, records
		)

	def simulate_moment_steps(
		self,
		program: Circuit,
		param_resolver: Any = None,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
	) -> Iterator[Any]:
		"""Yield the state after each moment, as Cirq's simulator does."""
		from ._sim_step import SparseSimulatorStep

		program = _resolve_program(program, param_resolver)
		qubits = _resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		state = _native.NativeState(len(qubits))
		if _is_basis_zero(initial_state):
			state.set_zero_state()
		else:
			state.set_amplitudes(vector_to_flat(_initial_state_vector(initial_state, 1 << len(qubits))))
		records: dict[str, list[list[int]]] = {}
		for moment in program.moments:
			self._apply_moment(state, moment, index, rng=random.Random(self._seed), measurements=records)
			yield SparseSimulatorStep(
				qubits=qubits, measurements=_copy_records(records), state_vector=vector_from_flat(state.amplitudes())
			)

	def compute_amplitudes(
		self,
		program: Circuit,
		bitstrings: Sequence[int],
		param_resolver: Any = None,
		qubit_order: Sequence[Qid] | None = None,
	) -> list[complex]:
		"""The amplitudes of the given basis states in the final state."""
		if is_numpy_array(bitstrings):
			shape = getattr(bitstrings, "shape", ())
			if len(shape) > 1:  # pragma: no cover - numpy extra
				raise ValueError(f"bitstrings must be a 1-dimensional array, got shape {shape}")
		state = self.simulate(program, param_resolver, qubit_order).state_vector()
		return [complex(state[index]) for index in bitstrings]

	def run(
		self,
		program: Circuit,
		param_resolver: Any = None,
		repetitions: int = 1,
		*,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
		seed: int | None = None,
	) -> Result:
		"""Run a circuit, sampling its measurements."""
		program = _resolve_program(program, param_resolver)
		if repetitions < 1:
			raise ValueError("repetitions must be at least 1")
		qubits = _resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		num_qubits = len(qubits)
		dim = 1 << num_qubits
		amplitudes = vector_to_flat(_initial_state_vector(initial_state, dim))
		rng = random.Random(self._seed if seed is None else seed)
		records: dict[str, list[list[int]]] = {}
		if _measurements_are_terminal(program) and _is_unitary_circuit(program):
			state = _native.NativeState(num_qubits)
			if _is_basis_zero(initial_state):
				state.set_zero_state()
			else:
				state.set_amplitudes(amplitudes)
			measurements: list[tuple[MeasurementGate, list[int]]] = []
			for moment in program.moments:
				for operation in moment.operations:
					gate = _operation_gate(operation)
					if isinstance(gate, MeasurementGate):
						measurements.append((gate, [index[q] for q in operation.qubits]))
						continue
					if gate is None or gate._cached_flat_unitary_() is None:
						raise TypeError(f"{operation!r} has no unitary and cannot be simulated")
					state.apply_gate(gate._cached_flat_unitary_(), [index[q] for q in operation.qubits])
			if not measurements:
				return Result(records, param_resolver)
			# Every measurement is terminal, so their joint distribution is the one
			# of the union of measured qubits. Sampling each measurement separately
			# from its own marginal would decorrelate them: two qubits of a Bell
			# state measured with two ``measure`` calls must still agree.
			union: list[int] = []
			for _, targets in measurements:
				for target in targets:
					if target not in union:
						union.append(target)
			positions = {target: position for position, target in enumerate(union)}
			probabilities = list(state.probabilities(union))
			if HAS_NUMPY:
				import numpy as _np

				np_rng = _np.random.default_rng(self._seed if seed is None else seed)
				cumulative = _np.cumsum(_np.asarray(probabilities, dtype=_np.float64))
				for outcome in _np.searchsorted(cumulative, np_rng.random(repetitions)).tolist():
					for gate, targets in measurements:
						bits = _outcome_of(targets, positions, int(outcome))
						_record_bits(gate, bits, len(targets), records)
				return Result(records, param_resolver)
			# The numpy branch above keeps an ``ndarray`` in ``cumulative``; the
			# pure-Python one uses the running totals, built once, not per shot.
			running = _cumulative_probabilities(probabilities)
			for _ in range(repetitions):
				outcome = _draw_cumulative(running, rng)
				for gate, targets in measurements:
					bits = _outcome_of(targets, positions, outcome)
					_record_bits(gate, bits, len(targets), records)
			return Result(records, param_resolver)
		for _ in range(repetitions):
			state = _native.NativeState(num_qubits)
			state.set_amplitudes(amplitudes)
			for moment in program.moments:
				self._apply_moment(state, moment, index, rng=rng, measurements=records)
		return Result(records, param_resolver)

	sample = run

	def _observable_expectation(self, result: SimulationResult, observable: Any) -> complex:
		"""Read an observable off the final state vector of a simulation."""
		return observable.expectation_from_state_vector(
			result.final_state_vector, result.qubit_map, check_preconditions=False
		)

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
			from ._circuits import CircuitOperation

			if isinstance(operation, CircuitOperation):
				for _ in range(operation.repetitions):
					for inner in operation.circuit.moments:
						self._apply_moment(state, inner, index, rng=rng, measurements=measurements)
				continue
			gate = _operation_gate(operation)
			targets = [index[q] for q in operation.qubits]
			if isinstance(gate, MeasurementGate):
				if rng is None or measurements is None:
					raise ValueError("simulate() cannot run a circuit with measurements; use run()")
				outcome = state.measure(targets, rng.random())
				_record_bits(gate, outcome, len(targets), measurements)
				continue
			if gate is None:
				raise TypeError(f"{operation!r} has no unitary and cannot be simulated")
			flat = gate._cached_flat_unitary_()
			if flat is None:
				if rng is None:
					raise ValueError(f"{gate!r} is a channel; simulate() needs a seed to pick a trajectory")
				_apply_channel(state, gate, targets, rng)
				continue
			state.apply_gate(flat, targets)


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
		"""The qubits of the density matrix, in order."""
		return self._qubits

	def density_matrix(self) -> Any:
		"""Return the density matrix (a zero-copy view when numpy is installed)."""
		if is_numpy_array(self._matrix):
			return self._matrix
		return [list(row) for row in self._matrix]

	def dirac_notation(self, decimals: int = 3) -> str:
		"""The state written in Dirac notation."""
		return StateVector(
			[self._matrix[i][i].real ** 0.5 + 0j for i in range(len(self._matrix))], self._qubits
		).dirac_notation(decimals)

	def __repr__(self) -> str:
		return f"vqubit.DensityMatrix({len(self._qubits)} qubits)"


class DensityMatrixTrialResult(SimulationTrialResultBase):
	"""The result of :meth:`DensityMatrixSimulator.simulate`."""

	def __init__(
		self,
		qubits: Sequence[Qid],
		density: DensityMatrix,
		params: Any = None,
		measurements: dict[str, list[list[int]]] | None = None,
	) -> None:
		self._qubits = tuple(qubits)
		self._final_density_matrix = density
		self._params = _as_resolver(params)
		self._measurements = measurements if measurements is not None else {}

	@property
	def params(self) -> Any:
		"""The resolver used for the simulation (empty when there was none)."""
		return self._params

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits of the final state, in order."""
		return self._qubits

	@property
	def qubit_map(self) -> dict[Qid, int]:
		"""Which index each qubit holds in the density matrix, as Cirq reports it."""
		return {qubit: index for index, qubit in enumerate(self._qubits)}

	@property
	def measurements(self) -> dict[str, list[list[int]]]:
		"""The measurement outcomes the simulation sampled (empty when none)."""
		return self._measurements

	@property
	def final_density_matrix(self) -> Any:
		"""The density matrix at the end of the simulation.

		An array (a nested list without numpy), like Cirq's; the qubits it is
		expressed in are :attr:`qubits` and :attr:`qubit_map`.
		"""
		return self._final_density_matrix.density_matrix()

	def density_matrix(self) -> Any:
		"""The final density matrix (the same value as the property)."""
		return self._final_density_matrix.density_matrix()

	def dirac_notation(self, decimals: int = 3) -> str:
		"""The diagonal of the density matrix written in Dirac notation."""
		return self._final_density_matrix.dirac_notation(decimals)

	def __repr__(self) -> str:
		return f"vqubit.DensityMatrixTrialResult({self._final_density_matrix!r})"


def _density_flat_to_nested(flat: Sequence[float], dim: int) -> list[list[complex]]:
	return [[complex(flat[2 * (i * dim + j)], flat[2 * (i * dim + j) + 1]) for j in range(dim)] for i in range(dim)]


class DensityMatrixSimulator(
	_SweepsAndBatches,
	_ExpectationValues,
	_SampledExpectationValues,
	SimulatorBase,
	SimulatesSamples,
	SimulatesFinalState,
	SimulatesIntermediateState,
):
	"""A simulator that evolves a density matrix, so it can run noisy circuits."""

	def __init__(self, *, seed: int | None = None) -> None:
		self._seed = seed

	def simulate(
		self,
		program: Circuit,
		param_resolver: Any = None,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
	) -> DensityMatrixTrialResult:
		"""Simulate a circuit and return its final density matrix."""
		program = _resolve_program(program, param_resolver)
		qubits = _resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		num_qubits = len(qubits)
		density = _native.NativeDensity(num_qubits)
		if _is_basis_zero(initial_state):
			density.set_zero_state()
		else:
			density.set_from_amplitudes(vector_to_flat(_initial_state_vector(initial_state, 1 << num_qubits)))
		# A measurement is sampled and collapses the state, as in Cirq's simulator.
		records: dict[str, list[list[int]]] = {}
		for moment in program.moments:
			self._apply_moment(density, moment, index, rng=random.Random(self._seed), measurements=records)
		view = native_complex_view(density)
		if view is not None:
			matrix = DensityMatrix(view.reshape(1 << num_qubits, 1 << num_qubits), qubits, owner=density)
		else:
			matrix = DensityMatrix(_density_flat_to_nested(density.density_matrix(), 1 << num_qubits), qubits)
		return DensityMatrixTrialResult(qubits, matrix, param_resolver, records)

	def simulate_moment_steps(
		self,
		program: Circuit,
		param_resolver: Any = None,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
	) -> Iterator[Any]:
		"""Yield the density matrix after each moment, as Cirq's simulator does."""
		from ._sim_step import DensityMatrixStepResult

		program = _resolve_program(program, param_resolver)
		qubits = _resolve_qubits(program, qubit_order)
		index = {q: i for i, q in enumerate(qubits)}
		num_qubits = len(qubits)
		dim = 1 << num_qubits
		density = _native.NativeDensity(num_qubits)
		if _is_basis_zero(initial_state):
			density.set_zero_state()
		else:
			density.set_from_amplitudes(vector_to_flat(_initial_state_vector(initial_state, dim)))
		records: dict[str, list[list[int]]] = {}
		for moment in program.moments:
			self._apply_moment(density, moment, index, rng=random.Random(self._seed), measurements=records)
			yield DensityMatrixStepResult(
				qubits=qubits,
				measurements=_copy_records(records),
				density_matrix=_density_flat_to_nested(density.density_matrix(), dim),
			)

	def run(
		self,
		program: Circuit,
		param_resolver: Any = None,
		repetitions: int = 1,
		*,
		qubit_order: Sequence[Qid] | None = None,
		initial_state: Any = 0,
		seed: int | None = None,
	) -> Result:
		"""Run a circuit, sampling its measurements."""
		program = _resolve_program(program, param_resolver)
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
		return Result(records, param_resolver)

	sample = run

	def _observable_expectation(self, result: DensityMatrixTrialResult, observable: Any) -> complex:
		"""Read an observable off the final density matrix of a simulation."""
		return observable.expectation_from_density_matrix(
			result.density_matrix(), result.qubit_map, check_preconditions=False
		)

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
			from ._circuits import CircuitOperation

			if isinstance(operation, CircuitOperation):
				for _ in range(operation.repetitions):
					for inner in operation.circuit.moments:
						self._apply_moment(density, inner, index, rng=rng, measurements=measurements)
				continue
			gate = _operation_gate(operation)
			targets = [index[q] for q in operation.qubits]
			if isinstance(gate, MeasurementGate):
				if rng is None or measurements is None:
					raise ValueError("simulate() cannot run a circuit with measurements; use run()")
				outcome = density.measure(targets, rng.random())
				_record_bits(gate, outcome, len(targets), measurements)
				continue
			if gate is None:
				raise TypeError(f"{operation!r} has neither a unitary nor Kraus operators")
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
