"""Abstract simulator interfaces, mirroring Cirq's ``cirq.sim`` protocols.

These are the marker/base classes the concrete simulators and simulation states
derive from. They carry no behaviour beyond the documented method names.
"""

from __future__ import annotations

from typing import Any

from ._qubits import Qid

CIRCUIT_LIKE = Any


class SimulatorBase:
	"""Base class for a simulator."""

	simulate: Any
	run: Any
	sample: Any


class SimulatesSamples:
	"""A simulator that can sample measurement outcomes."""

	run: Any
	sample: Any


class SimulatesAmplitudes:
	"""A simulator that can compute state amplitudes."""

	compute_amplitudes: Any


class SimulatesExpectationValues:
	"""A simulator that can compute expectation values of observables."""

	simulate_expectation_values: Any


class SimulatesFinalState:
	"""A simulator that can return the final state."""

	simulate: Any


class SimulatesIntermediateState:
	"""A simulator that can return the state after each moment."""

	simulate: Any


class SimulatesIntermediateStateVector(SimulatesIntermediateState):
	"""A simulator returning intermediate state vectors."""


class ClassicalStateSimulator:
	"""A simulator that tracks only classical (measurement) state."""

	simulate: Any


class SimulationStateBase:
	"""Base class for the mutable state of a step-by-step simulation."""

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits the state lives on."""
		raise NotImplementedError


class StepResultBase:
	"""Base class for the result of simulating a single step."""

	@property
	def qubits(self) -> tuple[Qid, ...]:
		"""The qubits of the state, in order."""
		raise NotImplementedError


class SimulationTrialResultBase:
	"""Base class for the result of a full simulation trial."""

	@property
	def params(self) -> Any:
		"""The parameter resolver the simulation used."""
		raise NotImplementedError


__all__ = [
	"CIRCUIT_LIKE",
	"ClassicalStateSimulator",
	"SimulatesAmplitudes",
	"SimulatesExpectationValues",
	"SimulatesFinalState",
	"SimulatesIntermediateState",
	"SimulatesIntermediateStateVector",
	"SimulatesSamples",
	"SimulationStateBase",
	"SimulationTrialResultBase",
	"SimulatorBase",
	"StepResultBase",
]
