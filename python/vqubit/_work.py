"""Sampling jobs and collectors, mirroring ``cirq.work``.

These let a driver request samples from a sampler towards some purpose (for
instance estimating the energy of a Hamiltonian) without caring where the
samples come from.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

from ._circuits import Circuit
from ._gates import MeasurementGate, measure
from ._pauli import PauliSum
from ._sim import Result
from ._study import Sampler, to_resolvers
from ._value import value_equality


@value_equality(unhashable=True)
class CircuitSampleJob:
	"""Describes a sampling task."""

	def __init__(self, circuit: Circuit, *, repetitions: int, tag: Any = None) -> None:
		self.circuit = circuit
		self.repetitions = repetitions
		self.tag = tag

	def _value_equality_values_(self) -> Any:
		return (self.circuit, self.repetitions, self.tag)

	def __repr__(self) -> str:
		return f"vqubit.CircuitSampleJob(tag={self.tag!r}, repetitions={self.repetitions!r}, circuit={self.circuit!r})"


def _flatten_jobs(tree: Any) -> Iterator[CircuitSampleJob]:
	if isinstance(tree, CircuitSampleJob):
		yield tree
	elif tree is not None:
		for item in tree:
			yield from _flatten_jobs(item)


class Collector:
	"""Collects data from a sampler, towards some purpose.

	Subclasses override :meth:`next_job` and :meth:`on_job_result`.
	"""

	def next_job(self) -> Any:
		"""Return the next sampling job, or ``None`` when there is none."""
		raise NotImplementedError

	def on_job_result(self, job: CircuitSampleJob, result: Result) -> None:
		"""Consume the result of a job this collector asked for."""
		raise NotImplementedError

	def collect(self, sampler: Sampler, *, concurrency: int = 2, max_total_samples: int | None = None) -> None:
		"""Collect needed samples from ``sampler``.

		Samples are taken one job at a time (``concurrency`` is Cirq's hint and is
		accepted for parity); ``max_total_samples`` stops dispatching once the
		budget is spent, even if the last batch handed out had jobs to spare.
		"""
		del concurrency
		remaining = float("inf") if max_total_samples is None else max_total_samples
		queued: list[CircuitSampleJob] = []
		while remaining > 0:
			if not queued:
				queued.extend(_flatten_jobs(self.next_job()))
			if not queued:
				break
			job = queued.pop(0)
			remaining -= job.repetitions
			self.on_job_result(job, sampler.run(job.circuit, repetitions=job.repetitions))

	async def collect_async(
		self, sampler: Sampler, *, concurrency: int = 2, max_total_samples: int | None = None
	) -> None:
		"""Asynchronous variant of :meth:`collect` (runs synchronously here)."""
		self.collect(sampler, concurrency=concurrency, max_total_samples=max_total_samples)


def _measurement_shapes(circuit: Circuit) -> dict[str, int]:
	shapes: dict[str, int] = {}
	for operation in circuit.all_operations():
		if isinstance(operation.gate, MeasurementGate):
			key = operation.gate.key
			shapes[key] = max(shapes.get(key, 0), len(operation.qubits))
	return shapes


class ZerosSampler(Sampler):
	"""A mock sampler for testing that immediately returns zeroes."""

	def __init__(self, device: Any = None) -> None:
		super().__init__()
		self.device = device

	def run(
		self,
		program: Circuit,
		repetitions: int = 1,
		*,
		qubit_order: Any = None,
		seed: int | None = None,
	) -> Result:
		"""Return all-zero outcomes, one per measured qubit."""
		if self.device is not None:
			self.device.validate_circuit(program)
		shapes = _measurement_shapes(program)
		records = {key: [[0] * size for _ in range(repetitions)] for key, size in shapes.items()}
		return Result(records)

	def run_sweep(self, program: Circuit, params: Any, repetitions: int = 1) -> list[Result]:
		"""Return all-zero outcomes for each set of parameters."""
		return [self.run(program, repetitions=repetitions) for _ in to_resolvers(params)]


def _circuit_plus_pauli_string_measurements(circuit: Circuit, pauli_string: Any) -> Circuit:
	"""A circuit that measures the given Pauli observable at the end."""
	operations = list(circuit.all_operations())
	operations.extend(pauli_string.to_z_basis_ops())
	operations.append(measure(*sorted(pauli_string.keys()), key="out"))
	return Circuit(operations)


class PauliSumCollector(Collector):
	"""Estimates the energy of a linear combination of Pauli observables."""

	def __init__(
		self,
		circuit: Circuit,
		observable: Any,
		*,
		samples_per_term: int,
		max_samples_per_job: int = 1000000,
	) -> None:
		observable = PauliSum.wrap(observable)
		self._circuit = circuit
		self._samples_per_job = max_samples_per_job
		self._pauli_coef_terms = [(p / p.coefficient, p.coefficient) for p in observable if p]
		self._identity_offset = 0j
		for p in observable:
			if not p:
				self._identity_offset += p.coefficient
		self._zeros: dict[Any, int] = {}
		self._ones: dict[Any, int] = {}
		self._samples_per_term = samples_per_term
		self._total_samples_requested = 0

	def next_job(self) -> CircuitSampleJob | None:
		"""Return the next chunk of samples needed for the energy."""
		index = self._total_samples_requested // self._samples_per_term
		if index >= len(self._pauli_coef_terms):
			return None
		pauli, _ = self._pauli_coef_terms[index]
		remaining = self._samples_per_term * (index + 1) - self._total_samples_requested
		amount = min(remaining, self._samples_per_job)
		self._total_samples_requested += amount
		return CircuitSampleJob(
			circuit=_circuit_plus_pauli_string_measurements(self._circuit, pauli),
			repetitions=amount,
			tag=pauli,
		)

	def on_job_result(self, job: CircuitSampleJob, result: Result) -> None:
		"""Add a chunk of samples to the energy estimate."""
		parities = result.histogram(key="out", fold_func=lambda bits: sum(bits) % 2)
		self._zeros[job.tag] = self._zeros.get(job.tag, 0) + parities.get(0, 0)
		self._ones[job.tag] = self._ones.get(job.tag, 0) + parities.get(1, 0)

	def estimated_energy(self) -> float | complex:
		"""Sum up the sampled expectations, weighted by their coefficients."""
		energy = 0j
		for pauli_string, coefficient in self._pauli_coef_terms:
			a = self._zeros.get(pauli_string, 0)
			b = self._ones.get(pauli_string, 0)
			if a + b:
				energy += coefficient * (a - b) / (a + b)
		energy = complex(energy) + self._identity_offset
		return energy.real if energy.imag == 0 else energy


__all__ = [
	"CircuitSampleJob",
	"Collector",
	"PauliSumCollector",
	"ZerosSampler",
]
