"""Tests for the sampling jobs and collectors behind ``vqubit.work``.

These let a driver ask a sampler for samples towards a purpose — estimating the
energy of a Hamiltonian, say — without caring where the samples come from.
"""

from __future__ import annotations

import asyncio
from typing import Any

import pytest

import vqubit


def bell_circuit() -> vqubit.Circuit:
	q0, q1 = vqubit.LineQubit.range(2)
	return vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))


def test_circuit_sample_job_is_a_value() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="m"))
	job = vqubit.CircuitSampleJob(circuit, repetitions=10, tag="a")
	assert job == vqubit.CircuitSampleJob(circuit, repetitions=10, tag="a")
	assert job != vqubit.CircuitSampleJob(circuit, repetitions=11, tag="a")
	assert job != vqubit.CircuitSampleJob(circuit, repetitions=10, tag="b")
	assert "tag='a'" in repr(job)


def test_zeros_sampler_returns_zeroes_for_every_key() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(
		vqubit.measure(q0, key="a"),
		vqubit.measure(q0, q1, key="b"),
	)
	sampler = vqubit.ZerosSampler()
	result = sampler.run(circuit, repetitions=3)
	assert result.measurements == {"a": [[0], [0], [0]], "b": [[0, 0], [0, 0], [0, 0]]}
	assert result.params == vqubit.ParamResolver({})


def test_zeros_sampler_run_sweep() -> None:
	circuit = vqubit.Circuit(vqubit.measure(vqubit.LineQubit(0), key="m"))
	results = vqubit.ZerosSampler().run_sweep(circuit, vqubit.Points([{"t": 0}, {"t": 1}]), repetitions=2)
	assert len(results) == 2
	assert all(result.measurements["m"] == [[0], [0]] for result in results)


def test_zeros_sampler_validates_against_a_device() -> None:
	circuit = vqubit.Circuit(vqubit.measure(vqubit.LineQubit(0), key="m"))

	class _StrictDevice:
		def validate_circuit(self, program: Any) -> None:
			raise ValueError("this device does not accept that circuit")

	with pytest.raises(ValueError):
		vqubit.ZerosSampler(device=_StrictDevice()).run(circuit)


def test_the_collector_interface_is_abstract() -> None:
	collector = vqubit.Collector()
	with pytest.raises(NotImplementedError):
		collector.next_job()
	with pytest.raises(NotImplementedError):
		collector.on_job_result(None, None)  # ty: ignore[invalid-argument-type]


class _RecordingCollector(vqubit.Collector):
	"""A collector that hands out a fixed tree of jobs and records the results."""

	def __init__(self, tree: Any) -> None:
		self._tree = tree
		self.jobs: list[vqubit.CircuitSampleJob] = []

	def next_job(self) -> Any:
		tree, self._tree = self._tree, None
		return tree

	def on_job_result(self, job: vqubit.CircuitSampleJob, result: Any) -> None:
		self.jobs.append(job)
		assert result.measurements["m"] == [[0]] * job.repetitions


def test_collector_collects_until_there_are_no_jobs_left() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="m"))
	tree = [
		vqubit.CircuitSampleJob(circuit, repetitions=1),
		[vqubit.CircuitSampleJob(circuit, repetitions=2), None],
	]
	collector = _RecordingCollector(tree)
	collector.collect(vqubit.ZerosSampler())
	assert [job.repetitions for job in collector.jobs] == [1, 2]


def test_collector_stops_at_max_total_samples() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="m"))
	jobs = [vqubit.CircuitSampleJob(circuit, repetitions=5) for _ in range(4)]
	collector = _RecordingCollector(jobs)
	collector.collect(vqubit.ZerosSampler(), max_total_samples=12)
	# The loop stops once the budget is spent: the third job takes it to -3.
	assert len(collector.jobs) == 3


def test_collector_collect_async_does_the_same() -> None:
	q0 = vqubit.LineQubit(0)
	circuit = vqubit.Circuit(vqubit.measure(q0, key="m"))
	collector = _RecordingCollector([vqubit.CircuitSampleJob(circuit, repetitions=1)])
	asyncio.run(collector.collect_async(vqubit.ZerosSampler(), concurrency=4))
	assert len(collector.jobs) == 1


def test_collector_with_no_jobs_does_nothing() -> None:
	collector = _RecordingCollector(None)
	collector.collect(vqubit.ZerosSampler())
	assert collector.jobs == []


def test_pauli_sum_collector_estimates_a_known_energy() -> None:
	# |0> is a +1 eigenstate of Z, so <Z> = 1 and the 0.25 identity offset adds up.
	q0 = vqubit.LineQubit(0)
	observable = 0.5 * vqubit.PauliString({q0: vqubit.Z}) + 0.25
	collector = vqubit.PauliSumCollector(vqubit.Circuit(vqubit.I(q0)), observable, samples_per_term=10)
	collector.collect(vqubit.ZerosSampler())
	assert collector.estimated_energy() == pytest.approx(0.75)


def test_pauli_sum_collector_chunks_large_jobs() -> None:
	q0 = vqubit.LineQubit(0)
	observable = vqubit.PauliString({q0: vqubit.Z})
	collector = vqubit.PauliSumCollector(
		vqubit.Circuit(vqubit.I(q0)), observable, samples_per_term=10, max_samples_per_job=4
	)
	collector.collect(vqubit.ZerosSampler())
	assert collector.estimated_energy() == pytest.approx(1.0)
	assert collector.next_job() is None


def test_pauli_sum_collector_is_zero_before_collecting() -> None:
	q0 = vqubit.LineQubit(0)
	collector = vqubit.PauliSumCollector(
		vqubit.Circuit(vqubit.I(q0)), vqubit.PauliString({q0: vqubit.Z}), samples_per_term=5
	)
	assert collector.estimated_energy() == 0


def test_pauli_sum_collector_measures_each_term_in_its_basis() -> None:
	# X and Y need a basis change before the Z measurement; on the Bell state
	# <XX> = 1 and <YY> = -1, so the energy of (XX - YY) / 2 is 1.
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	observable = (
		vqubit.PauliString({q0: vqubit.X, q1: vqubit.X}) - vqubit.PauliString({q0: vqubit.Y, q1: vqubit.Y})
	) / 2
	collector = vqubit.PauliSumCollector(circuit, observable, samples_per_term=8)
	collector.collect(vqubit.Sampler(seed=3))
	assert collector.estimated_energy() == pytest.approx(1.0)
