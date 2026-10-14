"""Running circuits and sweeping parameters.

The helpers here mirror Cirq's study module at its useful core: sampling a
circuit with a local simulator, and running it once per point of a parameter
sweep.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from typing import Any

from ._circuits import Circuit
from ._parameter import Parameter, ParamResolver
from ._protocols import resolve_parameters
from ._qubits import Qid
from ._sim import Result, Simulator


class Points:
	"""A sweep over an explicit list of parameter assignments."""

	def __init__(self, points: Mapping[Any, Any] | Sequence[Mapping[Any, Any]]) -> None:
		if isinstance(points, Mapping):
			self._points = [dict(points)]
		else:
			self._points = [dict(point) for point in points]
		if not self._points:
			raise ValueError("Points needs at least one parameter assignment")

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		return [ParamResolver(point) for point in self._points]

	def __iter__(self) -> Iterator[ParamResolver]:
		return iter(self.param_resolvers)

	def __len__(self) -> int:
		return len(self._points)

	def __repr__(self) -> str:
		return f"vqubit.Points({self._points!r})"


def _resolver_dict(value: Any) -> dict[str, float]:
	if isinstance(value, ParamResolver):
		return {key: float(item) for key, item in value.items()}
	if isinstance(value, Mapping):
		return {(key.name if isinstance(key, Parameter) else str(key)): float(item) for key, item in value.items()}
	raise TypeError(f"cannot read parameter values from {value!r}")


class Linspace:
	"""A sweep of evenly spaced values for a set of parameters.

	Like Cirq's ``Linspace``, the values start at ``start`` and stop before
	``stop``.
	"""

	def __init__(self, start: Any, stop: Any, length: int) -> None:
		if length < 1:
			raise ValueError("Linspace needs a length of at least 1")
		self._start = _resolver_dict(start)
		self._stop = _resolver_dict(stop)
		if self._start.keys() != self._stop.keys():
			raise ValueError("Linspace start and stop must name the same parameters")
		self._length = int(length)

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		resolvers: list[ParamResolver] = []
		for i in range(self._length):
			values = {
				name: self._start[name] + (self._stop[name] - self._start[name]) * i / self._length
				for name in self._start
			}
			resolvers.append(ParamResolver(values))
		return resolvers

	def __iter__(self) -> Iterator[ParamResolver]:
		return iter(self.param_resolvers)

	def __len__(self) -> int:
		return self._length

	def __repr__(self) -> str:
		return f"vqubit.Linspace({self._start!r}, {self._stop!r}, {self._length})"


def sample(
	program: Circuit,
	repetitions: int = 1,
	*,
	qubit_order: Sequence[Qid] | None = None,
	seed: int | None = None,
) -> Result:
	"""Sample ``program`` with the state-vector simulator."""
	return Simulator(seed=seed).run(program, repetitions, qubit_order=qubit_order)


def run(
	program: Circuit,
	repetitions: int = 1,
	*,
	qubit_order: Sequence[Qid] | None = None,
	seed: int | None = None,
) -> Result:
	"""Run ``program`` and return a :class:`~vqubit.Result`.

	This is an alias of :func:`sample`, matching Cirq, where running and
	sampling a circuit are the same thing.
	"""
	return sample(program, repetitions, qubit_order=qubit_order, seed=seed)


def run_sweep(
	program: Circuit,
	params: Points | Linspace | ParamResolver,
	repetitions: int = 1,
	*,
	qubit_order: Sequence[Qid] | None = None,
	seed: int | None = None,
) -> list[Result]:
	"""Run ``program`` once for each point of a parameter sweep."""
	if isinstance(params, ParamResolver):
		resolvers = [params]
	elif hasattr(params, "param_resolvers"):
		resolvers = list(params.param_resolvers)
	else:
		raise TypeError(f"cannot sweep over {params!r}")
	results: list[Result] = []
	for resolver in resolvers:
		resolved = resolve_parameters(program, resolver)
		results.append(sample(resolved, repetitions, qubit_order=qubit_order, seed=seed))
	return results


class Sampler:
	"""A thin sampler around a local simulator."""

	def __init__(self, *, seed: int | None = None) -> None:
		self._simulator = Simulator(seed=seed)

	def run(
		self,
		program: Circuit,
		repetitions: int = 1,
		*,
		qubit_order: Sequence[Qid] | None = None,
		seed: int | None = None,
	) -> Result:
		return self._simulator.run(program, repetitions, qubit_order=qubit_order, seed=seed)

	sample = run


__all__ = ["Linspace", "Points", "Sampler", "run", "run_sweep", "sample"]
