"""Running circuits and sweeping parameters.

The helpers here mirror Cirq's study module at its useful core: sampling a
circuit with a local simulator, the sweep algebra (``SingleSweep``, ``Zip``,
``Product``, ``Concat``), and running a circuit once per point of a sweep.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterator, Mapping, Sequence
from typing import Any

from ._circuits import Circuit
from ._parameter import Parameter, ParamResolver
from ._protocols import resolve_parameters
from ._qubits import Qid
from ._sim import Result, Simulator


class Sweep:
	"""A set of parameter assignments to run a circuit over."""

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		raise NotImplementedError

	def __iter__(self) -> Iterator[ParamResolver]:
		return iter(self.param_resolvers)

	def __len__(self) -> int:
		return len(self.param_resolvers)


def _merge(resolvers: Sequence[ParamResolver]) -> ParamResolver:
	merged: dict[str, Any] = {}
	for resolver in resolvers:
		merged.update(dict(resolver))
	return ParamResolver(merged)


class SingleSweep(Sweep):
	"""A sweep over the values of one parameter."""

	def __init__(self, key: Any, points: Sequence[Any]) -> None:
		self._key = key.name if hasattr(key, "name") else str(key)
		self._points = list(points)

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		return [ParamResolver({self._key: point}) for point in self._points]

	def __len__(self) -> int:
		return len(self._points)

	def __repr__(self) -> str:
		return f"vqubit.SingleSweep({self._key!r}, {self._points!r})"


class _Resolvers(Sweep):
	"""A sweep over an explicit list of resolvers."""

	def __init__(self, resolvers: Sequence[ParamResolver]) -> None:
		self._resolvers = list(resolvers)

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		return list(self._resolvers)


class Zip(Sweep):
	"""Runs several sweeps in step, one point of each at a time."""

	def __init__(self, *sweeps: Any) -> None:
		self._sweeps = [to_sweep(sweep) for sweep in sweeps]
		if not self._sweeps:
			raise ValueError("Zip needs at least one sweep")

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		all_resolvers = [sweep.param_resolvers for sweep in self._sweeps]
		length = len(all_resolvers[0])
		for resolvers in all_resolvers:
			if len(resolvers) != length:
				raise ValueError("Zip sweeps must have the same length")
		return [_merge([resolvers[i] for resolvers in all_resolvers]) for i in range(length)]


class Product(Sweep):
	"""Runs every combination of several sweeps."""

	def __init__(self, *sweeps: Any) -> None:
		self._sweeps = [to_sweep(sweep) for sweep in sweeps]
		if not self._sweeps:
			raise ValueError("Product needs at least one sweep")

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		all_resolvers = [sweep.param_resolvers for sweep in self._sweeps]
		return [_merge(combo) for combo in itertools.product(*all_resolvers)]


class Concat(Sweep):
	"""Runs several sweeps one after another."""

	def __init__(self, *sweeps: Any) -> None:
		self._sweeps = [to_sweep(sweep) for sweep in sweeps]
		if not self._sweeps:
			raise ValueError("Concat needs at least one sweep")

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		resolvers: list[ParamResolver] = []
		for sweep in self._sweeps:
			resolvers.extend(sweep.param_resolvers)
		return resolvers


def to_sweep(value: Any) -> Sweep:
	"""Return ``value`` as a :class:`Sweep`."""
	if isinstance(value, Sweep):
		return value
	if isinstance(value, ParamResolver):
		return _Resolvers([value])
	if isinstance(value, Mapping):
		values = list(value.values())
		if values and all(isinstance(v, Sequence) and not isinstance(v, str) for v in values):
			return Product(*[SingleSweep(key, points) for key, points in value.items()])
		return _Resolvers([ParamResolver(value)])
	raise TypeError(f"cannot read a sweep from {value!r}")


def to_resolvers(value: Any) -> list[ParamResolver]:
	"""Return the resolvers of ``value`` as a sweep."""
	return to_sweep(value).param_resolvers


class Points(Sweep):
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


class Linspace(Sweep):
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

	def __len__(self) -> int:
		return self._length

	def __repr__(self) -> str:
		return f"vqubit.Linspace({self._start!r}, {self._stop!r}, {self._length})"


class TrialResult(Result):
	"""The result of one point of a sweep: measurements and the parameters used."""

	def __init__(self, params: ParamResolver, measurements: dict[str, list[list[int]]]) -> None:
		super().__init__(measurements)
		self._params = params

	@property
	def params(self) -> ParamResolver:
		return self._params

	def __repr__(self) -> str:
		return f"vqubit.TrialResult(params={dict(self._params)!r})"


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
	params: Any,
	repetitions: int = 1,
	*,
	qubit_order: Sequence[Qid] | None = None,
	seed: int | None = None,
) -> list[Result]:
	"""Run ``program`` once for each point of a parameter sweep."""
	results: list[Result] = []
	for resolver in to_resolvers(params):
		resolved = resolve_parameters(program, resolver)
		results.append(sample(resolved, repetitions, qubit_order=qubit_order, seed=seed))
	return results


def run_batch(
	program: Circuit,
	params: Any,
	repetitions: int = 1,
	*,
	qubit_order: Sequence[Qid] | None = None,
	seed: int | None = None,
) -> list[TrialResult]:
	"""Run ``program`` once per point of a sweep, keeping the parameters."""
	results: list[TrialResult] = []
	for resolver in to_resolvers(params):
		resolved = resolve_parameters(program, resolver)
		result = sample(resolved, repetitions, qubit_order=qubit_order, seed=seed)
		results.append(TrialResult(resolver, result.measurements))
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


Sweepable = Any
ParamDictType = Any
ParamMappingType = Any
ParamResolverOrSimilarType = Any


def symbol(name: str) -> Parameter:
	"""Return a parameter symbol with the given name."""
	return Parameter(name)


class ListSweep(_Resolvers):
	"""A sweep over an explicit list of resolvers."""


class UnitSweep(Sweep):
	"""A sweep with a single, empty parameter assignment."""

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		return [ParamResolver({})]

	def __repr__(self) -> str:
		return "vqubit.UnitSweep()"


UNIT_SWEEP = UnitSweep()


class ZipLongest(Sweep):
	"""Runs several sweeps in step, padding the shorter ones with empty points."""

	def __init__(self, *sweeps: Any) -> None:
		self._sweeps = [to_sweep(sweep) for sweep in sweeps]
		if not self._sweeps:
			raise ValueError("ZipLongest needs at least one sweep")

	@property
	def param_resolvers(self) -> list[ParamResolver]:
		all_resolvers = [sweep.param_resolvers for sweep in self._sweeps]
		length = max(len(resolvers) for resolvers in all_resolvers)
		result: list[ParamResolver] = []
		for index in range(length):
			picked = [resolvers[index] if index < len(resolvers) else ParamResolver({}) for resolvers in all_resolvers]
			result.append(_merge(picked))
		return result


def to_sweeps(value: Any) -> list[Sweep]:
	"""Return ``value`` as a list of sweeps."""
	if isinstance(value, Sweep):
		return [value]
	if isinstance(value, Sequence) and not isinstance(value, str):
		return [to_sweep(item) for item in value]
	return [to_sweep(value)]


def flatten(sweepable: Any) -> list[ParamResolver]:
	"""Return the parameter resolvers of a sweepable object."""
	return to_resolvers(sweepable)


def flatten_with_params(sweepable: Any) -> tuple[list[ParamResolver], dict[str, Any]]:
	"""Return the resolvers and the parameter values of a sweepable object."""
	return to_resolvers(sweepable), {}


def flatten_with_sweep(sweepable: Any) -> tuple[list[ParamResolver], Sweep]:
	"""Return the resolvers and the sweep of a sweepable object."""
	return to_resolvers(sweepable), to_sweep(sweepable)


def dict_to_product_sweep(factors: Mapping[Any, Sequence[Any]]) -> Sweep:
	"""Build a product sweep from a mapping of parameters to value lists."""
	return Product(*[SingleSweep(key, points) for key, points in factors.items()])


def dict_to_zip_sweep(factors: Mapping[Any, Sequence[Any]]) -> Sweep:
	"""Build a zip sweep from a mapping of parameters to value lists."""
	return Zip(*[SingleSweep(key, points) for key, points in factors.items()])


def list_of_dicts_to_zip(points: Sequence[Mapping[Any, Any]]) -> Sweep:
	"""Build a zip sweep from a list of parameter-assignment dictionaries."""
	return ListSweep([ParamResolver(point) for point in points])


class ResultDict(dict):
	"""A dictionary of measurement results keyed by measurement key."""


class ExpressionMap(dict):
	"""A mapping from parameter names to their values."""

	def __getitem__(self, key: Any) -> Any:
		return dict.__getitem__(self, key.name if hasattr(key, "name") else key)


__all__ = [
	"UNIT_SWEEP",
	"Concat",
	"ExpressionMap",
	"Linspace",
	"ListSweep",
	"ParamDictType",
	"ParamMappingType",
	"ParamResolverOrSimilarType",
	"Points",
	"Product",
	"ResultDict",
	"Sampler",
	"SingleSweep",
	"Sweep",
	"Sweepable",
	"TrialResult",
	"UnitSweep",
	"Zip",
	"ZipLongest",
	"dict_to_product_sweep",
	"dict_to_zip_sweep",
	"flatten",
	"flatten_with_params",
	"flatten_with_sweep",
	"list_of_dicts_to_zip",
	"run",
	"run_batch",
	"run_sweep",
	"sample",
	"symbol",
	"to_resolvers",
	"to_sweep",
	"to_sweeps",
]
