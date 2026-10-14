"""A local engine.

Vendor engines (`cirq_google.Engine`, `cirq_ionq.Service`) submit circuits to
hardware over the network. This engine keeps the same shape but runs everything
locally on the simulators, so a script written against the API works offline
with no credentials.
"""

from __future__ import annotations

from collections.abc import Sequence

from ._circuits import Circuit
from ._devices import Device
from ._parameter import ParamResolver
from ._protocols import resolve_parameters
from ._qasm import qasm as _to_qasm
from ._qubits import Qid
from ._sim import Result, Simulator
from ._study import Points, Sampler


class Engine:
	"""A local stand-in for a vendor engine."""

	def __init__(self, *, device: Device | None = None, seed: int | None = None, project_id: str | None = None) -> None:
		self.device = device
		self.project_id = project_id
		self._seed = seed
		self._simulator = Simulator(seed=seed)

	def _prepare(self, program: Circuit, device: Device | None, param_resolver: ParamResolver | None) -> Circuit:
		resolved = resolve_parameters(program, param_resolver) if param_resolver is not None else program
		target = device or self.device
		if target is not None:
			target.validate_circuit(resolved)
		return resolved

	def run(
		self,
		program: Circuit,
		repetitions: int = 1,
		*,
		param_resolver: ParamResolver | None = None,
		device: Device | None = None,
		seed: int | None = None,
	) -> Result:
		"""Run ``program`` locally and return a :class:`~vqubit.Result`."""
		prepared = self._prepare(program, device, param_resolver)
		return self._simulator.run(prepared, repetitions, seed=seed)

	sample = run

	def run_sweep(
		self,
		program: Circuit,
		params: Points | ParamResolver,
		repetitions: int = 1,
		*,
		device: Device | None = None,
		seed: int | None = None,
	) -> list[Result]:
		"""Run ``program`` once per point of a parameter sweep."""
		if isinstance(params, ParamResolver):
			resolvers = [params]
		elif hasattr(params, "param_resolvers"):
			resolvers = list(params.param_resolvers)
		else:
			raise TypeError(f"cannot sweep over {params!r}")
		return [
			self.run(program, repetitions, param_resolver=resolver, device=device, seed=seed) for resolver in resolvers
		]

	def get_sampler(self, device: Device | None = None) -> Sampler:
		"""Return a :class:`~vqubit.Sampler` for this engine."""
		return Sampler(seed=self._seed)

	def qasm(self, program: Circuit) -> str:
		"""Return the OpenQASM 2.0 for ``program``."""
		return _to_qasm(program)


def get_engine(*, device: Device | None = None, seed: int | None = None, project_id: str | None = None) -> Engine:
	"""Return a local :class:`Engine`."""
	return Engine(device=device, seed=seed, project_id=project_id)


def sampler(
	program: Circuit, repetitions: int = 1, *, qubit_order: Sequence[Qid] | None = None, seed: int | None = None
) -> Result:
	"""Run ``program`` on a local sampler (a thin alias of ``vqubit.sample``)."""
	return Simulator(seed=seed).run(program, repetitions, qubit_order=qubit_order)


__all__ = ["Engine", "get_engine", "sampler"]
