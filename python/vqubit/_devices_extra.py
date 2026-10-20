"""Extra ``cirq.devices`` helpers: noise models, device metadata, topologies
and placement routines.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from ._circuits import Moment
from ._devices import Device, DeviceMetadata
from ._gates import Gate
from ._ops import Operation
from ._qubits import GridQubit, Qid
from ._value import value_equality

NOISE_MODEL_LIKE = Any


@value_equality
class OpIdentifier:
	"""Identifies a kind of operation by its gate type and qubit count."""

	def __init__(self, gate_type: type, num_qubits: int | None = None) -> None:
		self.gate_type = gate_type
		self.num_qubits = num_qubits

	def _value_equality_values_(self) -> Any:
		return (self.gate_type, self.num_qubits)

	def __repr__(self) -> str:
		return f"vqubit.OpIdentifier({self.gate_type.__name__}, num_qubits={self.num_qubits})"


class NoiseProperties:
	"""Base class describing a device's per-operation noise."""

	def gate_error(self, op_id: OpIdentifier) -> float:
		"""The error rate of one kind of gate."""
		raise NotImplementedError

	def readout_error(self, qubit: Qid) -> float:
		"""The readout error of one qubit."""
		raise NotImplementedError

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}()"


class SuperconductingQubitsNoiseProperties(NoiseProperties):
	"""Noise properties for a superconducting device."""

	def __init__(
		self,
		*,
		gate_times_ns: Mapping[OpIdentifier, float] | None = None,
		gate_error_1q: float = 0.001,
		gate_error_2q: float = 0.01,
		readout_error_rate: float = 0.02,
	) -> None:
		self.gate_times_ns = dict(gate_times_ns) if gate_times_ns is not None else {}
		self.gate_error_1q = float(gate_error_1q)
		self.gate_error_2q = float(gate_error_2q)
		self.readout_error_rate = float(readout_error_rate)

	def gate_error(self, op_id: OpIdentifier) -> float:
		"""The error rate the device declares for that gate."""
		if op_id.num_qubits is not None and op_id.num_qubits >= 2:
			return self.gate_error_2q
		return self.gate_error_1q

	def readout_error(self, qubit: Qid) -> float:
		"""The readout error the device declares for that qubit."""
		del qubit
		return self.readout_error_rate


class NoiseModel:
	"""Base class for a noise model."""

	def noisy_moment(self, moment: Moment, qubits: Sequence[Qid]) -> Moment | list[Moment]:
		"""The moment with the noise added: one moment, or a list of them."""
		raise NotImplementedError

	def noisy_moments(self, moments: Iterable[Moment], qubits: Sequence[Qid]) -> list[Moment]:
		"""Apply :meth:`noisy_moment` to each moment, flattening what it returns.

		A model may answer with one moment or with a list (the moment plus the
		noise that follows it), so the result is flattened here rather than
		leaving a nested list to the caller.
		"""
		result: list[Moment] = []
		for moment in moments:
			noisy = self.noisy_moment(moment, qubits)
			if isinstance(noisy, Moment):
				result.append(noisy)
			else:
				result.extend(noisy)
		return result

	def noisy_operation(self, operation: Operation) -> Operation:
		"""Return the operation with the model's noise applied."""
		return operation

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}()"


class ConstantQubitNoiseModel(NoiseModel):
	"""Applies a fixed noise gate to every qubit in every moment."""

	def __init__(self, qubit_noise_gate: Gate) -> None:
		self.qubit_noise_gate = qubit_noise_gate

	def noisy_moment(self, moment: Moment, qubits: Sequence[Qid]) -> list[Moment]:
		"""The moment followed by one noise gate on every qubit.

		The noise goes in its own moment: a qubit that the moment already uses
		cannot take a second operation in the same instant, and dropping the
		noise in that case would silently un-noise the circuit.
		"""
		return [moment, Moment(self.qubit_noise_gate.on(qubit) for qubit in qubits)]


class _NoNoiseModel(NoiseModel):
	"""The noise model that adds no noise."""

	def noisy_moment(self, moment: Moment, qubits: Sequence[Qid]) -> Moment:
		return moment


NO_NOISE = _NoNoiseModel()


class NoiseModelFromNoiseProperties(NoiseModel):
	"""A noise model derived from a device's noise properties."""

	def __init__(self, noise_properties: NoiseProperties) -> None:
		self.noise_properties = noise_properties

	def noisy_operation(self, operation: Operation) -> Operation:
		"""Return the operation with the device's noise applied."""
		return operation

	def noisy_moment(self, moment: Moment, qubits: Sequence[Qid]) -> list[Moment]:
		"""The moment followed by a depolarizing channel per operation.

		Each operation gets the device's gate error for its kind (one- or
		two-qubit), applied on the same qubits. Readout error is a measurement
		effect, not a gate, so it is left to the sampler.
		"""
		from ._channels import depolarize

		noise: list[Operation] = []
		for operation in moment.operations:
			error = self.noise_properties.gate_error(OpIdentifier(type(operation.gate), operation.num_qubits()))
			if error > 0:
				noise.append(depolarize(error, num_qubits=operation.num_qubits()).on(*operation.qubits))
		if not noise:
			return [moment]
		return [moment, Moment(noise)]


class GridDeviceMetadata(DeviceMetadata):
	"""Metadata for a grid device: its qubits and its allowed qubit pairs."""

	def __init__(
		self,
		qubits: Iterable[Qid],
		qubit_pairs: Iterable[tuple[Qid, Qid]],
		gate_durations: Mapping[OpIdentifier, Any] | None = None,
	) -> None:
		super().__init__(
			name="grid", qubit_set=qubits, gate_durations=dict(gate_durations) if gate_durations is not None else {}
		)
		self.qubit_pairs = frozenset(frozenset(pair) for pair in qubit_pairs)

	def __repr__(self) -> str:
		return f"vqubit.GridDeviceMetadata({len(self.qubit_set or [])} qubits)"


class _UnconstrainedDevice(Device):
	"""A device with no connectivity constraints."""

	def __init__(self) -> None:
		super().__init__([])

	def neighbors(self, qubit: Qid) -> list[Qid]:
		return []

	def validate_operation(self, operation: Operation) -> None:
		return None


UNCONSTRAINED_DEVICE = _UnconstrainedDevice()


class NamedTopology:
	"""A named qubit topology, described by its nodes and edges."""

	def __init__(self, name: str) -> None:
		self.name = name

	@property
	def nodes(self) -> frozenset[Any]:
		"""The qubits of the topology."""
		return frozenset()

	@property
	def edges(self) -> frozenset[frozenset[Any]]:
		"""The pairs of adjacent qubits."""
		return frozenset()

	@property
	def graph(self) -> dict[Any, set[Any]]:
		"""The adjacency of the topology, keyed by qubit."""
		adjacency: dict[Any, set[Any]] = {node: set() for node in self.nodes}
		for edge in self.edges:
			first, second = tuple(edge)
			adjacency[first].add(second)
			adjacency[second].add(first)
		return adjacency

	def __getitem__(self, node: Any) -> list[Any]:
		return list(self.graph.get(node, set()))

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}({self.name!r})"


class LineTopology(NamedTopology):
	"""A line topology with ``length`` nodes."""

	def __init__(self, length: int) -> None:
		super().__init__("line")
		self.length = int(length)

	@property
	def nodes(self) -> frozenset[Any]:
		"""The qubits of the topology."""
		return frozenset(range(self.length))

	@property
	def edges(self) -> frozenset[frozenset[Any]]:
		"""The pairs of adjacent qubits."""
		return frozenset(frozenset((i, i + 1)) for i in range(self.length - 1))


class TiltedSquareLattice(NamedTopology):
	"""A tilted square lattice of ``width`` by ``height`` grid qubits."""

	def __init__(self, width: int, height: int) -> None:
		super().__init__("tilted_square_lattice")
		self.width = int(width)
		self.height = int(height)

	@property
	def nodes(self) -> frozenset[Any]:
		"""The qubits of the lattice."""
		return frozenset(GridQubit(row, col) for row in range(self.height) for col in range(self.width))

	def _is_node(self, row: int, col: int) -> bool:
		return 0 <= row < self.height and 0 <= col < self.width

	@property
	def edges(self) -> frozenset[frozenset[Any]]:
		"""The pairs of adjacent qubits."""
		edges: set[frozenset[Any]] = set()
		for row in range(self.height):
			for col in range(self.width):
				for d_row, d_col in ((1, 1), (1, -1)):
					if self._is_node(row + d_row, col + d_col):
						edges.add(frozenset((GridQubit(row, col), GridQubit(row + d_row, col + d_col))))
		return frozenset(edges)


def _to_adjacency(graph: Any) -> dict[Any, set[Any]]:
	if isinstance(graph, NamedTopology):
		return graph.graph
	if isinstance(graph, Mapping):
		return {node: set(neighbors) for node, neighbors in graph.items()}
	if hasattr(graph, "nodes") and hasattr(graph, "neighbors"):
		return {node: set(graph.neighbors(node)) for node in graph.nodes()}
	raise TypeError(f"cannot read a graph from {graph!r}")


def is_valid_placement(
	logical_graph: Any,
	physical_qubits: Sequence[Qid],
	device_graph: Any,
) -> bool:
	"""Return whether ``physical_qubits`` places ``logical_graph`` validly."""
	adjacency = _to_adjacency(device_graph)
	nodes = list(_to_adjacency(logical_graph))
	if len(nodes) != len(physical_qubits):
		return False
	mapping = dict(zip(nodes, physical_qubits, strict=True))
	for node, neighbors in _to_adjacency(logical_graph).items():
		for neighbor in neighbors:
			if mapping[neighbor] not in adjacency.get(mapping[node], set()):
				return False
	return True


def get_placements(
	topology: NamedTopology,
	graph: Any,
	*,
	max_placements: int = 10,
	**kwargs: Any,
) -> list[dict[Any, Any]]:
	"""Find placements of ``graph`` onto ``topology`` (greedy backtracking)."""
	del kwargs
	adjacency = _to_adjacency(graph)
	nodes = list(adjacency)
	topology_nodes = list(topology.nodes)
	topology_adjacency = topology.graph
	if len(nodes) > len(topology_nodes):
		return []
	placements: list[dict[Any, Any]] = []

	def backtrack(index: int, mapping: dict[Any, Any], used: set[Any]) -> None:
		if len(placements) >= max_placements:
			return
		if index == len(nodes):
			placements.append(dict(mapping))
			return
		node = nodes[index]
		placed_neighbors = [mapping[n] for n in adjacency[node] if n in mapping]
		candidates = topology_nodes
		if placed_neighbors:
			common: set[Any] | None = None
			for placed in placed_neighbors:
				neighbors = topology_adjacency.get(placed, set())
				common = set(neighbors) if common is None else common & neighbors
			candidates = list(common or [])
		for candidate in candidates:
			if candidate in used:
				continue
			mapping[node] = candidate
			used.add(candidate)
			backtrack(index + 1, mapping, used)
			del mapping[node]
			used.discard(candidate)

	backtrack(0, {}, set())
	return placements


def draw_gridlike(graph: Any) -> str:
	"""Return a simple text drawing of a gridlike graph."""
	adjacency = _to_adjacency(graph)
	lines = []
	for node in sorted(adjacency, key=str):
		neighbors = sorted(adjacency[node], key=str)
		lines.append(f"{node}: {', '.join(str(n) for n in neighbors)}")
	return "\n".join(lines)


def draw_placements(
	logical_graph: Any,
	placements: Sequence[Mapping[Any, Any]],
	*,
	device_graph: Any = None,
	**kwargs: Any,
) -> str:
	"""Return a simple text drawing of the given placements."""
	del kwargs
	adjacency = _to_adjacency(logical_graph)
	lines = []
	for index, placement in enumerate(placements):
		pairs = []
		for node, neighbors in adjacency.items():
			for neighbor in neighbors:
				if node < neighbor:
					pairs.append(f"{placement[node]}-{placement[neighbor]}")
		lines.append(f"placement {index}: {', '.join(pairs)}")
	return "\n".join(lines)


__all__ = [
	"NOISE_MODEL_LIKE",
	"NO_NOISE",
	"UNCONSTRAINED_DEVICE",
	"ConstantQubitNoiseModel",
	"GridDeviceMetadata",
	"LineTopology",
	"NamedTopology",
	"NoiseModel",
	"NoiseModelFromNoiseProperties",
	"NoiseProperties",
	"OpIdentifier",
	"SuperconductingQubitsNoiseProperties",
	"TiltedSquareLattice",
	"draw_gridlike",
	"draw_placements",
	"get_placements",
	"is_valid_placement",
]
