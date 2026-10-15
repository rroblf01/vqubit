"""Qubit routing and mapping.

These map a circuit's logical qubits onto the physical qubits of a device,
inserting ``SWAP`` gates where a two-qubit gate's operands are not adjacent.
"""

from __future__ import annotations

import itertools
from collections import deque
from collections.abc import Mapping, Sequence
from typing import Any

from ._circuits import Circuit
from ._devices import Device
from ._gates import SWAP
from ._ops import Operation
from ._qubits import Qid


class AbstractInitialMapper:
	"""Base class for deciding the initial logical-to-physical mapping."""

	def __init__(self, *, device: Device | None = None, **kwargs: Any) -> None:
		del kwargs
		self.device = device

	def initial_mapping(self, circuit: Circuit) -> dict[Qid, Qid]:
		raise NotImplementedError

	def map_operations(self, operation: Operation, mapping: Mapping[Qid, Qid]) -> Operation:
		return operation.with_qubits(*[mapping[q] for q in operation.qubits])

	def update_mapping(self, operations: Sequence[Operation], mapping: dict[Qid, Qid]) -> dict[Qid, Qid]:
		return mapping


class HardCodedInitialMapper(AbstractInitialMapper):
	"""An initial mapper that uses a fixed logical-to-physical mapping."""

	def __init__(self, mapping: Mapping[Qid, Qid], *, device: Device | None = None) -> None:
		super().__init__(device=device)
		self._mapping = dict(mapping)

	def initial_mapping(self, circuit: Circuit) -> dict[Qid, Qid]:
		return dict(self._mapping)


class LineInitialMapper(AbstractInitialMapper):
	"""Places logical qubits onto the device's qubits in order of appearance."""

	def initial_mapping(self, circuit: Circuit) -> dict[Qid, Qid]:
		if self.device is None:
			raise ValueError("LineInitialMapper needs a device")
		logical = list(circuit.all_qubits())
		physical = sorted(self.device.qubits)
		if len(logical) > len(physical):
			raise ValueError("the circuit uses more qubits than the device has")
		return {logical_qubit: physical[index] for index, logical_qubit in enumerate(logical)}


class MappingManager:
	"""Manages the mapping of a circuit's qubits onto a device."""

	def __init__(self, device: Device) -> None:
		self.device = device
		self._mapping: dict[Qid, Qid] = {}

	@property
	def current_mapping(self) -> dict[Qid, Qid]:
		return dict(self._mapping)

	def initial_mapping(self, circuit: Circuit) -> dict[Qid, Qid]:
		self._mapping = LineInitialMapper(device=self.device).initial_mapping(circuit)
		return dict(self._mapping)

	def map_circuit(self, circuit: Circuit) -> Circuit:
		mapped, self._mapping = routed_circuit_with_mapping(circuit, LineInitialMapper(device=self.device), self.device)
		return mapped


class RouteCQC:
	"""Routes a circuit onto a device's connectivity, inserting SWAPs."""

	def __init__(self, device: Device) -> None:
		self.device = device

	def route_circuit(
		self, circuit: Circuit, *, initial_mapper: AbstractInitialMapper | None = None, **kwargs: Any
	) -> Circuit:
		mapped, _ = routed_circuit_with_mapping(circuit, initial_mapper, self.device, **kwargs)
		return mapped


def _shortest_path(device: Device, start: Qid, goal: Qid) -> list[Qid]:
	queue: deque[list[Qid]] = deque([[start]])
	visited = {start}
	while queue:
		path = queue.popleft()
		node = path[-1]
		if node == goal:
			return path
		for neighbor in device.neighbors(node):
			if neighbor not in visited:
				visited.add(neighbor)
				queue.append([*path, neighbor])
	return [start]


def routed_circuit_with_mapping(
	circuit: Circuit,
	initial_mapper: AbstractInitialMapper | None = None,
	device: Device | None = None,
	**kwargs: Any,
) -> tuple[Circuit, dict[Qid, Qid]]:
	"""Route ``circuit`` onto ``device`` and return it with the final mapping."""
	del kwargs
	if device is None and initial_mapper is not None:
		device = initial_mapper.device
	if device is None:
		raise ValueError("routed_circuit_with_mapping needs a device")
	mapper = initial_mapper if initial_mapper is not None else LineInitialMapper(device=device)
	mapping = mapper.initial_mapping(circuit)
	result = Circuit()
	for operation in circuit.all_operations():
		if len(operation.qubits) < 2:
			result.append(operation.with_qubits(*[mapping[q] for q in operation.qubits]))
			continue
		physical = [mapping[q] for q in operation.qubits]
		while not _all_adjacent(device, physical):
			path = _shortest_path(device, physical[-1], physical[0])
			if len(path) < 2:
				break
			first, second = path[0], path[1]
			result.append(SWAP(first, second))
			for logical, physical_qubit in list(mapping.items()):
				if physical_qubit == first:
					mapping[logical] = second
				elif physical_qubit == second:
					mapping[logical] = first
			physical = [mapping[q] for q in operation.qubits]
		result.append(operation.with_qubits(*[mapping[q] for q in operation.qubits]))
	return result, mapping


def _all_adjacent(device: Device, qubits: Sequence[Qid]) -> bool:
	return all(second in device.neighbors(first) for first, second in itertools.pairwise(qubits))


__all__ = [
	"AbstractInitialMapper",
	"HardCodedInitialMapper",
	"LineInitialMapper",
	"MappingManager",
	"RouteCQC",
	"routed_circuit_with_mapping",
]
