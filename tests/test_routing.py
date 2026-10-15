"""Tests for qubit routing and mapping."""

from __future__ import annotations

import vqubit


def test_line_initial_mapper() -> None:
	device = vqubit.LineDevice(3)
	q0, q1 = vqubit.LineQubit.range(2)
	mapping = vqubit.LineInitialMapper(device=device).initial_mapping(vqubit.Circuit(vqubit.CNOT(q0, q1)))
	assert len(mapping) == 2
	assert set(mapping.values()) <= set(device.qubits)


def test_hard_coded_initial_mapper() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	p0, p1 = vqubit.LineQubit(10), vqubit.LineQubit(11)
	mapper = vqubit.HardCodedInitialMapper({q0: p0, q1: p1})
	assert mapper.initial_mapping(vqubit.Circuit(vqubit.CNOT(q0, q1))) == {q0: p0, q1: p1}


def test_routed_circuit_inserts_swaps() -> None:
	device = vqubit.LineDevice(3)
	q0, q1, q2 = vqubit.LineQubit.range(3)
	# All three qubits are used, so q0 and q2 land on distant physical qubits.
	circuit = vqubit.Circuit(vqubit.X(q1), vqubit.CNOT(q0, q2))
	mapped, _ = vqubit.routed_circuit_with_mapping(circuit, vqubit.LineInitialMapper(device=device), device)
	assert len(list(mapped.all_operations())) >= 2
	for operation in mapped.all_operations():
		if len(operation.qubits) == 2:
			a, b = operation.qubits
			assert b in device.neighbors(a) or a in device.neighbors(b)


def test_route_cqc_and_mapping_manager() -> None:
	device = vqubit.LineDevice(3)
	q0, q1, q2 = vqubit.LineQubit.range(3)
	circuit = vqubit.Circuit(vqubit.CNOT(q0, q2), vqubit.CNOT(q1, q2))
	routed = vqubit.RouteCQC(device).route_circuit(circuit)
	for operation in routed.all_operations():
		if len(operation.qubits) == 2:
			a, b = operation.qubits
			assert b in device.neighbors(a) or a in device.neighbors(b)
	manager = vqubit.MappingManager(device)
	assert manager.map_circuit(circuit) is not None
	assert len(manager.current_mapping) == 3
