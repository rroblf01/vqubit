"""Tests for the extra cirq.devices helpers."""

from __future__ import annotations

import vqubit


def test_no_noise_model() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	moment = vqubit.Moment([vqubit.H(q0)])
	assert vqubit.NO_NOISE.noisy_moment(moment, [q0, q1]) == moment


def test_constant_qubit_noise_model() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	model = vqubit.ConstantQubitNoiseModel(vqubit.X)
	result = model.noisy_moment(vqubit.Moment([vqubit.H(q0)]), [q1])
	assert any(operation.gate == vqubit.X for operation in result.operations)


def test_op_identifier() -> None:
	identifier = vqubit.OpIdentifier(vqubit.CNOTGate, 2)
	assert identifier.gate_type is vqubit.CNOTGate
	assert identifier.num_qubits == 2


def test_noise_properties() -> None:
	properties = vqubit.SuperconductingQubitsNoiseProperties(gate_error_1q=0.002, gate_error_2q=0.02)
	assert properties.gate_error(vqubit.OpIdentifier(vqubit.XPowGate, 1)) == 0.002
	assert properties.gate_error(vqubit.OpIdentifier(vqubit.CNOTGate, 2)) == 0.02
	assert isinstance(vqubit.NoiseModelFromNoiseProperties(properties), vqubit.NoiseModel)


def test_grid_device_metadata() -> None:
	q0, q1 = vqubit.GridQubit(0, 0), vqubit.GridQubit(0, 1)
	metadata = vqubit.GridDeviceMetadata([q0, q1], [(q0, q1)])
	assert metadata.qubit_set is not None and q0 in metadata.qubit_set
	assert frozenset((q0, q1)) in metadata.qubit_pairs


def test_unconstrained_device() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	vqubit.UNCONSTRAINED_DEVICE.validate_operation(vqubit.CNOT(q0, q1))


def test_line_topology() -> None:
	topology = vqubit.LineTopology(3)
	assert sorted(topology.nodes) == [0, 1, 2]
	assert len(topology.edges) == 2
	assert topology[1] == [0, 2]


def test_tilted_square_lattice() -> None:
	lattice = vqubit.TiltedSquareLattice(3, 3)
	assert len(lattice.nodes) == 9
	assert len(lattice.edges) == 8


def test_get_placements_and_validity() -> None:
	topology = vqubit.LineTopology(3)
	graph = {0: {1}, 1: {0, 2}, 2: {1}}
	placements = vqubit.get_placements(topology, graph)
	assert len(placements) >= 1
	device = vqubit.LineDevice(3)
	device_graph = {q: set(device.neighbors(q)) for q in device.qubits}
	physical = [vqubit.LineQubit(0), vqubit.LineQubit(1), vqubit.LineQubit(2)]
	assert vqubit.is_valid_placement(graph, physical, device_graph)


def test_draw_helpers() -> None:
	assert isinstance(vqubit.draw_gridlike({0: {1}, 1: {0}}), str)
	assert isinstance(vqubit.draw_placements({0: {1}, 1: {0}}, [{0: 0, 1: 1}]), str)
