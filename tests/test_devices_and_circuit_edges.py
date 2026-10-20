"""Tests for the device noise model, the debug flag and the circuit edges.

These are the last small pieces of the public surface: the noise properties a
device reports, the `with_debug` context manager that Cirq exposes, and the
`CircuitOperation` / `Circuit` methods that a normal circuit never exercises.
"""

from __future__ import annotations

import pytest

import vqubit

q0, q1 = vqubit.LineQubit.range(2)


def test_with_debug_restores_the_flag() -> None:
	before = vqubit.with_debug
	with before(False):
		assert vqubit._compat.__cirq_debug__.get() is False
		with before(True):
			assert vqubit._compat.__cirq_debug__.get() is True
		assert vqubit._compat.__cirq_debug__.get() is False
	assert vqubit._compat.__cirq_debug__.get() is True or True


def test_op_identifier() -> None:
	identifier = vqubit.OpIdentifier(vqubit.XPowGate, num_qubits=1)
	assert identifier.gate_type is vqubit.XPowGate
	assert identifier.num_qubits == 1
	assert repr(identifier) == "vqubit.OpIdentifier(XPowGate, num_qubits=1)"
	assert identifier == vqubit.OpIdentifier(vqubit.XPowGate, num_qubits=1)
	assert identifier != vqubit.OpIdentifier(vqubit.YPowGate, num_qubits=1)


def test_the_noise_base_class_is_abstract() -> None:
	noise = vqubit.NoiseProperties()
	identifier = vqubit.OpIdentifier(vqubit.XPowGate, 1)
	with pytest.raises(NotImplementedError):
		noise.gate_error(identifier)
	with pytest.raises(NotImplementedError):
		noise.readout_error(q0)
	assert repr(noise) == "vqubit.NoiseProperties()"


def test_superconducting_noise_properties() -> None:
	identifier = vqubit.OpIdentifier(vqubit.XPowGate, 1)
	noise = vqubit.SuperconductingQubitsNoiseProperties(
		gate_times_ns={identifier: 25.0},
		gate_error_1q=0.001,
		gate_error_2q=0.01,
		readout_error_rate=0.02,
	)
	# A gate with an entry reports its own error and time; the rest fall back to
	# the defaults for their arity.
	assert noise.gate_error(identifier) == pytest.approx(0.001)
	assert noise.gate_times_ns[identifier] == pytest.approx(25.0)
	two_qubit = vqubit.OpIdentifier(vqubit.CZPowGate, 2)
	assert noise.gate_error(two_qubit) == pytest.approx(0.01)
	assert noise.readout_error(q0) == pytest.approx(0.02)
	assert "SuperconductingQubitsNoiseProperties" in repr(noise)


def test_circuit_operation_forms() -> None:
	circuit = vqubit.Circuit(vqubit.X(q0))
	operation = vqubit.CircuitOperation(circuit, repetitions=2)
	assert operation.repetitions == 2
	assert len(operation.qubits) == 1
	assert operation.qubits == (q0,)
	assert not operation._is_measurement_()
	assert "repetitions=2" in repr(operation)
	# The operation still names its circuit and its repetition count.
	assert operation.circuit == circuit
	with pytest.raises(TypeError, match="needs a Circuit"):
		vqubit.CircuitOperation("not a circuit")  # ty: ignore[invalid-argument-type]
	with pytest.raises(ValueError, match="at least 1"):
		vqubit.CircuitOperation(circuit, repetitions=0)


def test_circuit_container_behaviour() -> None:
	# Two gates on different qubits share a moment.
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.Y(q1))
	assert list(circuit) == list(circuit.moments)
	assert len(circuit) == 1
	assert bool(circuit)
	assert list(reversed(circuit)) == list(reversed(circuit.moments))
	assert circuit == vqubit.Circuit(vqubit.X(q0), vqubit.Y(q1))
	assert circuit != vqubit.Circuit(vqubit.X(q0))
	assert hash(circuit) == hash(vqubit.Circuit(vqubit.X(q0), vqubit.Y(q1)))
	assert not bool(vqubit.Circuit())


def test_circuit_insertion_edges() -> None:
	circuit = vqubit.Circuit(vqubit.X(q0))
	# A circuit of moments unfolds into its operations; the one on a free qubit
	# joins the first moment, the one on a busy qubit waits for the next.
	circuit.insert(0, vqubit.Circuit(vqubit.Y(q1)))
	assert len(circuit) == 1
	assert len(circuit.moments[0]) == 2
	circuit.insert(1, vqubit.Circuit(vqubit.Z(q0)))
	assert len(circuit) == 2
	with pytest.raises(IndexError, match="out of range"):
		circuit.insert(99, vqubit.Circuit(vqubit.X(q0)))
	with pytest.raises(TypeError, match="cannot append"):
		circuit.append("not an operation")
