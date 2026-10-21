"""Tests for the parameter expressions, the base device and the channel edges.

The expression tree has operations that only show up when you use them
(subtraction, unary minus, ``names``); the device base class and the Sycamore
gate have a few branches of their own; and the extra channels validate their
arguments. These tests exercise those paths with hand-checkable values.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

import vqubit
from tests.helpers import assert_matrix_close

# --- parameter expressions ---------------------------------------------------


def test_expression_names_and_negation() -> None:
	t = vqubit.Parameter("t")
	s = vqubit.Parameter("s")
	negated = -t
	assert negated.names() == {"t"}
	assert str(negated) == "(-t)"
	assert repr(negated) == "(-vqubit.Parameter('t'))"
	# Substitution resolves a negation to the negated number.
	assert negated.substitute(vqubit.ParamResolver({"t": 0.5})) == -0.5
	# Names come from both sides of a binary expression.
	assert (t + s).names() == {"t", "s"}
	assert (t - s).names() == {"t", "s"}
	assert (t * s / 2).names() == {"t", "s"}
	# Subtraction evaluates, and an unbound parameter stays as it is.
	assert (t - s).substitute(vqubit.ParamResolver({"t": 1.0, "s": 0.25})) == 0.75
	unbound = (t - s).substitute(vqubit.ParamResolver({"t": 1.0}))
	assert vqubit.is_parameterized(unbound)


def test_expression_evaluation_covers_every_operator() -> None:
	t = vqubit.Parameter("t")
	resolver = vqubit.ParamResolver({"t": 3.0})
	assert (t + 1).substitute(resolver) == 4.0
	assert (t - 1).substitute(resolver) == 2.0
	assert (1 - t).substitute(resolver) == -2.0
	assert (t * 2).substitute(resolver) == 6.0
	assert (2 * t).substitute(resolver) == 6.0
	assert (t / 3).substitute(resolver) == 1.0
	assert (3 / t).substitute(resolver) == 1.0
	assert (t**2).substitute(resolver) == 9.0
	assert (-t).substitute(resolver) == -3.0


def test_param_resolver_forms() -> None:
	resolver = vqubit.ParamResolver({"t": 0.5})
	assert resolver.value_of(vqubit.Parameter("t")) == 0.5
	assert resolver["t"] == 0.5
	assert list(resolver) == ["t"]
	assert len(resolver) == 1
	assert repr(resolver) == "vqubit.ParamResolver({'t': 0.5})"
	# An unbound name stays symbolic, and the kwargs form builds the map.
	assert vqubit.is_parameterized(resolver.value_of(vqubit.Parameter("s")))
	from_kwargs = vqubit.ParamResolver(t=0.25)
	assert from_kwargs.value_of(vqubit.Parameter("t")) == 0.25
	# A resolver is one point of a sweep.
	assert vqubit.to_resolvers(resolver) == [resolver]


# --- devices -----------------------------------------------------------------


def test_device_metadata_base_class() -> None:
	metadata = vqubit.DeviceMetadata(name="toy")
	assert metadata.name == "toy"
	assert repr(metadata) == "vqubit.DeviceMetadata(name='toy')"
	# The base class has no qubits of its own and no durations for them.
	assert metadata.qubit_set is None or tuple(metadata.qubit_set) == ()
	assert metadata.duration_of(vqubit.X(vqubit.LineQubit(0))) is None


def test_sycamore_gate_and_device() -> None:
	assert vqubit.SycamoreGate().num_qubits == 2
	assert vqubit.SycamoreGate()._diagram_symbols_([vqubit.LineQubit(0), vqubit.LineQubit(1)]) == ("Syc", "Syc")
	assert vqubit.SycamoreGate()._qid_shape_() == (2, 2)
	assert str(vqubit.SycamoreGate()) == "SycamoreGate"
	assert vqubit.has_unitary(vqubit.SycamoreGate())
	# The device is a grid with a metadata object of its own.
	device = vqubit.Sycamore()
	assert len(device.qubits) > 0
	assert device.metadata.name == "Sycamore"
	assert callable(device.validate_operation)


def test_adjacency_validation_rejects_a_non_adjacent_pair() -> None:
	# Two qubits the grid does not join cannot take a two-qubit gate.
	qubits = list(vqubit.GridQubit.square(2))
	metadata = vqubit.GridDeviceMetadata(qubits, [(qubits[0], qubits[1])])
	device = vqubit.GridDevice(qubits, metadata=metadata)
	with pytest.raises(ValueError, match="not adjacent"):
		device.validate_operation(vqubit.CNOT(qubits[0], qubits[3]))
	# The adjacent pair is accepted.
	device.validate_operation(vqubit.CNOT(qubits[0], qubits[1]))


# --- channels ----------------------------------------------------------------


def test_channel_edge_cases() -> None:
	# Amplitude damping with probability zero and one.
	quiet = vqubit.amplitude_damp(0.0)
	assert_matrix_close(np_kraus_sum(quiet), [[1, 0], [0, 1]])
	full = vqubit.amplitude_damp(1.0)
	result = vqubit.DensityMatrixSimulator().simulate(
		vqubit.Circuit(full.on(vqubit.LineQubit(0))), initial_state=[0, 1]
	)
	assert_matrix_close(np_asarray(result.final_density_matrix), [[1, 0], [0, 0]])
	with pytest.raises(ValueError):
		vqubit.amplitude_damp(1.5)


def test_phase_and_bit_flip_channels() -> None:
	for channel, reference in (
		(vqubit.bit_flip(1.0), [[0, 0], [0, 1]]),
		(vqubit.phase_flip(1.0), [[1, 0], [0, 0]]),
		(vqubit.bit_phase_flip(1.0), [[0, 0], [0, 1]]),
	):
		result = vqubit.DensityMatrixSimulator().simulate(
			vqubit.Circuit(channel.on(vqubit.LineQubit(0))), initial_state=[1, 0]
		)
		assert_matrix_close(np_asarray(result.final_density_matrix), reference)
	with pytest.raises(ValueError):
		vqubit.bit_flip(2.0)


def np_asarray(value: object) -> object:
	return value


def np_kraus_sum(channel: vqubit.Gate) -> Any:
	np = pytest.importorskip("numpy", exc_type=ImportError)

	operators = channel._kraus_()
	assert operators is not None
	total = np.zeros((len(operators[0]), len(operators[0])), dtype=complex)
	for operator in operators:
		array = np.asarray(operator)
		total = total + array.conj().T @ array
	return total


def test_rotation_angle_helpers() -> None:
	assert vqubit.canonicalize_half_turns(0.25) == pytest.approx(0.25)
	assert vqubit.canonicalize_half_turns(1.25) == pytest.approx(-0.75)
	assert vqubit.chosen_angle_to_half_turns(math.pi / 2) == pytest.approx(0.5)
	assert vqubit.state_vector_to_probabilities([1 / math.sqrt(2), 1j / math.sqrt(2)]) == pytest.approx([0.5, 0.5])
