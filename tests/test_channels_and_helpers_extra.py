"""Tests for the noisy channels, tags and value helpers.

The channels are checked the way a user checks a channel: its Kraus operators
are trace preserving, and applying it to a known state gives the state the
physics says it should.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import vqubit
from tests.helpers import assert_matrix_close

q0, q1, q2 = vqubit.LineQubit.range(3)


def density_after(channel: vqubit.Gate, state: object) -> np.ndarray:
	"""Apply a one-qubit channel to ``state`` and return the density matrix."""
	result = vqubit.DensityMatrixSimulator().simulate(vqubit.Circuit(channel.on(q0)), initial_state=state)
	return np.asarray(result.final_density_matrix.density_matrix())


def kraus_sum(channel: vqubit.Gate) -> np.ndarray:
	operators = channel._kraus_()
	assert operators is not None
	total = np.zeros((len(operators[0]), len(operators[0])), dtype=complex)
	for operator in operators:
		array = np.asarray(operator)
		total = total + array.conj().T @ array
	return total


# --- channels ----------------------------------------------------------------


def test_asymmetric_depolarizing_channel() -> None:
	channel = vqubit.AsymmetricDepolarizingChannel(p_x=0.1, p_y=0.2, p_z=0.3)
	assert channel.num_qubits == 1
	assert_matrix_close(kraus_sum(channel), np.eye(2))
	quiet = vqubit.AsymmetricDepolarizingChannel(p_x=0, p_y=0, p_z=0)
	quiet_ops = quiet._kraus_()
	assert quiet_ops is not None
	assert_matrix_close(np.asarray(quiet_ops[0]), np.eye(2))
	with pytest.raises(ValueError):
		vqubit.AsymmetricDepolarizingChannel(p_x=0.9, p_y=0.9, p_z=0)
	assert "AsymmetricDepolarizingChannel" in repr(channel)
	assert vqubit.asymmetric_depolarize(p_x=0.1) == vqubit.AsymmetricDepolarizingChannel(p_x=0.1)


def test_generalized_amplitude_damping() -> None:
	channel = vqubit.generalized_amplitude_damp(p=0.3, gamma=0.5)
	assert_matrix_close(kraus_sum(channel), np.eye(2))
	# Without damping the first Kraus operator is the identity.
	quiet = vqubit.generalized_amplitude_damp(p=1.0, gamma=0.0)
	quiet_ops = quiet._kraus_()
	assert quiet_ops is not None
	assert_matrix_close(np.asarray(quiet_ops[0]), np.eye(2))
	# Fully damped with the equilibrium in |0> sends |1> to |0>, as Cirq does.
	assert_matrix_close(density_after(vqubit.generalized_amplitude_damp(p=1.0, gamma=1.0), [0, 1]), [[1, 0], [0, 0]])
	# Half way there, the populations are even.
	assert_matrix_close(
		density_after(vqubit.generalized_amplitude_damp(p=0.5, gamma=1.0), [0, 1]), [[0.5, 0], [0, 0.5]]
	)


def test_kraus_channel_from_operators() -> None:
	# A bit flip with probability one is a single-operator Kraus channel.
	channel = vqubit.KrausChannel([vqubit.unitary(vqubit.X)], validate=True)
	assert_matrix_close(density_after(channel, [1, 0]), [[0, 0], [0, 1]])
	# The projection onto |0> alone is not trace preserving, and validation
	# says so instead of letting it through.
	with pytest.raises(ValueError, match="CPTP"):
		vqubit.KrausChannel([[[1, 0], [0, 0]]], validate=True)
	with pytest.raises(ValueError, match="square operator"):
		vqubit.KrausChannel([[[1, 0, 0], [0, 1, 0], [0, 0, 1]]])
	assert "KrausChannel" in repr(channel)


def test_mixed_unitary_channel() -> None:
	channel = vqubit.MixedUnitaryChannel([(0.5, vqubit.unitary(vqubit.I)), (0.5, vqubit.unitary(vqubit.X))])
	assert_matrix_close(kraus_sum(channel), np.eye(2))
	# A maximally mixed state stays maximally mixed: prepare it with a channel
	# and mix it again.
	# ``depolarize(0.75)`` is the fully depolarizing channel: any input becomes
	# the maximally mixed state.
	prepared = vqubit.depolarize(0.75)
	result = vqubit.DensityMatrixSimulator().simulate(
		vqubit.Circuit(prepared.on(q0), channel.on(q0)), initial_state=[1, 0]
	)
	assert_matrix_close(np.asarray(result.final_density_matrix.density_matrix()), [[0.5, 0], [0, 0.5]])
	# |0> becomes an even mixture of |0> and |1>.
	assert_matrix_close(density_after(channel, [1, 0]), [[0.5, 0], [0, 0.5]])
	with pytest.raises(ValueError, match="sum to 1"):
		vqubit.MixedUnitaryChannel([(0.5, vqubit.unitary(vqubit.I))])


def test_random_gate_channel() -> None:
	# With probability one the gate is always applied, with zero never.
	always = vqubit.RandomGateChannel(sub_gate=vqubit.X, probability=1.0)
	assert_matrix_close(density_after(always, [1, 0]), [[0, 0], [0, 1]])
	never = vqubit.RandomGateChannel(sub_gate=vqubit.X, probability=0.0)
	assert_matrix_close(density_after(never, [1, 0]), [[1, 0], [0, 0]])
	# Half and half mixes the two outcomes.
	half = vqubit.RandomGateChannel(sub_gate=vqubit.X, probability=0.5)
	assert_matrix_close(density_after(half, [1, 0]), [[0.5, 0], [0, 0.5]])
	assert "prob=0.5" in str(half)
	with pytest.raises(ValueError, match="not 0 <= probability"):
		vqubit.RandomGateChannel(sub_gate=vqubit.X, probability=1.5)


def test_state_preparation_channel() -> None:
	channel = vqubit.StatePreparationChannel([1, 0])
	# Whatever goes in comes out as |0>.
	assert_matrix_close(density_after(channel, [0, 1]), [[1, 0], [0, 0]])
	# A superposition target prepares that superposition.
	plus = vqubit.StatePreparationChannel([1 / math.sqrt(2), 1 / math.sqrt(2)])
	result = vqubit.DensityMatrixSimulator().simulate(vqubit.Circuit(plus.on(q0)), initial_state=[1, 0])
	assert_matrix_close(np.asarray(result.final_density_matrix.density_matrix()), [[0.5, 0.5], [0.5, 0.5]])
	with pytest.raises((ValueError, TypeError)):
		vqubit.StatePreparationChannel(np.eye(3))  # ty: ignore[invalid-argument-type]


def test_channel_validation() -> None:
	with pytest.raises(ValueError):
		vqubit.generalized_amplitude_damp(p=1.5, gamma=0.5)
	with pytest.raises(ValueError):
		vqubit.generalized_amplitude_damp(p=0.5, gamma=1.5)
	with pytest.raises(ValueError):
		vqubit.AsymmetricDepolarizingChannel(p_x=-0.1)


# --- tags --------------------------------------------------------------------


def test_tags_have_string_forms() -> None:
	virtual = vqubit.VirtualTag()
	routing = vqubit.RoutingSwapTag()
	assert str(virtual) == "virtual"
	assert repr(virtual) == "vqubit.VirtualTag()"
	assert str(routing) == "routing-swap"
	assert repr(routing) == "vqubit.RoutingSwapTag()"
	assert virtual == vqubit.VirtualTag()
	assert virtual != routing
	# A tagged operation carries its tags and still does the same thing.
	operation = vqubit.TaggedOperation(vqubit.X(q0), virtual)
	assert operation.tags == (virtual,)
	assert operation.untagged() == vqubit.X(q0)
	assert_matrix_close(vqubit.unitary(operation), vqubit.unitary(vqubit.X))
	with pytest.raises(ValueError, match="at least one tag"):
		vqubit.TaggedOperation(vqubit.X(q0))


# --- control values ----------------------------------------------------------


def test_control_value_formulas() -> None:
	first, second = vqubit.LineQubit.range(2)
	sum_of_products = vqubit.SumOfProducts([first, second], [[(first, True)], [(second, True)]])
	# "first OR second": anything except both false.
	assert (True, False) in sum_of_products
	assert (False, True) in sum_of_products
	assert (True, True) in sum_of_products
	assert (False, False) not in sum_of_products
	product_of_sums = vqubit.ProductOfSums([first, second], [[(first, True)], [(second, True)]])
	assert (True, True) in product_of_sums
	assert (True, False) not in product_of_sums
	assert sum_of_products == vqubit.SumOfProducts([first, second], [[(first, True)], [(second, True)]])
	assert "|" in sum_of_products.to_str()
	assert "&" in product_of_sums.to_str()
	assert "SumOfProducts" in repr(sum_of_products)
	with pytest.raises(ValueError, match="conditions length"):
		_ = (True,) in sum_of_products


def test_a_controlled_gate_with_a_formula_fires_on_every_assignment() -> None:
	control0, control1, target = vqubit.LineQubit.range(3)
	formula = vqubit.SumOfProducts(
		[control0, control1], [[(control0, True), (control1, False)], [(control0, False), (control1, True)]]
	)
	operation = vqubit.ControlledOperation([control0, control1], vqubit.X(target), formula)
	unitary = np.asarray(vqubit.unitary(operation))
	# The four control assignments, each |control0 control1 target>.
	for control_bits in ((0, 0), (0, 1), (1, 0), (1, 1)):
		index = (control_bits[0] << 2) | (control_bits[1] << 1)
		flipped = index | 1
		should_fire = control_bits in ((0, 1), (1, 0))
		assert abs(unitary[flipped][index]) == pytest.approx(1.0 if should_fire else 0.0), control_bits


def test_a_formula_without_satisfying_assignment_is_rejected() -> None:
	control = vqubit.LineQubit(0)
	with pytest.raises(ValueError, match="no satisfying"):
		vqubit.ControlledGate(vqubit.X, num_controls=1, control_values=vqubit.SumOfProducts([control], []))


# --- value helpers -----------------------------------------------------------


def test_bits_and_digits() -> None:
	assert vqubit.big_endian_bits_to_int([1, 0, 1]) == 5
	assert vqubit.big_endian_int_to_bits(5, bit_count=3) == [1, 0, 1]
	assert vqubit.big_endian_digits_to_int([1, 2], base=3) == 5
	assert vqubit.big_endian_int_to_digits(5, base=3) == [1, 2]
	assert vqubit.big_endian_int_to_digits(6, base=3) == [2, 0]


def test_angle_helpers() -> None:
	assert vqubit.canonicalize_half_turns(0.5) == pytest.approx(0.5)
	assert vqubit.canonicalize_half_turns(2.5) == pytest.approx(0.5)
	assert vqubit.chosen_angle_to_half_turns(math.pi) == pytest.approx(1.0)
	assert vqubit.chosen_angle_to_canonical_half_turns(3 * math.pi) == pytest.approx(1.0)


def test_probability_helpers() -> None:
	assert vqubit.state_vector_to_probabilities([1, 0]) == [1.0, 0.0]
	assert vqubit.validate_probability(0.5) == 0.5
	with pytest.raises(ValueError):
		vqubit.validate_probability(1.5)
	with pytest.raises(ValueError):
		vqubit.validate_probability(-0.1)
	assert vqubit.alternative(None, None, 3) == 3
	assert vqubit.alternative(None, default="d") == "d"


def test_linear_dict_and_duration() -> None:
	linear = vqubit.LinearDict({"a": 1})
	assert linear["a"] == 1
	linear["a"] += 2
	assert linear["a"] == 3
	assert dict(linear) == {"a": 3}
	assert len(linear) == 1
	duration = vqubit.Duration(picos=1000)
	assert duration.total_picos() == 1000
	assert vqubit.Timestamp(picos=5).raw_picos() == 5
