"""Tests for the extra cirq.value helpers."""

from __future__ import annotations

import math

import pytest

import vqubit


def test_big_endian_helpers() -> None:
	assert vqubit.big_endian_bits_to_int([1, 0, 1]) == 5
	assert vqubit.big_endian_int_to_bits(5, bit_count=3) == [1, 0, 1]
	assert vqubit.big_endian_digits_to_int([1, 2], base=3) == 5
	assert vqubit.big_endian_int_to_digits(5, base=3) == [1, 2]


def test_angle_helpers() -> None:
	assert vqubit.canonicalize_half_turns(3.5) == pytest.approx(-0.5)
	assert vqubit.chosen_angle_to_half_turns(math.pi) == pytest.approx(1)
	assert vqubit.chosen_angle_to_canonical_half_turns(3 * math.pi) == pytest.approx(1)


def test_validate_probability() -> None:
	assert vqubit.validate_probability(0.5) == 0.5
	with pytest.raises(ValueError):
		vqubit.validate_probability(1.5)


def test_alternative() -> None:
	assert vqubit.alternative(None, 3, 4) == 3
	assert vqubit.alternative(None, None, default="x") == "x"


def test_linear_dict() -> None:
	mapping = vqubit.LinearDict({"a": 1, "b": 0})
	assert dict(mapping) == {"a": 1}
	mapping["a"] = 0
	assert dict(mapping) == {}
	mapping["c"] = 2
	assert dict(mapping) == {"c": 2}


def test_duration_and_timestamp() -> None:
	duration = vqubit.Duration(picos=1000)
	assert duration.total_nanos() == 1
	assert (duration + duration).picos == 2000
	timestamp = vqubit.Timestamp(picos=5000)
	assert (timestamp + duration).raw_picos() == 6000
	assert (timestamp - vqubit.Timestamp(picos=4000)).picos == 1000


def test_periodic_value() -> None:
	# The value is normalized on the way in, as Cirq's is, so both read ``1`` and
	# hash alike.
	assert vqubit.PeriodicValue(3.0, 2.0) == vqubit.PeriodicValue(1.0, 2.0)
	assert vqubit.PeriodicValue(3.0, 2.0).value == 1.0
	assert hash(vqubit.PeriodicValue(3.0, 2.0)) == hash(vqubit.PeriodicValue(1.0, 2.0))
	assert vqubit.PeriodicValue(-0.5, 2.0) == vqubit.PeriodicValue(1.5, 2.0)
	assert vqubit.PeriodicValue(1.0, 2.0) != vqubit.PeriodicValue(1.0, 4.0)


def test_measurement_type() -> None:
	assert vqubit.MeasurementType.Z.value == "Z"
	assert vqubit.MeasurementType.P.value == "P"


def test_ket_constants_and_product_state() -> None:
	assert vqubit.KET_ZERO == [[1 + 0j], [0j]]
	state = vqubit.ProductState([vqubit.KET_PLUS, vqubit.KET_ZERO])
	vector = state.state_vector()
	assert len(vector) == 4
	assert vector[0] == pytest.approx(1 / math.sqrt(2))
	assert vector[1] == pytest.approx(0)


def test_classical_data_store() -> None:
	store = vqubit.ClassicalDataDictionaryStore(measurements={"m": [0, 1]})
	store.record(m=1)
	assert list(store.keys()) == ["m"]
	assert store.records("m") == [1]
	assert store.measurements("m") == [0, 1]


def test_conditions() -> None:
	assert isinstance(vqubit.BitMaskKeyCondition("m", 3), vqubit.BitMaskKeyCondition)
	assert isinstance(vqubit.SympyCondition("x & y"), vqubit.SympyCondition)
	assert vqubit.MEASUREMENT_KEY_SEPARATOR == "."
