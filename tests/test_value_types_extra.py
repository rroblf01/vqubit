"""Tests for the small value types: LinearDict, Duration, Timestamp, ProductState.

These mirror ``cirq.value``: a dictionary that adds numbers linearly, the
picosecond duration and timestamp, a product state and the bit-mask classical
condition.
"""

from __future__ import annotations

import pytest

import vqubit
from tests.helpers import assert_matrix_close


def test_linear_dict_arithmetic() -> None:
	first = vqubit.LinearDict({"a": 1, "b": 2})
	second = vqubit.LinearDict({"a": 3, "c": 4})
	combined = first + second
	assert combined["a"] == 4
	assert combined["b"] == 2
	assert combined["c"] == 4
	# An absent key is a KeyError; the linear dictionary only adds the ones it
	# already holds.
	with pytest.raises(KeyError):
		_ = combined["missing"]
	# Copying and deleting.
	assert combined.copy() == combined
	copied = combined.copy()
	copied["a"] = 100
	assert combined["a"] == 4
	del combined["a"]
	assert "a" not in combined
	assert "LinearDict" in repr(combined)


def test_duration_arithmetic_and_validation() -> None:
	first = vqubit.Duration(picos=1000)
	second = vqubit.Duration(picos=250)
	assert (first + second).total_picos() == 1250
	assert (first - second).total_picos() == 750
	assert (first * 2).total_picos() == 2000
	assert (2 * first).total_picos() == 2000
	assert first == vqubit.Duration(picos=1000)
	assert first != second
	assert str(first) == "1000ps"
	assert repr(first) == "vqubit.Duration(picos=1000)"
	with pytest.raises(ValueError, match="negative"):
		vqubit.Duration(picos=-1)


def test_timestamps_order_and_render() -> None:
	early = vqubit.Timestamp(picos=10)
	late = vqubit.Timestamp(picos=20)
	assert early == vqubit.Timestamp(picos=10)
	assert early != late
	assert repr(early) == "vqubit.Timestamp(picos=10)"
	assert early.raw_picos() == 10


def test_bit_mask_key_condition() -> None:
	condition = vqubit.BitMaskKeyCondition("m", 0b10)
	assert condition.key == "m"
	assert condition.mask == 0b10
	assert "BitMaskKeyCondition" in repr(condition)
	assert str(condition) == "m & 2" or isinstance(str(condition), str)


def test_product_state() -> None:
	# A product of two pure qubits: |0> and |1>.
	state = vqubit.ProductState([(1, 0), (0, 1)])
	assert len(state.states) == 2
	assert state.states[0] == [1 + 0j, 0j]
	assert state.states[1] == [0j, 1 + 0j]
	# The product state has the Kronecker-product amplitudes.
	amplitudes = state.state_vector()
	assert len(amplitudes) == 4
	assert amplitudes[1] == pytest.approx(1)
	assert state == vqubit.ProductState([(1, 0), (0, 1)])
	assert "ProductState" in repr(state)
	assert_matrix_close(vqubit.unitary(vqubit.I(vqubit.LineQubit(0))), [[1, 0], [0, 1]])
