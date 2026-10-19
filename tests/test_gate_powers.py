"""The plain gates *are* power gates, as in Cirq.

Cirq's ``CZ``, ``CNOT``, ``SWAP``, ``iSWAP`` and friends are power gates at
exponent 1 (``cirq.CZ`` is a ``CZPowGate``), so ``cirq.CZ ** 0.5`` gives a
square-root CZ, ``isinstance(cirq.CZ, cirq.CZPowGate)`` holds and
``cirq.rz(0.3) ** 0.5`` halves the angle. vqubit now does the same: the plain
gate names are aliases of the power-gate classes.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import pytest

import vqubit as vq
from tests.helpers import assert_matrix_close

# (radical, power-gate class, exponent that reproduces the radical)
POWER_CASES = [
	("CZ", "CZPowGate", 1.0),
	("CNOT", "CXPowGate", 1.0),
	("CX", "CXPowGate", 1.0),
	("SWAP", "SwapPowGate", 1.0),
	("ISWAP", "ISwapPowGate", 1.0),
	("SQRT_ISWAP", "ISwapPowGate", 0.5),
	("CCX", "CCXPowGate", 1.0),
	("TOFFOLI", "CCXPowGate", 1.0),
]


@pytest.mark.parametrize(("name", "power_class", "base"), POWER_CASES)
def test_a_plain_gate_is_an_instance_of_its_power_gate(name: str, power_class: str, base: float) -> None:
	gate = getattr(vq, name)
	assert isinstance(gate, getattr(vq, power_class))
	assert gate == getattr(vq, power_class)(exponent=base)


@pytest.mark.parametrize(("name", "power_class", "base"), POWER_CASES)
def test_a_plain_gate_raised_to_a_power_is_the_power_gate(name: str, power_class: str, base: float) -> None:
	gate = getattr(vq, name)
	expected = getattr(vq, power_class)(exponent=base * 0.5)
	assert gate**0.5 == expected
	assert isinstance(gate**0.5, getattr(vq, power_class))
	assert_matrix_close(vq.unitary(gate**0.5), vq.unitary(expected))


@pytest.mark.parametrize(("name", "power_class", "base"), POWER_CASES)
def test_the_plain_gate_is_that_power_gate_to_the_first(name: str, power_class: str, base: float) -> None:
	gate = getattr(vq, name)
	assert_matrix_close(vq.unitary(gate), vq.unitary(getattr(vq, power_class)(exponent=base)))


@pytest.mark.parametrize(("name", "power_class", "base"), POWER_CASES)
def test_power_one_round_trips(name: str, power_class: str, base: float) -> None:
	del base
	gate = getattr(vq, name)
	# A fresh object, but the same value -- exactly what ``cirq.CNOT ** 1`` does.
	assert gate**1 == gate
	assert gate**1 is not gate
	assert isinstance(gate**1, getattr(vq, power_class))


@pytest.mark.parametrize("name", ["CZ", "CNOT", "SWAP", "ISWAP", "SQRT_ISWAP", "CCX"])
def test_the_power_and_the_inverse_compose_to_the_identity(name: str) -> None:
	gate = getattr(vq, name)
	unitary = vq.unitary(gate)
	for quantity in (gate**-1, vq.inverse(gate)):
		product = _matmul(unitary, vq.unitary(quantity))
		assert_matrix_close(product, _identity(len(product)))


def test_identity_is_its_own_power() -> None:
	assert vq.I**0.5 == vq.I
	assert_matrix_close(vq.unitary(vq.I**0.5), vq.unitary(vq.I))


@pytest.mark.parametrize(("factory", "gate_class"), [(vq.rx, "RxGate"), (vq.ry, "RyGate"), (vq.rz, "RzGate")])
def test_rotations_scale_their_angle(factory: Callable[[float], vq.Gate], gate_class: str) -> None:
	gate = factory(0.3)
	assert gate**0.5 == getattr(vq, gate_class)(0.15)
	assert gate**2 == getattr(vq, gate_class)(0.6)
	assert_matrix_close(vq.unitary(gate**0.5), vq.unitary(getattr(vq, gate_class)(0.15)))


@pytest.mark.parametrize("name", ["FREDKIN", "CSWAP"])
def test_fredkin_is_not_powerable_like_in_cirq(name: str) -> None:
	# cirq.CSWAP ** 0.5 raises TypeError; vqubit says so explicitly.
	with pytest.raises(NotImplementedError):
		getattr(vq, name) ** 0.5


def test_a_symbolic_power_stays_symbolic_until_resolved() -> None:
	t: Any = vq.symbol("t")
	half = vq.rz(t) ** 0.5
	assert vq.is_parameterized(half)
	assert vq.resolve_parameters(half, vq.ParamResolver({t: math.pi / 2})) == vq.rz(math.pi / 4)
	scaled = vq.resolve_parameters(vq.CZ**t, vq.ParamResolver({t: 1}))
	assert_matrix_close(vq.unitary(scaled), vq.unitary(vq.CZ))


def test_iswap_inverse_is_iswap_inv_not_iswap() -> None:
	# Regression: _FixedGate used to return `type(self)()`, so the inverse of
	# iSWAP was iSWAP, and U @ U^-1 was not the identity.
	assert vq.inverse(vq.ISWAP) == vq.ISWAP_INV
	assert vq.inverse(vq.SQRT_ISWAP) == vq.SQRT_ISWAP_INV


def test_a_fractional_power_matches_the_analytic_matrix() -> None:
	# A square-root CZ squares back to CZ, which pins the branch cut.
	root = vq.unitary(vq.CZ**0.5)
	assert_matrix_close(_matmul(root, root), vq.unitary(vq.CZ))


def test_the_plain_gate_class_names_are_aliases() -> None:
	"""``vq.CZGate`` used to be a class of its own; it names the power gate now.

	The divergence this replaces is gone from the differences page, so keep the
	two in step the way that page says: this is the fold, and
	``vq.CZGate is vq.CZPowGate`` is what makes ``isinstance`` and ``==`` agree
	with Cirq.
	"""
	assert vq.CZGate is vq.CZPowGate
	assert vq.CNOTGate is vq.CXPowGate
	assert vq.SWAPGate is vq.SwapPowGate
	assert vq.CCXGate is vq.CCXPowGate
	assert vq.ISWAPGate is vq.ISwapPowGate
	# ``OpIdentifier`` takes a gate type, so the aliases have to stay usable.
	assert vq.OpIdentifier(vq.CNOTGate, 2).gate_type is vq.CXPowGate


def test_value_equality_uses_the_raw_exponent_not_the_cirq_period() -> None:
	"""The documented divergence: Cirq folds the exponent modulo the gate period.

	``cirq.CNOT ** -1 == cirq.CNOT`` and ``cirq.X ** 3 == cirq.X`` are both True,
	because Cirq compares the eigenphases modulo 2 (its value equality stores a
	``PeriodicValue``). vqubit compares the raw exponent, so both are False here.
	The git history and ``guide/cirq-differences`` record why; if this ever
	changes, update that page together with this test.
	"""
	assert vq.CNOT**-1 != vq.CNOT
	assert vq.X**3 != vq.X
	assert vq.ISWAP**-1 != vq.ISWAP
	assert_matrix_close(vq.unitary(vq.CNOT**-1), vq.unitary(vq.CNOT))
	assert_matrix_close(vq.unitary(vq.X**3), vq.unitary(vq.X))


def _matmul(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
	return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def _identity(dim: int) -> list[list[complex]]:
	return [[1.0 + 0j if i == j else 0j for j in range(dim)] for i in range(dim)]
