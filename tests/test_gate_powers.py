"""The plain gates answer ``**`` the way their Cirq counterparts do.

Cirq's ``CZ``, ``CNOT``, ``SWAP``, ``iSWAP`` and friends *are* power gates
(``cirq.CZ is cirq.CZPowGate(exponent=1)``), so ``cirq.CZ ** 0.5`` gives a
square-root CZ and ``cirq.rz(0.3) ** 0.5`` halves the angle. vqubit keeps the
fixed gates as their own classes, so ``**`` has to bridge to the power gate.
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
def test_a_plain_gate_raised_to_a_power_is_the_power_gate(name: str, power_class: str, base: float) -> None:
	gate = getattr(vq, name)
	expected = getattr(vq, power_class)(exponent=base * 0.5)
	assert gate**0.5 == expected
	assert_matrix_close(vq.unitary(gate**0.5), vq.unitary(expected))


@pytest.mark.parametrize(("name", "power_class", "base"), POWER_CASES)
def test_the_plain_gate_is_that_power_gate_to_the_first(name: str, power_class: str, base: float) -> None:
	gate = getattr(vq, name)
	assert_matrix_close(vq.unitary(gate), vq.unitary(getattr(vq, power_class)(exponent=base)))


@pytest.mark.parametrize(("name", "power_class", "base"), POWER_CASES)
def test_power_one_round_trips(name: str, power_class: str, base: float) -> None:
	del power_class, base
	assert getattr(vq, name) ** 1 is getattr(vq, name)


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


def test_the_documented_class_divergence_still_holds() -> None:
	"""``**`` bridges to the power gate; the plain gate is still its own class.

	This is the divergence documented in ``docs/*/guide/cirq-differences.md``
	("No `**` on operations"). If ``vq.CZ`` ever *becomes* a ``CZPowGate``
	(roadmap 0.3), delete this test together with that paragraph instead of
	letting the two drift apart.
	"""
	assert not isinstance(vq.CZ, vq.CZPowGate)
	assert vq.CZPowGate() != vq.CZ


def _matmul(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
	return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def _identity(dim: int) -> list[list[complex]]:
	return [[1.0 + 0j if i == j else 0j for j in range(dim)] for i in range(dim)]
