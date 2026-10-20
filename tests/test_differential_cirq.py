"""Differential tests against an installed Cirq.

The harness in ``benchmarks/parity.py`` compares the API shape; these tests
compare *answers* on a deterministic grid, so the equality rules and the
matrices are checked against cirq-core itself rather than against hand-written
expectations. They skip when Cirq is not installed (it is optional, and only the
dedicated CI job installs it).
"""

from __future__ import annotations

import itertools

import pytest

import vqubit

cirq = pytest.importorskip("cirq")

PHASED_XZ_VALUES = (0.0, 0.25, 0.5, 1.0, 1.5, 2.5)


def test_phased_xz_equality_matches_cirq_on_a_grid() -> None:
	# Every way of writing the same gate (plus two turns on any exponent) has to
	# compare as Cirq compares it: the canonical exponents, modulo their period.
	pairs = 0
	for x, z, a in itertools.product((0.0, 0.5, 1.0, 1.5), (0.0, 0.5, 1.5), (0.25, 0.5, 1.25)):
		first = vqubit.PhasedXZGate(x_exponent=x, z_exponent=z, axis_phase_exponent=a)
		cirq_first = cirq.PhasedXZGate(x_exponent=x, z_exponent=z, axis_phase_exponent=a)
		for x2, z2, a2 in itertools.product((x, x + 2, x - 2), (z, z + 2), (a, a + 2)):
			second = vqubit.PhasedXZGate(x_exponent=x2, z_exponent=z2, axis_phase_exponent=a2)
			cirq_second = cirq.PhasedXZGate(x_exponent=x2, z_exponent=z2, axis_phase_exponent=a2)
			assert (first == second) == (cirq_first == cirq_second), (x, z, a, x2, z2, a2)
			assert hash(first) == hash(second) if first == second else True
			pairs += 1
	assert pairs > 100


def test_phased_xz_canonical_forms_agree() -> None:
	# The canonical form is what equality sees, so the matrices must agree too.
	for x, z, a in itertools.product(PHASED_XZ_VALUES, (0.0, 0.5, 1.5), (0.25, 0.75)):
		gate = vqubit.PhasedXZGate(x_exponent=x, z_exponent=z, axis_phase_exponent=a)
		canonical = gate._canonical()
		cirq_gate = cirq.PhasedXZGate(x_exponent=x, z_exponent=z, axis_phase_exponent=a)
		cirq_canonical = cirq_gate._canonical()
		assert (canonical.x_exponent, canonical.z_exponent, canonical.axis_phase_exponent) == pytest.approx(
			(cirq_canonical.x_exponent, cirq_canonical.z_exponent, cirq_canonical.axis_phase_exponent)
		)
		assert vqubit.equal_up_to_global_phase(
			gate,
			vqubit.PhasedXZGate(
				x_exponent=canonical.x_exponent,
				z_exponent=canonical.z_exponent,
				axis_phase_exponent=canonical.axis_phase_exponent,
			),
		)


CLIFFORD_POWERS = (-2, -1, 0, 1, 2, 3, 4, 5)


@pytest.mark.parametrize("name", ["CZSWAP", "CXSWAP", "SWAP", "ISWAP", "CNOT", "CZ"])
def test_clifford_power_equality_matches_cirq(name: str) -> None:
	# Two Clifford gates are equal when their tableaux are, which for the
	# matrices means equality up to a global phase (and without floating-point
	# noise from the multiplications).
	our_gate = getattr(vqubit, name)
	cirq_gate = getattr(cirq, name)
	for power in CLIFFORD_POWERS:
		ours = our_gate**power if power != 1 else our_gate
		theirs = cirq_gate**power if power != 1 else cirq_gate
		assert (ours == our_gate) == (theirs == cirq_gate), (name, power)


def test_clifford_unitaries_match_cirq() -> None:
	import numpy as np

	for name in ("CZSWAP", "CXSWAP", "SWAP", "ISWAP", "CNOT", "CZ", "SQRT_ISWAP", "CCX"):
		ours = np.asarray(vqubit.unitary(getattr(vqubit, name)))
		theirs = np.asarray(cirq.unitary(getattr(cirq, name)))
		assert cirq.allclose_up_to_global_phase(ours, theirs), name
