"""Helpers for neutral-atom devices."""

from __future__ import annotations

from typing import Any

from ._gates import (
	CCXPowGate,
	CCZPowGate,
	CXPowGate,
	CZPowGate,
	IdentityGate,
	MeasurementGate,
	PhasedXPowGate,
	XPowGate,
	YPowGate,
	ZPowGate,
)
from ._gateset import AnyIntegerPowerGateFamily, Gateset, ParallelGateFamily

# The gateset ``cirq.neutral_atom_devices.neutral_atom_gateset`` builds: any
# integer power of ``CNOT``/``CCX``/``CZ``/``CCZ``, the ``X``/``Y``/``Z``/
# ``PhasedX`` families in parallel, plus measurement and identity. Note what it
# leaves out, because it is easy to get wrong: ``H`` is not native, and a
# *fractional* power of ``CZ`` (``vq.CZ ** 0.5``) is native in neither library.
_NATIVE_NEUTRAL_ATOM_GATESET = Gateset(
	[
		AnyIntegerPowerGateFamily(CXPowGate),
		AnyIntegerPowerGateFamily(CCXPowGate),
		AnyIntegerPowerGateFamily(CZPowGate),
		AnyIntegerPowerGateFamily(CCZPowGate),
		ParallelGateFamily(ZPowGate),
		ParallelGateFamily(XPowGate),
		ParallelGateFamily(YPowGate),
		ParallelGateFamily(PhasedXPowGate),
		MeasurementGate,
		IdentityGate,
	],
	name="neutral_atom_gateset",
)


def is_native_neutral_atom_gate(gate: Any) -> bool:
	"""Report whether ``gate`` is in the default neutral-atom gateset."""
	return gate in _NATIVE_NEUTRAL_ATOM_GATESET


def is_native_neutral_atom_op(operation: Any) -> bool:
	"""Report whether ``operation`` is in the default neutral-atom gateset."""
	return operation in _NATIVE_NEUTRAL_ATOM_GATESET


__all__ = ["is_native_neutral_atom_gate", "is_native_neutral_atom_op"]
