"""Helpers for neutral-atom devices."""

from __future__ import annotations

from typing import Any

from ._gates import CZPowGate

# ``CZGate`` and ``CZPowGate`` used to be listed together; they are the same
# class now (``vq.CZ`` *is* a ``CZPowGate``), so one entry says the same thing.
#
# Cirq's neutral-atom gateset is wider than this: it also takes integer powers
# of ``CZ``, ``CNOT`` and ``CCZ``, the ``X``/``Y``/``Z``/``PhasedX`` families and
# the measurement and identity gates, and rejects ``CZ ** 0.5``. That gap is
# recorded in ROADMAP.md rather than papered over here.
_NATIVE_NEUTRAL_ATOM_GATES = (CZPowGate,)


def is_native_neutral_atom_gate(gate: Any) -> bool:
	"""Report whether ``gate`` is native to a neutral-atom device."""
	return isinstance(gate, _NATIVE_NEUTRAL_ATOM_GATES)


def is_native_neutral_atom_op(operation: Any) -> bool:
	"""Report whether ``operation``'s gate is native to a neutral-atom device."""
	return is_native_neutral_atom_gate(getattr(operation, "gate", operation))


__all__ = ["is_native_neutral_atom_gate", "is_native_neutral_atom_op"]
