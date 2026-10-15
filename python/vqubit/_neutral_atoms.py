"""Helpers for neutral-atom devices."""

from __future__ import annotations

from typing import Any

from ._gates import CZGate
from ._gates_more import CZPowGate

_NATIVE_NEUTRAL_ATOM_GATES = (CZGate, CZPowGate)


def is_native_neutral_atom_gate(gate: Any) -> bool:
	"""Report whether ``gate`` is native to a neutral-atom device."""
	return isinstance(gate, _NATIVE_NEUTRAL_ATOM_GATES)


def is_native_neutral_atom_op(operation: Any) -> bool:
	"""Report whether ``operation``'s gate is native to a neutral-atom device."""
	return is_native_neutral_atom_gate(getattr(operation, "gate", operation))


__all__ = ["is_native_neutral_atom_gate", "is_native_neutral_atom_op"]
