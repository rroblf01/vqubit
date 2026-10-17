"""The ``vqubit.neutral_atoms`` namespace, mirroring ``cirq.neutral_atoms``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.neutral_atoms`` so a
deep import (``from cirq.neutral_atoms import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	is_native_neutral_atom_gate,
	is_native_neutral_atom_op,
)

__all__ = [
	"is_native_neutral_atom_gate",
	"is_native_neutral_atom_op",
]
