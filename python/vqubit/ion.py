"""The ``vqubit.ion`` namespace, mirroring ``cirq.ion``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.ion`` so a
deep import (``from cirq.ion import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	ms,
	two_qubit_matrix_to_ion_operations,
)

__all__ = [
	"ms",
	"two_qubit_matrix_to_ion_operations",
]
