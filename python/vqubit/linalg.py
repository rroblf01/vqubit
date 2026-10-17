"""The ``vqubit.linalg`` namespace, mirroring ``cirq.linalg``."""

from __future__ import annotations

from ._linalg_public import *  # noqa: F403
from ._linalg_public import __all__ as _public_all
from ._optimizers import deconstruct_single_qubit_matrix_into_angles
from ._qis import partial_trace
from ._vis import scatter_plot_normalized_kak_interaction_coefficients

__all__ = [
	*_public_all,
	"deconstruct_single_qubit_matrix_into_angles",
	"partial_trace",
	"scatter_plot_normalized_kak_interaction_coefficients",
]
