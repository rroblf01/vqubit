"""The ``vqubit.experiments`` namespace, mirroring ``cirq.experiments``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.experiments`` so a
deep import (``from cirq.experiments import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	TensoredConfusionMatrices,
	estimate_parallel_single_qubit_readout_errors,
	estimate_single_qubit_readout_errors,
	hog_score_xeb_fidelity_from_probabilities,
	linear_xeb_fidelity,
	linear_xeb_fidelity_from_probabilities,
	log_xeb_fidelity,
	log_xeb_fidelity_from_probabilities,
	measure_confusion_matrix,
	xeb_fidelity,
)

__all__ = [
	"TensoredConfusionMatrices",
	"estimate_parallel_single_qubit_readout_errors",
	"estimate_single_qubit_readout_errors",
	"hog_score_xeb_fidelity_from_probabilities",
	"linear_xeb_fidelity",
	"linear_xeb_fidelity_from_probabilities",
	"log_xeb_fidelity",
	"log_xeb_fidelity_from_probabilities",
	"measure_confusion_matrix",
	"xeb_fidelity",
]
