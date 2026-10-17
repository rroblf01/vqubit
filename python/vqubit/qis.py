"""The ``vqubit.qis`` namespace, mirroring ``cirq.qis``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.qis`` so a
deep import (``from cirq.qis import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	QUANTUM_STATE_LIKE,
	STATE_VECTOR_LIKE,
	CliffordTableau,
	QuantumState,
	QuantumStateRepresentation,
	StabilizerState,
	bloch_vector_from_state_vector,
	choi_to_kraus,
	choi_to_superoperator,
	density_matrix,
	density_matrix_from_state_vector,
	dirac_notation,
	entanglement_fidelity,
	eye_tensor,
	fidelity,
	kraus_to_choi,
	kraus_to_superoperator,
	one_hot,
	operation_to_choi,
	operation_to_superoperator,
	quantum_state,
	superoperator_to_choi,
	superoperator_to_kraus,
	to_valid_density_matrix,
	to_valid_state_vector,
	validate_density_matrix,
	validate_indices,
	validate_normalized_state_vector,
	validate_qid_shape,
	von_neumann_entropy,
)

__all__ = [
	"QUANTUM_STATE_LIKE",
	"STATE_VECTOR_LIKE",
	"CliffordTableau",
	"QuantumState",
	"QuantumStateRepresentation",
	"StabilizerState",
	"bloch_vector_from_state_vector",
	"choi_to_kraus",
	"choi_to_superoperator",
	"density_matrix",
	"density_matrix_from_state_vector",
	"dirac_notation",
	"entanglement_fidelity",
	"eye_tensor",
	"fidelity",
	"kraus_to_choi",
	"kraus_to_superoperator",
	"one_hot",
	"operation_to_choi",
	"operation_to_superoperator",
	"quantum_state",
	"superoperator_to_choi",
	"superoperator_to_kraus",
	"to_valid_density_matrix",
	"to_valid_state_vector",
	"validate_density_matrix",
	"validate_indices",
	"validate_normalized_state_vector",
	"validate_qid_shape",
	"von_neumann_entropy",
]
