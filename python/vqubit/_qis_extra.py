"""Extra ``cirq.qis`` helpers: validation, Choi/superoperator conversions and
the quantum-state wrapper types.

These build on the state helpers in :mod:`vqubit._qis`. A Choi matrix is the
``d**2`` by ``d**2`` matrix ``J[(a,b),(c,d)] = sum_k K_k[a,b] conj(K_k[c,d])``;
the superoperator ``S`` with ``S[(a,c),(b,d)] = J[(a,b),(c,d)]`` acts on
column-stacked density matrices.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from ._linalg import Matrix
from ._linalg_public import _hermitian_eig
from ._qis import partial_trace, to_valid_density_matrix, to_valid_state_vector
from ._value import value_equality

STATE_VECTOR_LIKE = Any
QUANTUM_STATE_LIKE = Any


def _product(values: Sequence[int]) -> int:
	result = 1
	for value in values:
		result *= value
	return result


def validate_qid_shape(qid_shape: Sequence[int], num_qubits: int | None = None) -> None:
	"""Raise ``ValueError`` if ``qid_shape`` is not a valid qid shape."""
	if num_qubits is not None and len(qid_shape) != num_qubits:
		raise ValueError(f"len(qid_shape) != num_qubits: {len(qid_shape)} != {num_qubits}")
	for dimension in qid_shape:
		if not isinstance(dimension, int) or dimension < 1:
			raise ValueError(f"qid dimensions must be positive integers, got {dimension!r}")


def validate_normalized_state_vector(
	state_vector: Sequence[complex],
	qid_shape: Sequence[int] | None = None,
	dtype: Any = None,
	atol: float = 1e-7,
) -> None:
	"""Raise ``ValueError`` if ``state_vector`` is not a normalized state."""
	del dtype
	dim = _product(qid_shape) if qid_shape is not None else len(state_vector)
	if len(state_vector) != dim:
		raise ValueError(f"state_vector has length {len(state_vector)}, expected {dim}")
	norm = sum(abs(complex(z)) ** 2 for z in state_vector)
	if abs(norm - 1.0) > atol:
		raise ValueError(f"state_vector is not normalized (norm squared {norm})")


def validate_density_matrix(
	density_matrix: Sequence[Sequence[complex]],
	qid_shape: Sequence[int] | None = None,
	dtype: Any = None,
	atol: float = 1e-7,
) -> None:
	"""Raise ``ValueError`` if ``density_matrix`` is not a valid density matrix."""
	del dtype
	dim = _product(qid_shape) if qid_shape is not None else len(density_matrix)
	if len(density_matrix) != dim or any(len(row) != dim for row in density_matrix):
		raise ValueError(f"density_matrix must be square of dimension {dim}")
	trace = sum(complex(density_matrix[i][i]) for i in range(dim))
	if abs(trace - 1.0) > atol:
		raise ValueError(f"density_matrix trace is {trace}, expected 1")


def validate_indices(indices: Sequence[int], shape: Sequence[int]) -> None:
	"""Raise ``ValueError`` if ``indices`` are not valid for ``shape``."""
	if len(indices) != len(shape):
		raise ValueError(f"len(indices)={len(indices)} != len(shape)={len(shape)}")
	for index, dimension in zip(indices, shape, strict=True):
		if not 0 <= index < dimension:
			raise ValueError(f"index {index} out of range for dimension {dimension}")


def density_matrix_from_state_vector(state_vector: Sequence[complex], indices: Sequence[int] | None = None) -> Any:
	"""Return ``|psi><psi|`` (optionally traced down to ``indices``)."""
	vector = [complex(z) for z in state_vector]
	rho = [[vector[i] * vector[j].conjugate() for j in range(len(vector))] for i in range(len(vector))]
	if indices is None:
		return rho
	return partial_trace(rho, indices)


def dirac_notation(state: Sequence[complex], decimals: int = 2) -> str:
	"""Return a Dirac-notation string for a state vector."""
	vector = [complex(z) for z in state]
	num_qubits = len(vector).bit_length() - 1
	terms = []
	for index, amplitude in enumerate(vector):
		magnitude = abs(amplitude)
		if magnitude < 10 ** (-decimals - 1):
			continue
		phase = amplitude / magnitude
		label = format(index, f"0{num_qubits}b")
		if abs(phase - 1) < 1e-9:
			term = f"{magnitude:.{decimals}f}|{label}>"
		elif abs(phase + 1) < 1e-9:
			term = f"-{magnitude:.{decimals}f}|{label}>"
		elif abs(phase - 1j) < 1e-9:
			term = f"{magnitude:.{decimals}f}i|{label}>"
		elif abs(phase + 1j) < 1e-9:
			term = f"-{magnitude:.{decimals}f}i|{label}>"
		else:
			term = f"({amplitude:.{decimals}f})|{label}>"
		terms.append(term)
	return " + ".join(terms) if terms else "0"


def one_hot(index: Sequence[int], shape: Sequence[int], value: complex = 1, dtype: Any = None) -> Any:
	"""Return a tensor that is ``value`` at ``index`` and zero elsewhere."""
	del dtype
	validate_indices(index, shape)

	def build(dims: Sequence[int], prefix: Sequence[int]) -> Any:
		if not dims:
			return complex(value)
		position = len(prefix)
		return [
			build(dims[1:], (*prefix, i)) if i == index[position] else _zeros_tensor(dims[1:]) for i in range(dims[0])
		]

	return build(shape, ())


def _zeros_tensor(dims: Sequence[int]) -> Any:
	if not dims:
		return 0j
	return [_zeros_tensor(dims[1:]) for _ in range(dims[0])]


def eye_tensor(shape: Sequence[int], *, dtype: Any = None, split_index: int | None = None) -> Any:
	"""Return the identity tensor of the given shape (1 where all indices match)."""
	del dtype, split_index

	def build(dims: Sequence[int], prefix: Sequence[int]) -> Any:
		if not dims:
			return 1 + 0j if len(set(prefix)) <= 1 else 0j
		return [build(dims[1:], (*prefix, i)) for i in range(dims[0])]

	return build(shape, ())


def entanglement_fidelity(state1: Any, state2: Any, qid_shape: Sequence[int] | None = None) -> float:
	"""Return the entanglement fidelity between two states."""
	del qid_shape
	from ._qis import fidelity

	return fidelity(state1, state2) ** 2


def kraus_to_choi(kraus_ops: Sequence[Sequence[Sequence[complex]]]) -> Matrix:
	"""Return the Choi matrix of a set of Kraus operators."""
	operators = [[[complex(z) for z in row] for row in op] for op in kraus_ops]
	d = len(operators[0])
	choi = [[0j] * (d * d) for _ in range(d * d)]
	for op in operators:
		for a in range(d):
			for b in range(d):
				for c in range(d):
					for e in range(d):
						choi[a * d + b][c * d + e] += op[a][b] * op[c][e].conjugate()
	return choi


def choi_to_kraus(choi: Sequence[Sequence[complex]], atol: float = 1e-8) -> list[Matrix]:
	"""Return Kraus operators from a Choi matrix."""
	matrix = [[complex(z) for z in row] for row in choi]
	d = round(math.sqrt(len(matrix)))
	values, vectors = _hermitian_eig(matrix)
	result: list[Matrix] = []
	for value, vector in zip(values, vectors, strict=True):
		if value > atol:
			scale = math.sqrt(value)
			result.append([[scale * vector[i * d + j] for j in range(d)] for i in range(d)])
	return result


def kraus_to_superoperator(kraus_ops: Sequence[Sequence[Sequence[complex]]]) -> Matrix:
	"""Return the superoperator of a set of Kraus operators."""
	operators = [[[complex(z) for z in row] for row in op] for op in kraus_ops]
	d = len(operators[0])
	superoperator = [[0j] * (d * d) for _ in range(d * d)]
	for op in operators:
		for a in range(d):
			for b in range(d):
				for c in range(d):
					for e in range(d):
						superoperator[a * d + c][b * d + e] += op[a][b] * op[c][e].conjugate()
	return superoperator


def choi_to_superoperator(choi: Sequence[Sequence[complex]]) -> Matrix:
	"""Return the superoperator of a Choi matrix."""
	return _transpose_axes([[complex(z) for z in row] for row in choi])


def superoperator_to_choi(superoperator: Sequence[Sequence[complex]]) -> Matrix:
	"""Return the Choi matrix of a superoperator."""
	return _transpose_axes([[complex(z) for z in row] for row in superoperator])


def _transpose_axes(matrix: Matrix) -> Matrix:
	d = round(math.sqrt(len(matrix)))
	tensor = [[[[0j] * d for _ in range(d)] for _ in range(d)] for _ in range(d)]
	for r in range(d * d):
		for c in range(d * d):
			tensor[r // d][r % d][c // d][c % d] = matrix[r][c]
	# Transpose (0,2,1,3).
	result = [[0j] * (d * d) for _ in range(d * d)]
	for a in range(d):
		for b in range(d):
			for c in range(d):
				for e in range(d):
					result[a * d + c][b * d + e] = tensor[a][b][c][e]
	return result


def superoperator_to_kraus(superoperator: Sequence[Sequence[complex]], atol: float = 1e-8) -> list[Matrix]:
	"""Return Kraus operators from a superoperator."""
	return choi_to_kraus(superoperator_to_choi(superoperator), atol=atol)


def operation_to_choi(operation: Any) -> Matrix:
	"""Return the Choi matrix of an operation."""
	from ._protocols import kraus

	return kraus_to_choi(kraus(operation))


def operation_to_superoperator(operation: Any) -> Matrix:
	"""Return the superoperator of an operation."""
	from ._protocols import kraus

	return kraus_to_superoperator(kraus(operation))


class QuantumStateRepresentation:
	"""Base class for the concrete representations of a quantum state."""

	def state_vector(self) -> list[complex]:
		raise NotImplementedError

	def density_matrix(self) -> Matrix:
		raise NotImplementedError


class QuantumState(QuantumStateRepresentation):
	"""A quantum state, backed by a state vector or a density matrix."""

	def __init__(self, state: Any, *, qid_shape: Sequence[int] | None = None, atol: float = 1e-7) -> None:
		if isinstance(state, QuantumStateRepresentation):
			self._representation: QuantumStateRepresentation | None = state
			self._vector: list[complex] | None = None
			self._density: Matrix | None = None
			return
		self._representation = None
		if _is_density(state):
			self._vector = None
			self._density = to_valid_density_matrix(state, qid_shape or (2,) * (len(state).bit_length() - 1), atol=atol)
		else:
			self._density = None
			shape = qid_shape or (2,) * ((len(state) if not isinstance(state, int) else 1).bit_length() - 1)
			self._vector = to_valid_state_vector(state, shape, atol=atol)

	def state_vector(self) -> list[complex]:
		if self._vector is not None:
			return list(self._vector)
		if self._representation is not None:
			return self._representation.state_vector()
		rho = self._density
		assert rho is not None
		return [complex(math.sqrt(max(rho[i][i].real, 0.0))) for i in range(len(rho))]

	def density_matrix(self) -> Matrix:
		if self._density is not None:
			return [list(row) for row in self._density]
		if self._representation is not None:
			return self._representation.density_matrix()
		vector = self._vector
		assert vector is not None
		return density_matrix_from_state_vector(vector)


def _is_density(value: Any) -> bool:
	return isinstance(value, Sequence) and len(value) > 0 and isinstance(value[0], Sequence)


def quantum_state(
	state: Any, *, qid_shape: Sequence[int] | None = None, dtype: Any = None, atol: float = 1e-7
) -> QuantumState:
	"""Return a :class:`QuantumState` for ``state``."""
	del dtype
	return QuantumState(state, qid_shape=qid_shape, atol=atol)


@value_equality
class StabilizerState(QuantumStateRepresentation):
	"""A stabilizer state, backed by a Clifford tableau."""

	def __init__(self, clifford_tableau: Any) -> None:
		self._tableau = clifford_tableau

	@property
	def clifford_tableau(self) -> Any:
		return self._tableau

	def state_vector(self) -> list[complex]:
		return list(self._tableau.state_vector())

	def density_matrix(self) -> Matrix:
		return density_matrix_from_state_vector(self.state_vector())

	def _value_equality_values_(self) -> Any:
		return (tuple(self.state_vector()),)

	def __repr__(self) -> str:
		return f"vqubit.StabilizerState({self._tableau!r})"


__all__ = [
	"QUANTUM_STATE_LIKE",
	"STATE_VECTOR_LIKE",
	"QuantumState",
	"QuantumStateRepresentation",
	"StabilizerState",
	"choi_to_kraus",
	"choi_to_superoperator",
	"density_matrix_from_state_vector",
	"dirac_notation",
	"entanglement_fidelity",
	"eye_tensor",
	"kraus_to_choi",
	"kraus_to_superoperator",
	"one_hot",
	"operation_to_choi",
	"operation_to_superoperator",
	"quantum_state",
	"superoperator_to_choi",
	"superoperator_to_kraus",
	"validate_density_matrix",
	"validate_indices",
	"validate_normalized_state_vector",
	"validate_qid_shape",
]
