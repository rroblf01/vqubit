"""Cross-entropy benchmarking and readout-error estimation.

These mirror the numerical helpers Cirq exposes under ``cirq.experiments``.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from typing import Any

from ._linalg import Matrix
from ._qubits import Qid


def linear_xeb_fidelity_from_probabilities(hilbert_space_dimension: int, probabilities: Sequence[float]) -> float:
	"""Linear XEB fidelity estimate from measured probabilities."""
	if hilbert_space_dimension <= 1:
		raise ValueError("hilbert_space_dimension must be greater than 1")
	if not probabilities:
		return 0.0
	mean = sum(probabilities) / len(probabilities)
	return (hilbert_space_dimension * mean - 1) / (hilbert_space_dimension - 1)


def log_xeb_fidelity_from_probabilities(hilbert_space_dimension: int, probabilities: Sequence[float]) -> float:
	"""Log XEB fidelity estimate from measured probabilities."""
	if not probabilities:
		return 0.0
	mean_log = sum(math.log(p) for p in probabilities if p > 0) / len(probabilities)
	return mean_log + math.log(hilbert_space_dimension)


def _bitstrings_to_probabilities(
	bitstrings: Iterable[int], hilbert_space_dimension: int, bitstrings_are_probabilities: bool
) -> list[float]:
	values = list(bitstrings)
	if bitstrings_are_probabilities:
		return [float(value) for value in values]
	return [value / hilbert_space_dimension for value in values]


def linear_xeb_fidelity(
	bitstrings: Iterable[int],
	hilbert_space_dimension: int,
	*,
	bitstrings_are_probabilities: bool = False,
) -> float:
	"""Linear XEB fidelity estimate from bitstrings or probabilities."""
	probabilities = _bitstrings_to_probabilities(bitstrings, hilbert_space_dimension, bitstrings_are_probabilities)
	return linear_xeb_fidelity_from_probabilities(hilbert_space_dimension, probabilities)


def log_xeb_fidelity(
	bitstrings: Iterable[int],
	hilbert_space_dimension: int,
	*,
	bitstrings_are_probabilities: bool = False,
) -> float:
	"""Log XEB fidelity estimate from bitstrings or probabilities."""
	probabilities = _bitstrings_to_probabilities(bitstrings, hilbert_space_dimension, bitstrings_are_probabilities)
	return log_xeb_fidelity_from_probabilities(hilbert_space_dimension, probabilities)


def xeb_fidelity(
	bitstrings: Iterable[int],
	hilbert_space_dimension: int,
	*,
	bitstrings_are_probabilities: bool = False,
) -> float:
	"""Linear XEB fidelity estimate (the default XEB estimator)."""
	return linear_xeb_fidelity(
		bitstrings, hilbert_space_dimension, bitstrings_are_probabilities=bitstrings_are_probabilities
	)


def hog_score_xeb_fidelity_from_probabilities(hilbert_space_dimension: int, probabilities: Sequence[float]) -> float:
	"""HOG-score XEB fidelity estimate from probabilities."""
	if not probabilities:
		return 0.0
	total = sum(probabilities)
	return total / (total + (hilbert_space_dimension - 1) * probabilities[0]) if probabilities[0] else 0.0


def estimate_single_qubit_readout_errors(
	measurement_results: Any,
) -> tuple[list[float], list[float]]:
	"""Estimate per-qubit readout errors from measurements of ``|0...0>`` and ``|1...1>``."""
	measurements = measurement_results.measurements
	zero_rows = measurements["0"]
	one_rows = measurements["1"]
	num_qubits = len(zero_rows[0])
	p1_given_0 = [sum(row[q] for row in zero_rows) / len(zero_rows) for q in range(num_qubits)]
	p0_given_1 = [1 - sum(row[q] for row in one_rows) / len(one_rows) for q in range(num_qubits)]
	return p0_given_1, p1_given_0


def estimate_parallel_single_qubit_readout_errors(
	measurement_results: Any,
) -> tuple[list[float], list[float]]:
	"""Estimate per-qubit readout errors when all qubits are measured together."""
	return estimate_single_qubit_readout_errors(measurement_results)


class TensoredConfusionMatrices:
	"""Per-qubit confusion matrices, tensored into a single matrix."""

	def __init__(
		self,
		confusion_matrices: Sequence[Sequence[Sequence[float]]],
		qubits: Sequence[Qid],
		repetitions: int,
	) -> None:
		self._confusion_matrices = [[list(row) for row in matrix] for matrix in confusion_matrices]
		self.qubits = tuple(qubits)
		self.repetitions = int(repetitions)

	@property
	def confusion_matrices(self) -> list[list[list[float]]]:
		return [[list(row) for row in matrix] for matrix in self._confusion_matrices]

	@property
	def confusion_matrix(self) -> Matrix:
		result: Matrix = [[1 + 0j]]
		for matrix in self._confusion_matrices:
			result = _kron_real(result, matrix)
		return result

	def __repr__(self) -> str:
		return f"vqubit.TensoredConfusionMatrices({len(self.qubits)} qubits)"


def _kron_real(a: Sequence[Sequence[Any]], b: Sequence[Sequence[float]]) -> Matrix:
	rows = len(a)
	cols = len(a[0])
	b_rows = len(b)
	b_cols = len(b[0])
	result: Matrix = [[0j] * (cols * b_cols) for _ in range(rows * b_rows)]
	for i in range(rows):
		for j in range(cols):
			for k in range(b_rows):
				for m in range(b_cols):
					result[i * b_rows + k][j * b_cols + m] = complex(a[i][j]) * b[k][m]
	return result


def measure_confusion_matrix(
	sampler: Any, qubits: Sequence[Qid], repetitions: int = 1000, **kwargs: Any
) -> TensoredConfusionMatrices:
	"""Measure the readout confusion matrix of ``qubits`` on ``sampler``."""
	from ._circuits import Circuit
	from ._gates import X, measure

	confusion_matrices: list[list[list[float]]] = []
	for qubit in qubits:
		matrix = [[0.0, 0.0], [0.0, 0.0]]
		for prepared in (0, 1):
			circuit = Circuit()
			if prepared == 1:
				circuit.append(X(qubit))
			circuit.append(measure(qubit, key="m"))
			result = sampler.run(circuit, repetitions=repetitions, **kwargs)
			rows = result.measurements["m"]
			for observed in (0, 1):
				matrix[observed][prepared] = sum(1 for row in rows if row[0] == observed) / len(rows)
		confusion_matrices.append(matrix)
	return TensoredConfusionMatrices(confusion_matrices, qubits, repetitions)


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
