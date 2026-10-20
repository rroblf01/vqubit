"""Tests for the cross-entropy benchmarking and readout-error helpers.

These are numerical utilities, so the tests pin down the numbers a reader would
check by hand (a uniform distribution is zero fidelity, a peaked one is one) and
the cases the helpers refuse.
"""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close


def test_linear_xeb_fidelity_from_probabilities() -> None:
	# The estimator works on the probabilities of the *sampled* bitstrings: a
	# uniform sample carries no fidelity, a single always-sampled one all of it.
	dimension = 4
	assert vqubit.linear_xeb_fidelity_from_probabilities(dimension, [0.25] * dimension) == pytest.approx(0.0)
	assert vqubit.linear_xeb_fidelity_from_probabilities(dimension, [1.0]) == pytest.approx(1.0)
	assert vqubit.linear_xeb_fidelity_from_probabilities(dimension, []) == 0.0
	with pytest.raises(ValueError):
		vqubit.linear_xeb_fidelity_from_probabilities(1, [1.0])


def test_log_xeb_fidelity_from_probabilities() -> None:
	# log(dimension * geometric mean) is zero for the uniform case.
	dimension = 4
	assert vqubit.log_xeb_fidelity_from_probabilities(dimension, [0.25] * dimension) == pytest.approx(0.0, abs=1e-12)
	assert vqubit.log_xeb_fidelity_from_probabilities(dimension, []) == 0.0
	# The zeros are skipped rather than taking the log of zero.
	assert vqubit.log_xeb_fidelity_from_probabilities(dimension, [1.0, 0.0, 0.0, 0.0]) == pytest.approx(math.log(4.0))


def test_the_xeb_helpers_take_bitstrings_or_probabilities() -> None:
	dimension = 4
	probabilities: list[float] = [0.25] * 4
	# Bitstrings are converted to probabilities by dividing by the dimension.
	assert vqubit.linear_xeb_fidelity([0, 0, 0, 0], dimension) == pytest.approx(-1 / 3)
	assert vqubit.linear_xeb_fidelity(probabilities, dimension, bitstrings_are_probabilities=True) == pytest.approx(0.0)
	assert vqubit.log_xeb_fidelity(probabilities, dimension, bitstrings_are_probabilities=True) == pytest.approx(
		0.0, abs=1e-12
	)
	assert vqubit.xeb_fidelity(probabilities, dimension, bitstrings_are_probabilities=True) == pytest.approx(0.0)


def test_hog_score() -> None:
	assert vqubit.hog_score_xeb_fidelity_from_probabilities(4, []) == 0.0
	# total / (total + (d - 1) * first probability)
	assert vqubit.hog_score_xeb_fidelity_from_probabilities(4, [1.0, 0.0]) == pytest.approx(0.25)
	assert vqubit.hog_score_xeb_fidelity_from_probabilities(1, [1.0]) == pytest.approx(1.0)
	assert vqubit.hog_score_xeb_fidelity_from_probabilities(4, [0.0, 0.5]) == 0.0


def _readout_result(zeros: list[list[int]], ones: list[list[int]]) -> vqubit.Result:
	return vqubit.Result({"0": zeros, "1": ones})


def test_estimate_single_qubit_readout_errors() -> None:
	# Qubit 0 never misreads; qubit 1 reads 1 as 1 always and 0 as 1 half the time.
	result = _readout_result([[0, 0], [0, 1]], [[0, 1], [0, 1]])
	p0_given_1, p1_given_0 = vqubit.estimate_single_qubit_readout_errors(result)
	assert p1_given_0 == pytest.approx([0.0, 0.5])
	assert p0_given_1 == pytest.approx([1.0, 0.0])


def test_estimate_parallel_single_qubit_readout_errors() -> None:
	result = _readout_result([[0, 0]], [[1, 0]])
	p0_given_1, p1_given_0 = vqubit.estimate_parallel_single_qubit_readout_errors(result)
	assert p1_given_0 == [0.0, 0.0]
	assert p0_given_1 == [0.0, 1.0]


def test_tensored_confusion_matrices() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	identity = [[1.0, 0.0], [0.0, 1.0]]
	bit_flip = [[0.9, 0.1], [0.1, 0.9]]
	tensored = vqubit.TensoredConfusionMatrices([identity, bit_flip], [q0, q1], repetitions=10)
	assert tensored.confusion_matrices == [identity, bit_flip]
	assert tensored.qubits == (q0, q1)
	assert tensored.repetitions == 10
	# The tensor product of the two: the corner entries are the products.
	assert_matrix_close(
		tensored.confusion_matrix,
		[
			[0.9, 0.1, 0.0, 0.0],
			[0.1, 0.9, 0.0, 0.0],
			[0.0, 0.0, 0.9, 0.1],
			[0.0, 0.0, 0.1, 0.9],
		],
	)
	assert "2 qubits" in repr(tensored)


def test_measure_confusion_matrix_of_an_ideal_sampler() -> None:
	q0 = vqubit.LineQubit(0)
	confusion = vqubit.measure_confusion_matrix(vqubit.Simulator(seed=0), [q0], repetitions=20)
	assert_matrix_close(confusion.confusion_matrix, [[1.0, 0.0], [0.0, 1.0]])
	assert confusion.repetitions == 20


def test_measure_confusion_matrix_of_a_noisy_sampler() -> None:
	"""A sampler that always reports 1 shows up as an off-diagonal matrix."""

	class _StuckSampler:
		def run(self, program: vqubit.Circuit, repetitions: int = 1, **kwargs: object) -> vqubit.Result:
			return vqubit.Result({"m": [[1]] * repetitions})

	q0 = vqubit.LineQubit(0)
	confusion = vqubit.measure_confusion_matrix(_StuckSampler(), [q0], repetitions=5)
	assert_matrix_close(confusion.confusion_matrix, [[0.0, 0.0], [1.0, 1.0]])
