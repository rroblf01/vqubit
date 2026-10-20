"""Tests for the dependency-free code paths.

The published wheel runs without numpy or sympy, so the pure-Python fallbacks
are the code most users execute. These tests force them (by hiding the fast
paths) and check that they agree with the fast ones, and they cover the small
matrix/vector helpers that the fast path skips.
"""

from __future__ import annotations

import pytest

import vqubit
from tests.helpers import assert_matrix_close, assert_vector_close
from vqubit import _compat, _linalg

MATRIX = [[1 + 2j, 3 - 4j], [5j, 6 + 0j]]
VECTOR = [1 + 2j, 3 - 4j]


def test_flat_conversions_without_the_fast_path(monkeypatch: pytest.MonkeyPatch) -> None:
	# Hide every numpy shortcut: the pure-Python loops have to take over.
	monkeypatch.setattr(_compat, "flat_from_matrix", lambda matrix: None)
	monkeypatch.setattr(_compat, "matrix_from_flat", lambda flat, dim: None)
	monkeypatch.setattr(_compat, "flat_from_vector", lambda state: None)
	monkeypatch.setattr(_compat, "vector_from_flat", lambda flat: None)

	flat = _linalg.to_flat(MATRIX)
	assert flat == [1.0, 2.0, 3.0, -4.0, 0.0, 5.0, 6.0, 0.0]
	assert_matrix_close(_linalg.from_flat(flat, 2), MATRIX)
	assert _linalg.vector_to_flat(VECTOR) == [1.0, 2.0, 3.0, -4.0]
	assert_vector_close(_linalg.vector_from_flat([1.0, 2.0, 3.0, -4.0]), VECTOR)


def test_flat_conversions_with_the_fast_path() -> None:
	# The numpy path has a minimum size below which it is not worth the round
	# trip, so exercise both sides of it.
	small = [[1 + 0j, 0j], [0j, 1 + 0j]]
	assert _linalg.to_flat(small) == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
	large = [[complex(i, j) for j in range(4)] for i in range(4)]
	flat = _linalg.to_flat(large)
	assert_matrix_close(_linalg.from_flat(flat, 4), large)


def test_apply_gate_validates_the_dimension() -> None:
	state = [1.0, 0.0]
	with pytest.raises(ValueError, match="does not match"):
		_linalg.apply_gate(state, vqubit.unitary(vqubit.CNOT), [0])
	# And applies a two-qubit gate when the dimensions agree.
	assert_vector_close(_linalg.apply_gate([1, 0, 0, 0], vqubit.unitary(vqubit.CNOT), [0, 1]), [1, 0, 0, 0])


def test_linalg_predicates_and_helpers() -> None:
	assert _linalg.is_unitary([[0, 1], [1, 0]])
	assert not _linalg.is_unitary([[1, 0], [0, 2]])
	assert_matrix_close(_linalg.conjugate([[1 + 1j, 2 - 3j]]), [[1 - 1j, 2 + 3j]])
	assert_matrix_close(_linalg.transpose([[1, 2], [3, 4]]), [[1, 3], [2, 4]])
	assert _linalg.det([[1, 2], [3, 4]]) == pytest.approx(-2)


def test_target_gateset_rejects_a_non_unitary_single_qubit_gate() -> None:
	# The gatesets only accept unitary single-qubit gates for merging.
	q0 = vqubit.LineQubit(0)
	channel = vqubit.amplitude_damp(0.1)
	gateset = vqubit.CZTargetGateset()
	with pytest.raises(TypeError, match="unitary single-qubit gate"):
		gateset.compile_operation(channel.on(q0))
	# A unitary single-qubit gate compiles into a PhasedXZ.
	compiled = gateset.compile_operation(vqubit.H(q0))
	assert len(compiled) == 1
	assert vqubit.linalg.allclose_up_to_global_phase(vqubit.unitary(vqubit.Circuit(compiled)), vqubit.unitary(vqubit.H))
