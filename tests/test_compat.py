"""Tests for the optional numpy/sympy interoperability layer."""

from __future__ import annotations

import vqubit
from vqubit import _compat


def test_compat_flags_are_booleans() -> None:
	assert isinstance(_compat.HAS_NUMPY, bool)
	assert isinstance(_compat.HAS_SYMPY, bool)


def test_unitary_return_type_follows_compat() -> None:
	result = vqubit.unitary(vqubit.X)
	if _compat.HAS_NUMPY:
		import importlib

		np = importlib.import_module("numpy")
		assert isinstance(result, np.ndarray)
	else:
		assert isinstance(result, list)


def test_matrix_values_are_correct_in_either_mode() -> None:
	matrix = vqubit.unitary(vqubit.X)
	assert matrix[0][1] == 1
	assert matrix[1][0] == 1
	assert matrix[0][0] == 0


def test_kraus_return_type_follows_compat() -> None:
	operators = vqubit.kraus(vqubit.X)
	assert len(operators) == 1
	if _compat.HAS_NUMPY:
		import importlib

		np = importlib.import_module("numpy")
		assert isinstance(operators[0], np.ndarray)
	else:
		assert isinstance(operators[0], list)


def test_linalg_accepts_either_return_type() -> None:
	matrix = vqubit.unitary(vqubit.X)
	assert vqubit.is_unitary(matrix)
	assert vqubit.is_hermitian(matrix)
	assert vqubit.kron(matrix, matrix)[0][3] == 1
