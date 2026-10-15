"""Optional interoperability with numpy and sympy.

vqubit has no runtime dependencies: the core works with plain Python lists and
its own :class:`~vqubit.Parameter`. When ``numpy`` (or ``sympy``) happens to be
installed, the public API can hand back ``numpy.ndarray`` objects (and accept
``sympy`` symbols) so code written against Cirq runs unchanged. This module is
the single place that knows whether those packages are importable.

Nothing here is required: :data:`HAS_NUMPY` and :data:`HAS_SYMPY` are ``False``
on a bare install, and every helper degrades to the dependency-free behavior.
"""

from __future__ import annotations

import importlib
from typing import Any

_np: Any = None
_sympy: Any = None

try:  # pragma: no cover - exercised only when numpy is installed
	_np = importlib.import_module("numpy")
	HAS_NUMPY = True
except ImportError:  # pragma: no cover - the default, dependency-free path
	HAS_NUMPY = False

try:  # pragma: no cover - exercised only when sympy is installed
	_sympy = importlib.import_module("sympy")
	HAS_SYMPY = True
except ImportError:  # pragma: no cover - the default, dependency-free path
	HAS_SYMPY = False


def is_numpy_array(value: Any) -> bool:
	"""Report whether ``value`` is a numpy array (always ``False`` without numpy)."""
	return HAS_NUMPY and isinstance(value, _np.ndarray)


def asarray(matrix: Any) -> Any:
	"""Return ``matrix`` as a numpy array when numpy is available, else a list."""
	rows = [[complex(z) for z in row] for row in matrix]
	if HAS_NUMPY:
		return _np.array(rows, dtype=complex)
	return rows


def asvector(vector: Any) -> Any:
	"""Return ``vector`` as a numpy array when numpy is available, else a list."""
	entries = [complex(z) for z in vector]
	if HAS_NUMPY:
		return _np.array(entries, dtype=complex)
	return entries


def to_nested_list(matrix: Any) -> list[list[complex]]:
	"""Return a nested list of ``complex`` from a list or a numpy array."""
	if is_numpy_array(matrix):
		return [[complex(z) for z in row] for row in matrix.tolist()]
	return [[complex(z) for z in row] for row in matrix]


__all__ = [
	"HAS_NUMPY",
	"HAS_SYMPY",
	"asarray",
	"asvector",
	"is_numpy_array",
	"to_nested_list",
]
