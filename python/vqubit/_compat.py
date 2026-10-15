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

import contextlib
import contextvars
import importlib
from collections.abc import Generator
from typing import Any

# A Cirq-compatible debug flag: when False, validations may be skipped for speed.
__cirq_debug__: contextvars.ContextVar[bool] = contextvars.ContextVar("__cirq_debug__", default=__debug__)


@contextlib.contextmanager
def with_debug(value: bool) -> Generator[None, None, None]:
	"""Set ``__cirq_debug__`` for the duration of the context, then restore it."""
	token = __cirq_debug__.set(value)
	try:
		yield
	finally:
		__cirq_debug__.reset(token)


_np: Any = None
_sympy: Any = None

# Matrices up to this side use the pure-Python flatten, which beats the numpy
# round trip for a handful of entries.
_NUMPY_FLATTEN_MIN = 8

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
	if HAS_NUMPY:
		if isinstance(matrix, _np.ndarray):
			return matrix
		return _np.array([[complex(z) for z in row] for row in matrix], dtype=complex)
	return [[complex(z) for z in row] for row in matrix]


def asvector(vector: Any) -> Any:
	"""Return ``vector`` as a numpy array when numpy is available, else a list."""
	if HAS_NUMPY:
		if isinstance(vector, _np.ndarray):
			return vector
		return _np.array([complex(z) for z in vector], dtype=complex)
	return [complex(z) for z in vector]


def _real(value: Any) -> float:
	return float(value.real) if isinstance(value, complex) else float(value)


def as_real_vector(vector: Any) -> Any:
	"""Return ``vector`` as a real numpy array (or list), like Cirq's ``kak_vector``."""
	values = [_real(z) for z in vector]
	if HAS_NUMPY:
		return _np.array(values, dtype=float)
	return values


def as_real_matrix(rows: Any) -> Any:
	"""Return a list of rows as a real numpy array (or list of lists)."""
	values = [[_real(z) for z in row] for row in rows]
	if HAS_NUMPY:
		return _np.array(values, dtype=float)
	return values


def to_nested_list(matrix: Any) -> list[list[complex]]:
	"""Return a nested list of ``complex`` from a list or a numpy array."""
	if is_numpy_array(matrix):
		return [[complex(z) for z in row] for row in matrix.tolist()]
	return [[complex(z) for z in row] for row in matrix]


def flat_from_matrix(matrix: Any) -> list[float] | None:
	"""Flatten a complex matrix to the native interleaved form, fast with numpy.

	Small gate matrices are left to the pure-Python fallback, which is faster
	than paying the numpy round trip for a handful of entries.
	"""
	if not HAS_NUMPY or len(matrix) <= _NUMPY_FLATTEN_MIN:
		return None
	return _np.asarray(matrix, dtype=_np.complex128).ravel().view(_np.float64).tolist()


def flat_from_vector(vector: Any) -> list[float] | None:
	"""Flatten a complex vector to the native interleaved form, fast with numpy."""
	if not HAS_NUMPY or len(vector) <= _NUMPY_FLATTEN_MIN:
		return None
	return _np.asarray(vector, dtype=_np.complex128).ravel().view(_np.float64).tolist()


def matrix_from_flat(flat: Any, dim: int) -> list[list[complex]] | None:
	"""Rebuild a nested complex matrix from the flat buffer, fast with numpy."""
	if not HAS_NUMPY:
		return None
	return _np.asarray(flat, dtype=_np.float64).view(_np.complex128).reshape(dim, dim).tolist()


def vector_from_flat(flat: Any) -> list[complex] | None:
	"""Rebuild a complex vector from the flat buffer, fast with numpy."""
	if not HAS_NUMPY:
		return None
	return _np.asarray(flat, dtype=_np.float64).view(_np.complex128).tolist()


def native_complex_view(obj: Any) -> Any:
	"""Return a zero-copy numpy complex view of a native ``@[vc_buffer]`` object."""
	if not HAS_NUMPY:
		return None
	view = memoryview(obj)
	if view.format == "B":
		view = view.cast("d")
	return _np.frombuffer(view, dtype=_np.float64).view(_np.complex128)


__all__ = [
	"HAS_NUMPY",
	"HAS_SYMPY",
	"asarray",
	"asvector",
	"flat_from_matrix",
	"flat_from_vector",
	"is_numpy_array",
	"matrix_from_flat",
	"native_complex_view",
	"to_nested_list",
	"vector_from_flat",
]
