"""The public linear-algebra helpers of :mod:`vqubit`.

These mirror the functions Cirq exposes under ``cirq.linalg`` (and re-exports at
the top level). Matrices are nested lists of :class:`complex`; the heavy
operations reuse the native kernels in :mod:`vqubit._linalg`. When numpy is
installed the results come back as ``numpy.ndarray`` (see :mod:`vqubit._compat`).
"""

from __future__ import annotations

import cmath
import itertools
import math
import random
from collections.abc import Callable, Sequence
from typing import Any

import vqubit_native as _native

from . import _linalg as _la
from ._compat import as_real_matrix, as_real_vector, asarray, asvector, is_numpy_array, to_nested_list

Matrix = list[list[complex]]
Vector = list[complex]

PAULI_BASIS: dict[str, Matrix] = {
	"I": [[1 + 0j, 0j], [0j, 1 + 0j]],
	"X": [[0j, 1 + 0j], [1 + 0j, 0j]],
	"Y": [[0j, -1j], [1j, 0j]],
	"Z": [[1 + 0j, 0j], [0j, -1 + 0j]],
}

CONTROL_TAG: Matrix = [[complex(float("nan")), 0j], [0j, 1 + 0j]]


def _eye(n: int) -> Matrix:
	return [[1 + 0j if i == j else 0j for j in range(n)] for i in range(n)]


def _zeros(rows: int, cols: int) -> Matrix:
	return [[0j] * cols for _ in range(rows)]


def _mm(a: Any, b: Any) -> Matrix:
	"""Matrix (or matrix-vector) product for arbitrary shapes."""
	a = _as_matrix_or_vector(a)
	b = _as_matrix_or_vector(b)
	if b and not isinstance(b[0], (list, tuple)):
		return [sum(a[i][j] * b[j] for j in range(len(b))) for i in range(len(a))]
	inner = len(b)
	cols = len(b[0])
	return [[sum(a[i][t] * b[t][j] for t in range(inner)) for j in range(cols)] for i in range(len(a))]


def _transpose(m: Any) -> Matrix:
	rows = len(m)
	cols = len(m[0]) if rows else 0
	return [[complex(m[i][j]) for i in range(rows)] for j in range(cols)]


def _conj_t(m: Any) -> Matrix:
	rows = len(m)
	cols = len(m[0]) if rows else 0
	return [[complex(m[j][i]).conjugate() for j in range(rows)] for i in range(cols)]


def _allclose(a: Any, b: Any, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
	a = _as_matrix_or_vector(a)
	b = _as_matrix_or_vector(b)
	if len(a) != len(b):
		return False
	for row_a, row_b in zip(a, b, strict=True):
		left: Any = row_a
		right: Any = row_b
		if isinstance(left, list):
			if not _allclose(left, right, rtol, atol):
				return False
		elif abs(left - right) > atol + rtol * abs(right):
			return False
	return True


# --- tolerance ---------------------------------------------------------------


def all_near_zero(a: Any, *, atol: float = 1e-8) -> bool:
	"""Return whether every element of ``a`` is near zero."""
	if isinstance(a, (list, tuple)):
		return all(all_near_zero(x, atol=atol) for x in a)
	return abs(a) <= atol


def all_near_zero_mod(a: Any, period: float, *, atol: float = 1e-8) -> bool:
	"""Return whether every element of ``a`` is near a multiple of ``period``."""
	if isinstance(a, (list, tuple)):
		return all(all_near_zero_mod(x, period, atol=atol) for x in a)
	return abs((a + period / 2) % period - period / 2) <= atol


# --- predicates --------------------------------------------------------------


def is_diagonal(matrix: Sequence[Sequence[complex]], *, atol: float = 1e-8) -> bool:
	"""Return whether ``matrix`` is approximately diagonal."""
	m = to_nested_list(matrix)
	if any(len(row) != len(m) for row in m):
		return False
	return all(abs(m[i][j]) <= atol for i in range(len(m)) for j in range(len(m)) if i != j)


def is_hermitian(matrix: Sequence[Sequence[complex]], *, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
	"""Return whether ``matrix`` is square and equal to its adjoint."""
	m = to_nested_list(matrix)
	return (not m or len(m) == len(m[0])) and _allclose(m, _conj_t(m), rtol, atol)


def is_orthogonal(matrix: Sequence[Sequence[complex]], *, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
	"""Return whether ``matrix`` is square, real, and its transpose is its inverse."""
	m = to_nested_list(matrix)
	if not m or len(m) != len(m[0]) or any(z.imag != 0 for row in m for z in row):
		return not m
	return _allclose(_mm(m, _transpose(m)), _eye(len(m)), rtol, atol)


def is_special_orthogonal(matrix: Sequence[Sequence[complex]], *, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
	"""Return whether ``matrix`` is orthogonal with determinant one."""
	m = to_nested_list(matrix)
	return is_orthogonal(m, rtol=rtol, atol=atol) and (len(m) == 0 or abs(_la.det(m) - 1) <= atol + rtol)


def is_unitary(matrix: Sequence[Sequence[complex]], *, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
	"""Return whether ``matrix`` is square and its adjoint is its inverse."""
	m = to_nested_list(matrix)
	if not m:
		return True
	if len(m) != len(m[0]):
		return False
	return _allclose(_mm(m, _conj_t(m)), _eye(len(m)), rtol, atol)


def is_special_unitary(matrix: Sequence[Sequence[complex]], *, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
	"""Return whether ``matrix`` is unitary with unit determinant."""
	m = to_nested_list(matrix)
	return is_unitary(m, rtol=rtol, atol=atol) and (len(m) == 0 or abs(_la.det(m) - 1) <= atol + rtol)


def is_normal(matrix: Sequence[Sequence[complex]], *, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
	"""Return whether ``matrix`` commutes with its adjoint."""
	m = to_nested_list(matrix)
	return matrix_commutes(m, _conj_t(m), rtol=rtol, atol=atol)


def is_cptp(*, kraus_ops: Sequence[Sequence[Sequence[complex]]], rtol: float = 1e-5, atol: float = 1e-8) -> bool:
	"""Return whether the Kraus operators define a CPTP map."""
	total: Matrix | None = None
	for k in kraus_ops:
		term = _mm(_conj_t(k), k)
		total = (
			term if total is None else [[total[i][j] + term[i][j] for j in range(len(term))] for i in range(len(term))]
		)
	if total is None:
		return True
	return _allclose(total, _eye(len(total)), rtol, atol)


def matrix_commutes(
	m1: Sequence[Sequence[complex]], m2: Sequence[Sequence[complex]], *, rtol: float = 1e-5, atol: float = 1e-8
) -> bool:
	"""Return whether two square matrices commute."""
	a = to_nested_list(m1)
	b = to_nested_list(m2)
	if len(a) != len(a[0]) or len(a) != len(b):
		return False
	return _allclose(_mm(a, b), _mm(b, a), rtol, atol)


def allclose_up_to_global_phase(
	a: Sequence[Sequence[complex]],
	b: Sequence[Sequence[complex]],
	*,
	rtol: float = 1e-5,
	atol: float = 1e-8,
	equal_nan: bool = False,
) -> bool:
	"""Return whether ``a`` equals ``b`` times some global phase."""
	del equal_nan
	x = to_nested_list(a)
	y = to_nested_list(b)
	if len(x) != len(y) or len(x[0]) != len(y[0]):
		return False
	x2, y2 = match_global_phase(x, y)
	return _allclose(x2, y2, rtol, atol)


def match_global_phase(a: Sequence[Sequence[complex]], b: Sequence[Sequence[complex]]) -> tuple[Matrix, Matrix]:
	"""Phase ``a`` and ``b`` so they agree on the phase of one entry."""
	x = to_nested_list(a)
	y = to_nested_list(b)
	if len(x) != len(y) or len(x[0]) != len(y[0]) or not x:
		return x, y

	best = (0, 0)
	best_mag = -1.0
	for i in range(len(y)):
		for j in range(len(y[0])):
			mag = abs(y[i][j])
			if mag > best_mag:
				best_mag = mag
				best = (i, j)

	def dephase(v: complex) -> complex:
		r, i = v.real, v.imag
		if i == 0:
			return -1 + 0j if r < 0 else 1 + 0j
		if r == 0:
			return 1j if i < 0 else -1j
		return cmath.exp(-1j * math.atan2(i, r))

	pa = dephase(x[best[0]][best[1]])
	pb = dephase(y[best[0]][best[1]])
	return (
		[[z * pa for z in row] for row in x],
		[[z * pb for z in row] for row in y],
	)


# --- combinators -------------------------------------------------------------


def _is_matrix(value: Any) -> bool:
	if is_numpy_array(value):
		return value.ndim == 2
	return isinstance(value, (list, tuple)) and bool(value) and isinstance(value[0], (list, tuple))


def _as_matrix_or_vector(value: Any) -> Any:
	result: Any
	if is_numpy_array(value):
		# The two branches have different types, and mypy rejects the ternary.
		if value.ndim == 2:  # noqa: SIM108
			result = to_nested_list(value)
		else:
			result = [complex(z) for z in value]
	elif isinstance(value, (list, tuple)) and value and isinstance(value[0], (list, tuple)):
		result = to_nested_list(value)
	else:
		result = [complex(z) for z in value]
	return result


def _as_vector(value: Any) -> list[complex]:
	"""Turn a flat sequence (or numpy vector) into a list of complex numbers."""
	if is_numpy_array(value):
		return [complex(z) for z in value]
	return [complex(z) for z in value]


def dot(*values: Any) -> Any:
	"""Compute the matrix/vector product of a sequence of values, in order."""
	if not values:
		raise ValueError("vqubit.dot must be called with arguments")
	if len(values) == 1:
		value = values[0]
		return asarray(value) if _is_matrix(value) else asvector(value)
	result: Any = _as_matrix_or_vector(values[0])
	for value in values[1:]:
		result = _mm(result, _as_matrix_or_vector(value))
	return asarray(result) if _is_matrix(result) else asvector(result)


def kron(*factors: Any, shape_len: int = 2) -> Any:
	"""Compute the Kronecker product of a sequence of matrices or scalars."""
	del shape_len
	product: Any = None
	for factor in factors:
		if _is_matrix(factor):
			matrix = to_nested_list(factor)
			if product is None:
				product = matrix
			elif isinstance(product, list):
				product = _la.kron(product, matrix)
			else:
				# A scalar seen first scales the matrix that follows it.
				product = [[entry * product for entry in row] for row in matrix]
		else:
			scalar = complex(factor)
			if product is None:
				product = scalar
			elif isinstance(product, list):
				product = [[entry * scalar for entry in row] for row in product]
			else:
				product = product * scalar
	if product is None:
		return asarray([[1 + 0j]])
	return asarray(product) if isinstance(product, list) else product


def kron_with_controls(*factors: Any) -> Any:
	"""Kronecker product where ``CONTROL_TAG`` entries become identity blocks."""
	product = to_nested_list(kron(*factors))
	for i in range(len(product)):
		for j in range(len(product)):
			if cmath.isnan(product[i][j]):
				product[i][j] = 1 + 0j if i == j else 0j
	return asarray(product)


def block_diag(*blocks: Any) -> Any:
	"""Concatenate square matrices into a block-diagonal matrix."""
	mats = [to_nested_list(b) for b in blocks]
	for b in mats:
		if b and len(b) != len(b[0]):
			raise ValueError("Blocks must be square.")
	if not any(mats):
		return asarray([])
	n = sum(len(b) for b in mats)
	result = _zeros(n, n)
	offset = 0
	for b in mats:
		for i in range(len(b)):
			for j in range(len(b)):
				result[offset + i][offset + j] = b[i][j]
		offset += len(b)
	return asarray(result)


# --- operator spaces ---------------------------------------------------------


def hilbert_schmidt_inner_product(m1: Any, m2: Any) -> complex:
	"""Return the Hilbert-Schmidt inner product ``tr(m1^dagger m2)``."""
	a = to_nested_list(m1)
	b = to_nested_list(m2)
	return sum(a[i][j].conjugate() * b[i][j] for i in range(len(a)) for j in range(len(a[0])))


def expand_matrix_in_orthogonal_basis(m: Any, basis: dict[str, Any]) -> dict[str, complex]:
	"""Return the coefficients of ``m`` in an orthogonal matrix basis."""
	return {
		name: hilbert_schmidt_inner_product(b, m) / hilbert_schmidt_inner_product(b, b) for name, b in basis.items()
	}


def matrix_from_basis_coefficients(expansion: dict[str, Any], basis: dict[str, Any]) -> Any:
	"""Return the linear combination of basis matrices with given coefficients."""
	some = to_nested_list(next(iter(basis.values())))
	result = _zeros(len(some), len(some[0]))
	for name, coefficient in expansion.items():
		b = to_nested_list(basis[name])
		for i in range(len(result)):
			for j in range(len(result[0])):
				result[i][j] += complex(coefficient) * b[i][j]
	return asarray(result)


def kron_bases(*bases: dict[str, Any], repeat: int = 1) -> dict[str, Matrix]:
	"""Return the tensor product of the given matrix bases."""
	product: dict[str, Matrix] = {"": [[1 + 0j]]}
	for basis in list(bases) * repeat:
		product = {
			name1 + name2: _la.kron(to_nested_list(matrix1), to_nested_list(matrix2))
			for name1, matrix1 in product.items()
			for name2, matrix2 in basis.items()
		}
	return product


def pow_pauli_combination(
	ai: complex, ax: complex, ay: complex, az: complex, exponent: int
) -> tuple[complex, complex, complex, complex]:
	"""Return ``(bi, bx, by, bz)`` with ``(ai I + ax X + ay Y + az Z)^exponent``."""
	if exponent == 0:
		return 1 + 0j, 0j, 0j, 0j
	v = cmath.sqrt(ax * ax + ay * ay + az * az)
	s = (ai + v) ** exponent
	t = (ai - v) ** exponent
	ci = (s + t) / 2
	cxyz = exponent * ai ** (exponent - 1) if abs(s - t) < 1e-15 else (s - t) / (2 * v)
	return ci, cxyz * ax, cxyz * ay, cxyz * az


# --- transformations ---------------------------------------------------------


def reflection_matrix_pow(reflection_matrix: Sequence[Sequence[complex]], exponent: float) -> Any:
	"""Raise a matrix with two opposing eigenvalues to a power."""
	m = to_nested_list(reflection_matrix)
	squared_phase = sum(m[i][0] * m[0][i] for i in range(len(m)))
	phase = cmath.sqrt(squared_phase)
	i_mat = [[phase if i == j else 0j for j in range(len(m))] for i in range(len(m))]
	pos = [[(i_mat[i][j] + m[i][j]) * 0.5 for j in range(len(m))] for i in range(len(m))]
	neg = [[(i_mat[i][j] - m[i][j]) * 0.5 for j in range(len(m))] for i in range(len(m))]
	pos_factor = phase ** (exponent - 1)
	neg_factor = pos_factor * complex(-1) ** exponent
	result = [[pos[i][j] * pos_factor + neg[i][j] * neg_factor for j in range(len(m))] for i in range(len(m))]
	return asarray(result)


def to_special(u: Sequence[Sequence[complex]]) -> Any:
	"""Convert a unitary to a special unitary (unit determinant)."""
	m = to_nested_list(u)
	correction = _la.det(m) ** (-1 / len(m))
	return asarray([[z * correction for z in row] for row in m])


def state_vector_kronecker_product(t1: Sequence[complex], t2: Sequence[complex]) -> Any:
	"""Merge two state vectors into a single unified state vector."""
	a = list(t1)
	b = list(t2)
	return asvector([x * y for x in a for y in b])


def density_matrix_kronecker_product(t1: Sequence[Sequence[complex]], t2: Sequence[Sequence[complex]]) -> Any:
	"""Merge two density matrices into a single unified density matrix."""
	a = to_nested_list(t1)
	b = to_nested_list(t2)
	rows_b = len(b)
	cols_b = len(b[0]) if b else 0
	rows_a = len(a)
	cols_a = len(a[0]) if a else 0
	return asarray(
		[
			[a[i // rows_b][j // cols_b] * b[i % rows_b][j % cols_b] for j in range(cols_a * cols_b)]
			for i in range(rows_a * rows_b)
		]
	)


def _shape_of(tensor: Any) -> tuple[int, ...]:
	if is_numpy_array(tensor):
		return tuple(int(size) for size in tensor.shape)
	shape: list[int] = []
	node = tensor
	while isinstance(node, list):
		shape.append(len(node))
		node = node[0] if node else None
	return tuple(shape)


def _flatten(tensor: Any) -> tuple[list[Any], tuple[int, ...]]:
	shape = _shape_of(tensor)
	flat: list[Any] = []

	def rec(node: Any, depth: int) -> None:
		if depth == len(shape):
			flat.append(node)
			return
		for child in node:
			rec(child, depth + 1)

	rec(tensor, 0)
	return flat, shape


def _unflatten(flat: Sequence[Any], shape: tuple[int, ...]) -> Any:
	if not shape:
		return flat[0]
	step = math.prod(shape[1:])
	return [_unflatten(flat[i * step : (i + 1) * step], shape[1:]) for i in range(shape[0])]


def _strides(shape: Sequence[int]) -> list[int]:
	strides = [1] * len(shape)
	for i in range(len(shape) - 2, -1, -1):
		strides[i] = strides[i + 1] * shape[i + 1]
	return strides


def targeted_left_multiply(
	left_matrix: Sequence[Sequence[complex]],
	right_target: Any,
	target_axes: Sequence[int],
	out: Any = None,
) -> Any:
	"""Left-multiply the given axes of a tensor by a matrix."""
	left = to_nested_list(left_matrix)
	flat, shape = _flatten(right_target)
	k = len(target_axes)
	strides = _strides(shape)
	other_axes = [i for i in range(len(shape)) if i not in target_axes]
	result = [0j] * len(flat)
	for other_vals in itertools.product(*[range(shape[a]) for a in other_axes]):
		base = sum(v * strides[a] for a, v in zip(other_axes, other_vals, strict=True))
		for w_out in range(1 << k):
			acc = 0j
			for w_in in range(1 << k):
				coeff = left[w_out][w_in]
				if coeff == 0:
					continue
				idx = base
				for pos, a in enumerate(target_axes):
					idx += ((w_in >> (k - 1 - pos)) & 1) * strides[a]
				acc += coeff * flat[idx]
			out_idx = base
			for pos, a in enumerate(target_axes):
				out_idx += ((w_out >> (k - 1 - pos)) & 1) * strides[a]
			result[out_idx] = acc
	if len(shape) == 1:
		return asvector(result)
	return asarray(_unflatten(result, shape))


def targeted_conjugate_about(
	tensor: Any,
	target: Any,
	indices: Sequence[int],
	conj_indices: Sequence[int] | None = None,
	buffer: Any = None,
	out: Any = None,
) -> Any:
	"""Compute ``tensor · target · tensor^dagger`` contracted on the given axes."""
	del buffer, out
	if conj_indices is None:
		half = len(_shape_of(target)) // 2
		conj_indices = [i + half for i in indices]
	first = targeted_left_multiply(tensor, target, indices)
	conj_tensor = _conj_elementwise(tensor)
	return targeted_left_multiply(conj_tensor, first, conj_indices)


def _conj_elementwise(tensor: Any) -> Any:
	if is_numpy_array(tensor):
		return tensor.conjugate()
	if isinstance(tensor, list):
		return [_conj_elementwise(x) for x in tensor]
	return complex(tensor).conjugate()


def apply_matrix_to_slices(
	target: Sequence[complex],
	matrix: Sequence[Sequence[complex]],
	slices: Sequence[Any],
	*,
	out: Any = None,
) -> Any:
	"""Left-multiply a matrix onto the given slices of a flat array."""
	del out
	m = to_nested_list(matrix)
	values = [target[s] for s in slices]
	result = list(target)
	for i, s_i in enumerate(slices):
		result[s_i] = sum(m[i][j] * values[j] for j in range(len(values)))
	return asvector(result)


def sub_state_vector(
	state_vector: Any,
	keep_indices: Sequence[int],
	*,
	default: Any = None,
	atol: float = 1e-6,
) -> Any:
	"""Factor ``state_vector`` into the part on ``keep_indices`` and the rest.

	If the state is a product ``|x> ⊗ |y>`` with ``|x>`` living on the kept
	qubits, the returned vector is ``|x>`` (up to a global phase, which is
	randomised so callers do not come to depend on it). If the state is
	entangled across the cut, ``default`` is returned instead — this mirrors
	``cirq.sub_state_vector``, including the tolerance on the reduced state's
	largest eigenvalue.
	"""

	def flatten(node: Any) -> list[Any]:
		if isinstance(node, (list, tuple)) or is_numpy_array(node):
			return [item for child in node for item in flatten(child)]
		return [node]

	vector = [complex(z) for z in flatten(state_vector)]
	size = len(vector)
	if size == 0 or (size & (size - 1)) != 0:
		raise ValueError(f"Input state_vector of size {size} does not represent a state over qubits.")
	n_qubits = size.bit_length() - 1
	keep = [int(index) for index in keep_indices]
	if len(set(keep)) != len(keep):
		raise ValueError(f"keep_indices were {keep} but must be unique.")
	if any(index < 0 or index >= n_qubits for index in keep):
		raise ValueError(f"keep_indices {keep} are an invalid subset of the input state vector.")
	norm = math.sqrt(sum(abs(z) ** 2 for z in vector))
	if abs(norm - 1) > 1e-6:
		raise ValueError("Input state must be normalized.")
	keep_dims = 1 << len(keep)
	if not keep:
		return asvector([complex(1)])
	remaining = [index for index in range(n_qubits) if index not in set(keep)]
	rest_dims = 1 << len(remaining)
	permuted = _permute_qubits(vector, [*sorted(keep), *remaining], n_qubits)
	columns = [permuted[i * rest_dims : (i + 1) * rest_dims] for i in range(keep_dims)]
	# The reduced density matrix on the kept qubits; a product state has it
	# rank one, with eigenvalue 1.
	rho = [
		[sum(columns[i][k] * columns[j][k].conjugate() for k in range(rest_dims)) for j in range(keep_dims)]
		for i in range(keep_dims)
	]
	values, vectors = _hermitian_eig(rho)
	index = max(range(len(values)), key=lambda i: values[i])
	if abs(values[index] - 1) > atol:
		return default
	factored = [vectors[index][i] for i in range(keep_dims)]
	phase = cmath.exp(2j * math.pi * random.random())
	return asvector([phase * z for z in factored])


def _permute_qubits(vector: Sequence[complex], order: Sequence[int], n_qubits: int) -> Vector:
	result = [0j] * len(vector)
	for idx in range(len(vector)):
		bits = [(idx >> (n_qubits - 1 - i)) & 1 for i in range(n_qubits)]
		new_idx = 0
		for _position, qubit in enumerate(order):
			new_idx = (new_idx << 1) | bits[qubit]
		result[new_idx] = vector[idx]
	return result


def partial_trace_of_state_vector_as_mixture(
	state_vector: Sequence[complex], keep_indices: Sequence[int], *, atol: float = 1e-8
) -> tuple[tuple[float, Any], ...]:
	"""Return a mixture representing ``state_vector`` traced down to ``keep_indices``."""
	vector = [complex(z) for z in state_vector]
	n_qubits = int(math.log2(len(vector)))
	keep = list(keep_indices)
	rho = _partial_trace_density(vector, keep, n_qubits)
	vals, vecs = _hermitian_eig(rho)
	mixed = [(float(v.real), asvector(vec)) for v, vec in zip(vals, vecs, strict=True) if abs(v) > atol]
	return tuple(mixed)


def _partial_trace_density(vector: Vector, keep: Sequence[int], n_qubits: int) -> Matrix:
	dim = 1 << len(keep)
	remaining = [i for i in range(n_qubits) if i not in set(keep)]
	rho = _zeros(dim, dim)
	for a in range(dim):
		for b in range(dim):
			acc = 0j
			for r in range(1 << len(remaining)):
				ia = _interleave_bits(a, r, keep, remaining, n_qubits)
				ib = _interleave_bits(b, r, keep, remaining, n_qubits)
				acc += vector[ia] * vector[ib].conjugate()
			rho[a][b] = acc
	return rho


def _interleave_bits(keep_val: int, rest_val: int, keep: Sequence[int], remaining: Sequence[int], n_qubits: int) -> int:
	bits = [0] * n_qubits
	for position, qubit in enumerate(keep):
		bits[qubit] = (keep_val >> (len(keep) - 1 - position)) & 1
	for position, qubit in enumerate(remaining):
		bits[qubit] = (rest_val >> (len(remaining) - 1 - position)) & 1
	idx = 0
	for bit in bits:
		idx = (idx << 1) | bit
	return idx


def _hermitian_eig(matrix: Sequence[Sequence[complex]]) -> tuple[list[float], list[Vector]]:
	"""Eigen-decompose a Hermitian matrix via a real embedding and ``jacobi_eig``.

	The embedding ``[[Re, -Im], [Im, Re]]`` is real symmetric and has each
	eigenvalue of ``matrix`` doubled. Returns the (real) eigenvalues and the
	eigenvectors as columns.
	"""
	m = to_nested_list(matrix)
	n = len(m)
	size = 2 * n
	embed = [[0.0] * size for _ in range(size)]
	for i in range(n):
		for j in range(n):
			re = m[i][j].real
			im = m[i][j].imag
			embed[i][j] = re
			embed[i][n + j] = -im
			embed[n + i][j] = im
			embed[n + i][n + j] = re
	flat = [embed[i][j] for i in range(size) for j in range(size)]
	result = _native.jacobi_eig(flat, size)
	values = list(result["values"])
	columns = list(result["vectors"])
	order = sorted(range(size), key=lambda index: values[index])
	eigenvalues: list[float] = []
	vectors: list[Vector] = []
	i = 0
	while i < size:
		j = i
		while j < size and abs(values[order[j]] - values[order[i]]) < 1e-8:
			j += 1
		group = order[i:j]
		candidates = [
			[complex(columns[k * size + index], columns[(n + k) * size + index]) for k in range(n)] for index in group
		]
		basis: list[Vector] = []
		for candidate in candidates:
			w = candidate[:]
			for b in basis:
				projection = sum(w[k] * b[k].conjugate() for k in range(n))
				w = [w[k] - projection * b[k] for k in range(n)]
			norm = math.sqrt(sum(abs(z) ** 2 for z in w))
			if norm > 1e-8:
				basis.append([z / norm for z in w])
		eigenvalues.extend([values[order[i]]] * len(basis))
		vectors.extend(basis)
		i = j
	return eigenvalues, vectors


# --- diagonalization ---------------------------------------------------------


def _eigh_real(matrix: Sequence[Sequence[complex]]) -> tuple[list[float], list[list[float]]]:
	n = len(matrix)
	flat = [float(matrix[i][j].real) for i in range(n) for j in range(n)]
	result = _native.jacobi_eig(flat, n)
	values = list(result["values"])
	raw = result["vectors"]
	vectors = [[raw[i * n + j] for i in range(n)] for j in range(n)]
	return values, vectors


def diagonalize_real_symmetric_matrix(
	matrix: Sequence[Sequence[complex]],
	*,
	rtol: float = 1e-5,
	atol: float = 1e-8,
	check_preconditions: bool = True,
) -> Any:
	"""Return an orthogonal matrix ``P`` with ``P^T matrix P`` diagonal."""
	m = to_nested_list(matrix)
	if check_preconditions and (
		any(z.imag != 0 for row in m for z in row) or not is_hermitian(m, rtol=rtol, atol=atol)
	):
		raise ValueError("Input must be real and symmetric.")
	_, vectors = _eigh_real(m)
	n = len(m)
	return asarray([[vectors[j][i] for j in range(n)] for i in range(n)])


def _contiguous_groups(length: int, comparator: Callable[[int, int], bool]) -> list[tuple[int, int]]:
	result = []
	start = 0
	while start < length:
		past = start + 1
		while past < length and comparator(start, past):
			past += 1
		result.append((start, past))
		start = past
	return result


def diagonalize_real_symmetric_and_sorted_diagonal_matrices(
	symmetric_matrix: Sequence[Sequence[complex]],
	diagonal_matrix: Sequence[Sequence[complex]],
	*,
	rtol: float = 1e-5,
	atol: float = 1e-8,
	check_preconditions: bool = True,
) -> Any:
	"""Return an orthogonal matrix diagonalizing both commuting matrices."""
	sym = to_nested_list(symmetric_matrix)
	diag = to_nested_list(diagonal_matrix)
	if check_preconditions:
		if any(z.imag != 0 for row in sym for z in row) or not is_hermitian(sym, rtol=rtol, atol=atol):
			raise ValueError("symmetric_matrix must be real symmetric.")
		if not is_diagonal(diag, atol=atol) or not matrix_commutes(diag, sym, rtol=rtol, atol=atol):
			raise ValueError("Given matrices must commute and be diagonal/symmetric.")
	n = len(sym)
	p = _zeros(n, n)
	for i in range(n):
		p[i][i] = 1.0

	def similar(i: int, j: int) -> bool:
		return abs(diag[i][i] - diag[j][j]) <= atol + rtol * abs(diag[j][j])

	for start, end in _contiguous_groups(n, similar):
		block = [row[start:end] for row in sym[start:end]]
		_, vectors = _eigh_real(block)
		for i in range(end - start):
			for j in range(end - start):
				p[start + i][start + j] = vectors[j][i]
	return asarray(p)


def _real_svd(mat: Any) -> tuple[list[list[Any]], list[float], list[list[Any]]]:
	# The matrices reach here as complex numbers with a zero imaginary part.
	mat = [[float(complex(entry).real) for entry in row] for row in mat]
	m = len(mat)
	n = len(mat[0]) if m else 0
	if m == 0 or n == 0:
		return _zeros(0, 0), [], _zeros(0, 0)
	ata = [[sum(mat[k][i] * mat[k][j] for k in range(m)) for j in range(n)] for i in range(n)]
	vals, vecs = _eigh_real(ata)
	pairs = sorted(zip(vals, vecs, strict=True), key=lambda p: -p[0])
	s = [math.sqrt(max(v, 0.0)) for v, _ in pairs]
	v_mat = [[pairs[j][1][i] for j in range(n)] for i in range(n)]
	# Build the left singular vectors one column at a time. A zero singular value
	# (a rank-deficient matrix) gets any unit vector orthogonal to the previous
	# ones, so the result stays an orthonormal basis.
	columns: list[list[float]] = []
	for j in range(n):
		if s[j] > 1e-12:
			column = [sum(mat[i][k] * v_mat[k][j] for k in range(n)) / s[j] for i in range(m)]
		else:
			column = [value.real for value in _first_orthogonal_direction(columns, m)]
			columns.append(column)
			continue
		for previous in columns:
			projection = sum(column[i] * previous[i] for i in range(m))
			column = [column[i] - projection * previous[i] for i in range(m)]
		norm = math.sqrt(sum(value * value for value in column))
		columns.append([value / norm for value in column] if norm > 1e-300 else column)
	u_mat = [[columns[j][i] for j in range(n)] for i in range(m)]
	return u_mat, s, [[v_mat[j][i] for j in range(n)] for i in range(n)]


def bidiagonalize_real_matrix_pair_with_symmetric_products(
	mat1: Sequence[Sequence[complex]],
	mat2: Sequence[Sequence[complex]],
	*,
	rtol: float = 1e-5,
	atol: float = 1e-8,
	check_preconditions: bool = True,
) -> tuple[Any, Any]:
	"""Return orthogonal ``(L, R)`` making both matrices diagonal."""
	a = to_nested_list(mat1)
	b = to_nested_list(mat2)
	if check_preconditions:
		if any(z.imag != 0 for row in a for z in row):
			raise ValueError("mat1 must be real.")
		if any(z.imag != 0 for row in b for z in row):
			raise ValueError("mat2 must be real.")
		if len(a) != len(a[0]) or len(a) != len(b) or any(len(row) != len(b) for row in b):
			raise ValueError("mat1 and mat2 must be square matrices of the same size.")
		if not _allclose(_mm(a, _transpose(b)), _transpose(_mm(a, _transpose(b))), rtol, atol):
			raise ValueError("mat1 @ mat2.T must be symmetric.")
		if not _allclose(_mm(_transpose(a), b), _transpose(_mm(_transpose(a), b)), rtol, atol):
			raise ValueError("mat1.T @ mat2 must be symmetric.")
	base_left, s, base_right = _real_svd([[z.real for z in row] for row in a])
	rank = len(s)
	while rank > 0 and s[rank - 1] <= atol:
		rank -= 1
	base_diag = [[s[i] if i == j else 0.0 for j in range(rank)] for i in range(rank)]
	semi = _mm(_transpose(base_left), [[z.real for z in row] for row in b])
	semi = _mm(semi, _transpose(base_right))
	overlap = [row[:rank] for row in semi[:rank]]
	overlap_adjust = to_nested_list(
		diagonalize_real_symmetric_and_sorted_diagonal_matrices(
			overlap, base_diag, rtol=rtol, atol=atol, check_preconditions=False
		)
	)
	extra = [row[rank:] for row in semi[rank:]]
	extra_left, _, extra_right = _real_svd(extra)
	left_adjust = to_nested_list(block_diag(overlap_adjust, extra_left))
	right_adjust = to_nested_list(block_diag(_transpose(overlap_adjust), extra_right))
	left = _mm(_transpose(left_adjust), _transpose(base_left))
	right = _mm(_transpose(base_right), _transpose(right_adjust))
	return asarray(left), asarray(right)


def bidiagonalize_unitary_with_special_orthogonals(
	mat: Sequence[Sequence[complex]],
	*,
	rtol: float = 1e-5,
	atol: float = 1e-8,
	check_preconditions: bool = True,
) -> tuple[Any, Any, Any]:
	"""Return ``(L, d, R)`` with ``L mat R = diag(d)`` and ``L, R`` special orthogonal."""
	m = to_nested_list(mat)
	if check_preconditions and not is_unitary(m, rtol=rtol, atol=atol):
		raise ValueError("matrix must be unitary.")
	pair = bidiagonalize_real_matrix_pair_with_symmetric_products(
		[[z.real for z in row] for row in m],
		[[z.imag for z in row] for row in m],
		rtol=rtol,
		atol=atol,
		check_preconditions=False,
	)
	left = to_nested_list(pair[0])
	right = to_nested_list(pair[1])
	if _la.det(left).real < 0:
		for j in range(len(left[0])):
			left[0][j] = -left[0][j]
	if _la.det(right).real < 0:
		for i in range(len(right)):
			right[i][0] = -right[i][0]
	diag = _mm(_mm(left, m), right)
	return (
		asarray(left),
		asarray([[diag[i][i] if i == j else 0j for j in range(len(diag))] for i in range(len(diag))]),
		asarray(right),
	)


# --- eigen-decomposition and KAK --------------------------------------------


def _qr(a: Matrix) -> tuple[Matrix, Matrix]:
	n = len(a)
	q = _zeros(n, n)
	r = _zeros(n, n)
	for j in range(n):
		v = [a[i][j] for i in range(n)]
		for k in range(j):
			r[k][j] = sum(q[i][k].conjugate() * v[i] for i in range(n))
			v = [v[i] - r[k][j] * q[i][k] for i in range(n)]
		norm = math.sqrt(sum(abs(x) ** 2 for x in v))
		r[j][j] = norm
		if norm > 1e-300:
			for i in range(n):
				q[i][j] = v[i] / norm
	return q, r


def _normalize_vector(vector: Sequence[complex]) -> Vector:
	norm = math.sqrt(sum(abs(z) ** 2 for z in vector))
	if norm <= 1e-300:
		return [0j] * len(vector)
	return [z / norm for z in vector]


def _null_space_basis(matrix: Sequence[Sequence[complex]], *, atol: float = 1e-9) -> list[Vector]:
	"""Return an orthonormal basis of the null space of ``matrix``.

	Gauss-Jordan elimination with partial pivoting, then Gram-Schmidt on the
	free columns. The matrices we factor are ``A - λI`` for a normal ``A``, so
	the basis vectors are already orthogonal in exact arithmetic.
	"""
	a = [[complex(entry) for entry in row] for row in matrix]
	n = len(a)
	pivots: list[int] = []
	row = 0
	for col in range(n):
		best = max(range(row, n), key=lambda r: abs(a[r][col]))
		if abs(a[best][col]) <= atol:
			continue
		a[row], a[best] = a[best], a[row]
		pivot = a[row][col]
		a[row] = [entry / pivot for entry in a[row]]
		for r in range(n):
			if r != row and abs(a[r][col]) > atol:
				factor = a[r][col]
				a[r] = [entry - factor * reference for entry, reference in zip(a[r], a[row], strict=True)]
		pivots.append(col)
		row += 1
		if row == n:
			break
	free = [col for col in range(n) if col not in set(pivots)]
	basis: list[Vector] = []
	for column in free:
		vector = [0j] * n
		vector[column] = 1 + 0j
		for r, pivot_column in enumerate(pivots):
			vector[pivot_column] = -a[r][column]
		for previous in basis:
			projection = sum(vector[i] * previous[i].conjugate() for i in range(n))
			vector = [vector[i] - projection * previous[i] for i in range(n)]
		if math.sqrt(sum(abs(z) ** 2 for z in vector)) > atol:
			basis.append(_normalize_vector(vector))
	return basis


def _characteristic_coefficients(matrix: Sequence[Sequence[complex]]) -> list[complex]:
	"""The coefficients of ``det(xI - A)``, via Leverrier's recursion."""
	n = len(matrix)
	a = [[complex(entry) for entry in row] for row in matrix]
	traces: list[complex] = []
	power = _eye(n)
	for _ in range(n):
		power = _mm(power, a)
		traces.append(sum(power[i][i] for i in range(n)))
	coefficients = [1 + 0j]
	for k in range(1, n + 1):
		total = sum(coefficients[k - i] * traces[i - 1] for i in range(1, k + 1))
		coefficients.append(-total / k)
	return coefficients


def _polynomial_roots(coefficients: Sequence[complex], iterations: int = 400) -> list[complex]:
	"""The roots of a monic polynomial by Durand-Kerner iteration.

	Unlike a QR sweep, this does not stall when every root has the same modulus
	(a permutation matrix, say), and for the small matrices this library
	decomposes it converges in a handful of steps.
	"""
	n = len(coefficients) - 1
	if n <= 0:
		return []
	if n == 1:
		return [-coefficients[1] / coefficients[0]]
	radius = max(1.0, abs(coefficients[-1]) ** (1.0 / n))
	roots = [cmath.exp(2j * math.pi * (k + 0.25) / n) * radius for k in range(n)]

	def evaluate(value: complex) -> complex:
		result = complex(coefficients[0])
		for coefficient in coefficients[1:]:
			result = result * value + coefficient
		return result

	for _ in range(iterations):
		converged = True
		for i in range(n):
			denominator = 1 + 0j
			for j in range(n):
				if j != i:
					denominator *= roots[i] - roots[j]
			if abs(denominator) <= 1e-300:
				continue
			step = evaluate(roots[i]) / denominator
			roots[i] -= step
			if abs(step) > 1e-13 * max(1.0, abs(roots[i])):
				converged = False
		if converged:
			break
	return roots


def _complex_eig(
	matrix: Sequence[Sequence[complex]], iterations: int = 3000, tol: float = 1e-12
) -> tuple[Vector, Matrix]:
	"""Eigen-decompose a normal matrix.

	The eigenvalues are the roots of the characteristic polynomial, and each
	eigenvector comes from the null space of ``A - λI``. For a normal matrix the
	null spaces of different eigenvalues are orthogonal, so the resulting basis
	is unitary even when an eigenvalue repeats; a repeated eigenvalue whose
	space is already spanned takes any direction left by Gram-Schmidt.
	"""
	del iterations, tol
	m = to_nested_list(matrix)
	n = len(m)
	if n == 0:
		return [], []
	if n == 1:
		return [complex(m[0][0])], [[1 + 0j]]
	values = _polynomial_roots(_characteristic_coefficients(m))
	basis: list[Vector] = []
	for value in values:
		shifted = [[complex(m[i][j]) - (value if i == j else 0j) for j in range(n)] for i in range(n)]
		chosen: Vector | None = None
		for candidate in _null_space_basis(shifted):
			orthogonal = all(
				abs(sum(candidate[i] * previous[i].conjugate() for i in range(n))) < 1e-6 for previous in basis
			)
			if orthogonal:
				chosen = candidate
				break
		if chosen is None:
			chosen = _first_orthogonal_direction(basis, n)
		basis.append(chosen)
	# ``unitary_eig`` hands back a matrix whose columns are the eigenvectors.
	return values, [[basis[j][i] for j in range(n)] for i in range(n)]


def _first_orthogonal_direction(basis: Sequence[Any], n: int) -> Vector:
	"""A unit vector orthogonal to ``basis`` (the standard basis, then Gram-Schmidt)."""
	for index in range(n):
		candidate = [0j] * n
		candidate[index] = 1 + 0j
		for previous in basis:
			projection = sum(candidate[i] * previous[i].conjugate() for i in range(n))
			candidate = [candidate[i] - projection * previous[i] for i in range(n)]
		if math.sqrt(sum(abs(z) ** 2 for z in candidate)) > 1e-8:
			return _normalize_vector(candidate)
	return [0j] * n


def unitary_eig(
	matrix: Sequence[Sequence[complex]], check_preconditions: bool = True, atol: float = 1e-8
) -> tuple[Any, Any]:
	"""Return the eigenvalues and unitary eigenbasis of a normal matrix."""
	m = to_nested_list(matrix)
	if check_preconditions and not is_normal(m, atol=atol):
		raise ValueError("Input must correspond to a normal matrix.")
	if is_hermitian(m, atol=atol):
		hermitian_values, hermitian_vectors = _hermitian_eig(m)
		n = len(m)
		return asvector(hermitian_values), asarray([[hermitian_vectors[j][i] for j in range(n)] for i in range(n)])
	complex_values, complex_vectors = _complex_eig(m)
	return asvector(complex_values), asarray(complex_vectors)


def map_eigenvalues(
	matrix: Sequence[Sequence[complex]], func: Callable[[complex], complex], *, atol: float = 1e-8
) -> Any:
	"""Apply ``func`` to the eigenvalues of a normal (e.g. Hermitian) matrix."""
	m = to_nested_list(matrix)
	n = len(m)
	if is_hermitian(m, atol=atol):
		values, vectors = _hermitian_eig(m)
	else:
		values, basis = unitary_eig(m, atol=atol)
		values = list(values)
		vectors = [[basis[i][j] for i in range(n)] for j in range(n)]
	total = _zeros(n, n)
	for value, vec in zip(values, vectors, strict=True):
		for i in range(n):
			for j in range(n):
				total[i][j] += func(complex(value)) * vec[i] * vec[j].conjugate()
	return asarray(total)


def kron_factor_4x4_to_2x2s(
	matrix: Sequence[Sequence[complex]], rtol: float = 1e-5, atol: float = 1e-8
) -> tuple[complex, Matrix, Matrix]:
	"""Split a 4x4 matrix ``U = g * kron(A, B)`` into ``g``, ``A`` and ``B``."""
	nested = to_nested_list(matrix)
	if len(nested) != 4 or any(len(row) != 4 for row in nested):
		raise ValueError(f"Expected a 4x4 matrix, got shape {len(nested)}x{len(nested[0]) if nested else 0}.")
	result = _native.kron_factor(_la.to_flat(nested))
	global_phase = complex(result["g_re"], result["g_im"])
	f1 = _la.from_flat(result["a"], 2)
	f2 = _la.from_flat(result["b"], 2)
	del rtol, atol
	return global_phase, f1, f2


_MAGIC: Matrix = [
	[math.sqrt(0.5) + 0j, 0j, 0j, 1j * math.sqrt(0.5)],
	[0j, 1j * math.sqrt(0.5), math.sqrt(0.5) + 0j, 0j],
	[0j, 1j * math.sqrt(0.5), -math.sqrt(0.5) + 0j, 0j],
	[math.sqrt(0.5) + 0j, 0j, 0j, -1j * math.sqrt(0.5)],
]
_MAGIC_CONJ_T = _conj_t(_MAGIC)
_KAK_GAMMA: Matrix = [
	[0.25 + 0j, 0.25 + 0j, 0.25 + 0j, 0.25 + 0j],
	[0.25 + 0j, 0.25 + 0j, -0.25 + 0j, -0.25 + 0j],
	[-0.25 + 0j, 0.25 + 0j, -0.25 + 0j, 0.25 + 0j],
	[0.25 + 0j, -0.25 + 0j, -0.25 + 0j, 0.25 + 0j],
]


def so4_to_magic_su2s(
	mat: Sequence[Sequence[complex]],
	*,
	rtol: float = 1e-5,
	atol: float = 1e-8,
	check_preconditions: bool = True,
) -> tuple[Matrix, Matrix]:
	"""Return ``A, B`` in SU(2) with ``Mag^H kron(A, B) Mag = mat``."""
	m = to_nested_list(mat)
	if len(m) != 4 or any(len(row) != 4 for row in m):
		raise ValueError("mat must be 4x4 special orthogonal.")
	if check_preconditions and not is_special_orthogonal(m, rtol=rtol, atol=atol):
		raise ValueError("mat must be 4x4 special orthogonal.")
	ab = _mm(_mm(_MAGIC, m), _MAGIC_CONJ_T)
	_, a, b = kron_factor_4x4_to_2x2s(ab, rtol, atol)
	return a, b


def axis_angle(single_qubit_unitary: Sequence[Sequence[complex]]) -> Any:
	"""Decompose a single-qubit unitary into an axis, angle, and global phase."""
	u = to_nested_list(single_qubit_unitary)
	if not is_unitary(u, atol=1e-8):
		raise ValueError("Input must be a 2x2 unitary.")
	a, b = u[0]
	c, d = u[1]
	wp = (a + d) / 2
	xp = (b + c) / 2j
	yp = (b - c) / 2
	zp = (a - d) / 2j
	p = max((wp, xp, yp, zp), key=abs)
	p = p / abs(p)
	w = min(1.0, max(-1.0, (wp / p).real))
	x = (xp / p).real
	y = (yp / p).real
	z = (zp / p).real
	angle = -2 * math.acos(w)
	n = math.sqrt(x * x + y * y + z * z)
	if n < 1e-7:
		return AxisAngleDecomposition(global_phase=p, angle=0.0, axis=(1.0, 0.0, 0.0))
	return AxisAngleDecomposition(axis=(x / n, y / n, z / n), angle=angle, global_phase=p).canonicalize()


class AxisAngleDecomposition:
	"""A unitary written as ``g * exp(-i theta/2 (x X + y Y + z Z))``."""

	def __init__(self, *, angle: float, axis: tuple[float, float, float], global_phase: complex) -> None:
		norm = math.sqrt(sum(c * c for c in axis))
		if abs(norm - 1) > 1e-8:
			raise ValueError("Axis vector must be normalized.")
		self.global_phase = complex(global_phase)
		self.axis = tuple(float(c) for c in axis)
		self.angle = float(angle)

	def canonicalize(self, atol: float = 1e-8) -> AxisAngleDecomposition:
		"""Return the standard representative with the same unitary."""
		angle = self.angle
		x, y, z = self.axis
		p = self.global_phase
		if x + y + z < 0:
			x, y, z, angle = -x, -y, -z, -angle
		if abs(angle) >= math.pi * 2:
			angle %= math.pi * 4
		while angle <= -math.pi + atol:
			angle += math.pi * 2
			p = -p
		while angle > math.pi + atol:
			angle -= math.pi * 2
			p = -p
		return AxisAngleDecomposition(axis=(x, y, z), angle=angle, global_phase=p)

	def matrix(self) -> Any:
		"""The rotation matrix this decomposition describes."""
		x, y, z = self.axis
		xm = [[0j, 1 + 0j], [1 + 0j, 0j]]
		ym = [[0j, -1j], [1j, 0j]]
		zm = [[1 + 0j, 0j], [0j, -1 + 0j]]
		c = math.cos(-self.angle / 2)
		s = math.sin(-self.angle / 2)
		result = _zeros(2, 2)
		for i in range(2):
			for j in range(2):
				result[i][j] = (
					c * (1 + 0j if i == j else 0j) + 1j * s * (x * xm[i][j] + y * ym[i][j] + z * zm[i][j])
				) * self.global_phase
		return asarray(result)

	_unitary_ = matrix

	def __repr__(self) -> str:
		return (
			f"vqubit.AxisAngleDecomposition(angle={self.angle!r}, axis={self.axis!r}, "
			f"global_phase={self.global_phase!r})"
		)


class KakDecomposition:
	"""A two-qubit operation ``g (a1⊗a0) exp(i(x XX + y YY + z ZZ)) (b1⊗b0)``."""

	def __init__(
		self,
		*,
		global_phase: complex = 1 + 0j,
		single_qubit_operations_before: tuple[Any, Any] | None = None,
		interaction_coefficients: tuple[float, ...],
		single_qubit_operations_after: tuple[Any, Any] | None = None,
	) -> None:
		self.global_phase = complex(global_phase)
		self.single_qubit_operations_before = (
			tuple(to_nested_list(m) for m in single_qubit_operations_before)
			if single_qubit_operations_before
			else (_eye(2), _eye(2))
		)
		self.interaction_coefficients = tuple(float(c) for c in interaction_coefficients)
		self.single_qubit_operations_after = (
			tuple(to_nested_list(m) for m in single_qubit_operations_after)
			if single_qubit_operations_after
			else (_eye(2), _eye(2))
		)

	def matrix(self) -> Any:
		"""The two-qubit unitary this decomposition describes."""
		before = _la.kron(self.single_qubit_operations_before[0], self.single_qubit_operations_before[1])
		after = _la.kron(self.single_qubit_operations_after[0], self.single_qubit_operations_after[1])
		x, y, z = self.interaction_coefficients
		xm = [[0j, 1 + 0j], [1 + 0j, 0j]]
		ym = [[0j, -1j], [1j, 0j]]
		zm = [[1 + 0j, 0j], [0j, -1 + 0j]]

		def interaction(m: Matrix, c: float) -> Matrix:
			return to_nested_list(map_eigenvalues(_la.kron(m, m), lambda v: cmath.exp(1j * v * c)))

		combined = _mm(interaction(zm, z), _mm(interaction(ym, y), interaction(xm, x)))
		result = _mm(_mm(after, combined), before)
		return asarray([[z0 * self.global_phase for z0 in row] for row in result])

	_unitary_ = matrix

	def __repr__(self) -> str:
		return f"vqubit.KakDecomposition(interaction_coefficients={self.interaction_coefficients!r}, global_phase={self.global_phase!r})"


def kak_canonicalize_vector(x: float, y: float, z: float, atol: float = 1e-9) -> KakDecomposition:
	"""Canonicalize an ``XX/YY/ZZ`` interaction into the Weyl chamber."""
	phase = [1 + 0j]
	left = [_eye(2), _eye(2)]
	right = [_eye(2), _eye(2)]
	v = [x, y, z]

	flippers = [
		[[0j, 1j], [1j, 0j]],
		[[0j, 1 + 0j], [-1 + 0j, 0j]],
		[[1j, 0j], [0j, -1j]],
	]
	swappers = [
		[[1j * math.sqrt(0.5), -1j * 1j * math.sqrt(0.5)], [1j * 1j * math.sqrt(0.5), -1j * math.sqrt(0.5)]],
		[[1j * math.sqrt(0.5), 1j * math.sqrt(0.5)], [1j * math.sqrt(0.5), -1j * math.sqrt(0.5)]],
		[[0j, (1 - 1j) * 1j * math.sqrt(0.5)], [(1 + 1j) * 1j * math.sqrt(0.5), 0j]],
	]

	def matpow(m: Matrix, k: int) -> Matrix:
		result = _eye(len(m))
		for _ in range(k):
			result = _mm(result, m)
		return result

	def shift(k: int, step: int) -> None:
		v[k] += step * math.pi / 2
		phase[0] *= 1j**step
		right[0] = _mm(matpow(flippers[k], step % 4), right[0])
		right[1] = _mm(matpow(flippers[k], step % 4), right[1])

	def negate(k1: int, k2: int) -> None:
		v[k1] *= -1
		v[k2] *= -1
		phase[0] *= -1
		s = flippers[3 - k1 - k2]
		left[1] = _mm(left[1], s)
		right[1] = _mm(s, right[1])

	def swap(k1: int, k2: int) -> None:
		v[k1], v[k2] = v[k2], v[k1]
		s = swappers[3 - k1 - k2]
		left[0] = _mm(left[0], s)
		left[1] = _mm(left[1], s)
		right[0] = _mm(s, right[0])
		right[1] = _mm(s, right[1])

	def canonical_shift(k: int) -> None:
		while v[k] <= -math.pi / 4:
			shift(k, 1)
		while v[k] > math.pi / 4:
			shift(k, -1)

	def sort() -> None:
		if abs(v[0]) < abs(v[1]):
			swap(0, 1)
		if abs(v[1]) < abs(v[2]):
			swap(1, 2)
		if abs(v[0]) < abs(v[1]):
			swap(0, 1)

	canonical_shift(0)
	canonical_shift(1)
	canonical_shift(2)
	sort()
	if v[0] < 0:
		negate(0, 2)
	if v[1] < 0:
		negate(1, 2)
	canonical_shift(2)
	if v[0] > math.pi / 4 - atol and v[2] < 0:
		shift(0, -1)
		negate(0, 2)

	return KakDecomposition(
		global_phase=phase[0],
		single_qubit_operations_after=(left[1], left[0]),
		interaction_coefficients=(v[0], v[1], v[2]),
		single_qubit_operations_before=(right[1], right[0]),
	)


def kak_decomposition(
	matrix: Sequence[Sequence[complex]], *, rtol: float = 1e-5, atol: float = 1e-8, check_preconditions: bool = True
) -> KakDecomposition:
	"""Decompose a two-qubit unitary into single-qubit ops and XX/YY/ZZ interactions."""
	from ._kak import kak_decomposition as _tuple_decomposition

	m = to_nested_list(matrix)
	if check_preconditions and (len(m) != 4 or not is_unitary(m, rtol=rtol, atol=atol)):
		raise ValueError("Input must correspond to a 4x4 unitary matrix.")
	phase, (x, y, z), (a1, a2), (b1, b2) = _tuple_decomposition(m, atol=atol)
	inner = kak_canonicalize_vector(x, y, z, atol=atol)
	b1 = _mm(inner.single_qubit_operations_before[0], b1)
	b2 = _mm(inner.single_qubit_operations_before[1], b2)
	a1 = _mm(a1, inner.single_qubit_operations_after[0])
	a2 = _mm(a2, inner.single_qubit_operations_after[1])
	return KakDecomposition(
		interaction_coefficients=inner.interaction_coefficients,
		global_phase=phase * inner.global_phase,
		single_qubit_operations_before=(b1, b2),
		single_qubit_operations_after=(a1, a2),
	)


def _canonicalize_kak_vector(k_vec: list[float], atol: float) -> list[float]:
	v = [(c + math.pi / 4) % (math.pi / 2) - math.pi / 4 for c in k_vec]
	order = sorted(range(3), key=lambda i: abs(v[i]), reverse=True)
	v = [v[i] for i in order]
	if v[0] < 0:
		v[0] *= -1
		v[2] *= -1
	if v[1] < 0:
		v[1] *= -1
		v[2] *= -1
	if abs(v[0] - math.pi / 4) <= atol and v[2] < 0:
		v[2] *= -1
	return v


def kak_vector(
	unitary: Any,
	*,
	rtol: float = 1e-5,
	atol: float = 1e-8,
	check_preconditions: bool = True,
) -> Any:
	"""Compute the canonical KAK vector of a two-qubit unitary."""
	if is_numpy_array(unitary):
		unitary = to_nested_list(unitary)
	if not unitary:
		return asarray([])
	first = unitary[0]
	if isinstance(first[0], (list, tuple)) or is_numpy_array(first[0]):
		return as_real_matrix([_kak_vector_one(u, rtol, atol, check_preconditions) for u in unitary])  # type: ignore[arg-type]
	return as_real_vector(_kak_vector_one(unitary, rtol, atol, check_preconditions))  # type: ignore[arg-type]


def _kak_vector_one(matrix: Any, rtol: float, atol: float, check_preconditions: bool) -> Any:
	from ._kak import kak_decomposition as _tuple_decomposition

	m = to_nested_list(matrix)
	if check_preconditions and (len(m) != 4 or not is_unitary(m, rtol=rtol, atol=atol)):
		raise ValueError("Input must correspond to a 4x4 unitary matrix.")
	_, (x, y, z), _, _ = _tuple_decomposition(m, atol=atol)
	return _canonicalize_kak_vector([x, y, z], atol)


def num_cnots_required(u: Sequence[Sequence[complex]], atol: float = 1e-8) -> int:
	"""Return the minimum number of CNOT/CZ gates needed for a two-qubit unitary."""
	m = to_nested_list(u)
	if len(m) != 4 or len(m[0]) != 4:
		raise ValueError("Expected unitary of shape (4,4).")
	special = to_nested_list(to_special(m))
	yy = [[0j, 0j, 0j, -1 + 0j], [0j, 0j, 1 + 0j, 0j], [0j, 1 + 0j, 0j, 0j], [-1 + 0j, 0j, 0j, 0j]]
	g = _mm(_mm(_mm(special, yy), _transpose(special)), yy)
	a3 = -sum(g[i][i] for i in range(4))
	if abs(a3 - 4) < atol or abs(a3 + 4) < atol:
		return 0
	a2 = (a3 * a3 - sum(_mm(g, g)[i][i] for i in range(4))) / 2
	if abs(a3) < atol and abs(a2 - 2) < atol:
		return 1
	if abs(a3.imag) < atol:
		return 2
	return 3


def extract_right_diag(u: Sequence[Sequence[complex]], atol: float = 1e-8) -> Matrix:
	"""A diagonal unitary ``D`` so that ``u @ D`` needs only two CNOTs.

	This is ``cirq.extract_right_diag`` (Proposition V.2 of
	quant-ph/0308033): it only makes sense for a unitary that needs three
	CNOTs, which is what the two-qubit diagonal-plus-CZ decomposition uses.
	"""
	m = to_nested_list(u)
	if len(m) != 4 or any(len(row) != 4 for row in m):
		raise ValueError(f"Expected unitary of shape (4,4), got {len(m)}x{len(m[0]) if m else 0}.")
	yy: Matrix = [
		[0j, 0j, 0j, -1 + 0j],
		[0j, 0j, 1 + 0j, 0j],
		[0j, 1 + 0j, 0j, 0j],
		[-1 + 0j, 0j, 0j, 0j],
	]
	special = to_nested_list(to_special(m))
	transposed = _transpose(special)
	gamma = _mm(_mm(_mm(transposed, yy), special), yy)
	diagonal = [gamma[i][i] for i in range(4)]
	k = (diagonal[0] + diagonal[3] - diagonal[1] - diagonal[2]).real
	psi = math.atan2(sum(diagonal).imag, k)
	f = cmath.exp(1j * psi)
	del atol
	return [[1 + 0j, 0j, 0j, 0j], [0j, f, 0j, 0j], [0j, 0j, f, 0j], [0j, 0j, 0j, 1 + 0j]]


def slice_for_qubits_equal_to(
	target_qubit_axes: Sequence[int],
	little_endian_qureg_value: int = 0,
	*,
	big_endian_qureg_value: int = 0,
	num_qubits: int | None = None,
	qid_shape: tuple[int, ...] | None = None,
) -> tuple[Any, ...]:
	"""Return an index selecting a subset of a state tensor."""
	if qid_shape is not None or num_qubits is not None:
		shape_specified: tuple[int, ...] = qid_shape if qid_shape is not None else (2,) * (num_qubits or 0)
		if num_qubits is None:
			num_qubits = len(shape_specified)
		if num_qubits != len(shape_specified):
			raise ValueError("len(qid_shape) != num_qubits")
		qid_shape = shape_specified
	if little_endian_qureg_value and big_endian_qureg_value:
		raise ValueError("Specify exactly one of little_endian_qureg_value or big_endian_qureg_value.")
	out_size_specified = num_qubits is not None
	out_size = num_qubits if out_size_specified else max(target_qubit_axes, default=-1) + 1
	assert out_size is not None
	result: list[Any] = [slice(None)] * out_size
	if not out_size_specified:
		result.append(Ellipsis)
	shape = qid_shape or (2,) * out_size  # type: ignore[assignment]
	target_shape = [shape[i] for i in target_qubit_axes]
	if big_endian_qureg_value:
		digits = _int_to_digits(big_endian_qureg_value, target_shape)
	else:
		digits = _int_to_digits(little_endian_qureg_value, target_shape[::-1])[::-1]
	for axis, digit in zip(target_qubit_axes, digits, strict=True):
		result[axis] = digit
	return tuple(result)


def _int_to_digits(value: int, bases: Sequence[int]) -> list[int]:
	digits = [0] * len(bases)
	for i in range(len(bases) - 1, -1, -1):
		digits[i] = value % bases[i]
		value //= bases[i]
	return digits


__all__ = [
	"CONTROL_TAG",
	"PAULI_BASIS",
	"AxisAngleDecomposition",
	"KakDecomposition",
	"all_near_zero",
	"all_near_zero_mod",
	"allclose_up_to_global_phase",
	"apply_matrix_to_slices",
	"axis_angle",
	"bidiagonalize_real_matrix_pair_with_symmetric_products",
	"bidiagonalize_unitary_with_special_orthogonals",
	"block_diag",
	"density_matrix_kronecker_product",
	"diagonalize_real_symmetric_and_sorted_diagonal_matrices",
	"diagonalize_real_symmetric_matrix",
	"dot",
	"expand_matrix_in_orthogonal_basis",
	"extract_right_diag",
	"hilbert_schmidt_inner_product",
	"is_cptp",
	"is_diagonal",
	"is_hermitian",
	"is_normal",
	"is_orthogonal",
	"is_special_orthogonal",
	"is_special_unitary",
	"is_unitary",
	"kak_canonicalize_vector",
	"kak_decomposition",
	"kak_vector",
	"kron",
	"kron_bases",
	"kron_factor_4x4_to_2x2s",
	"kron_with_controls",
	"map_eigenvalues",
	"match_global_phase",
	"matrix_commutes",
	"matrix_from_basis_coefficients",
	"num_cnots_required",
	"partial_trace_of_state_vector_as_mixture",
	"pow_pauli_combination",
	"reflection_matrix_pow",
	"slice_for_qubits_equal_to",
	"so4_to_magic_su2s",
	"state_vector_kronecker_product",
	"sub_state_vector",
	"targeted_conjugate_about",
	"targeted_left_multiply",
	"to_special",
	"unitary_eig",
]
