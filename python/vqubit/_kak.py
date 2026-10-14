"""KAK (Cartan) decomposition and two-qubit synthesis.

Any two-qubit unitary can be written, up to a global phase, as a single-qubit
gate on each qubit, an interaction ``exp(i(x XX + y YY + z ZZ))``, and another
single-qubit gate on each qubit. That is the KAK decomposition, and it is the
key to compiling an arbitrary two-qubit gate into single-qubit gates plus an
entangler.

The decomposition is computed in the *magic basis*, where local gates become
real orthogonal matrices and the interaction becomes diagonal. Because the
matrix ``U^T U`` there is a symmetric unitary, its real and imaginary parts are
commuting real symmetric matrices, which a Jacobi eigenvalue routine can
diagonalize together -- so no complex eigensolver is needed.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Sequence

from ._gates import CNOTGate, CZGate, H, Rx, Rz, S
from ._linalg import Matrix, det, matmul, real_symmetric_eig, transpose
from ._ops import Operation
from ._qubits import Qid

_SQRT_HALF = math.sqrt(0.5)

# Cirq's magic basis. Local gates become real orthogonal in it.
_MAGIC: Matrix = [
	[_SQRT_HALF + 0j, 0j, 0j, 1j * _SQRT_HALF],
	[0j, 1j * _SQRT_HALF, _SQRT_HALF + 0j, 0j],
	[0j, 1j * _SQRT_HALF, -_SQRT_HALF + 0j, 0j],
	[_SQRT_HALF + 0j, 0j, 0j, -1j * _SQRT_HALF],
]
_MAGIC_DAG = [[_MAGIC[j][i].conjugate() for j in range(4)] for i in range(4)]

_SIGMA: dict[str, Matrix] = {
	"X": [[0j, 1 + 0j], [1 + 0j, 0j]],
	"Y": [[0j, -1j], [1j, 0j]],
	"Z": [[1 + 0j, 0j], [0j, -1 + 0j]],
}


def _kron4(a: Matrix, b: Matrix) -> Matrix:
	return [[a[i // 2][j // 2] * b[i % 2][j % 2] for j in range(4)] for i in range(4)]


# Diagonal matrices (in the magic basis) of XX, YY, ZZ, used to read the
# interaction angles back out of the diagonal factor.
_CARTAN = {
	key: [
		[matmul(matmul(_MAGIC_DAG, _kron4(_SIGMA[key], _SIGMA[key])), _MAGIC)[i][j] for j in range(4)] for i in range(4)
	]
	for key in ("X", "Y", "Z")
}


def _simultaneous_real_diag(
	real_part: Sequence[Sequence[float]], imag_part: Sequence[Sequence[float]]
) -> tuple[list[complex], list[list[float]]]:
	"""Jointly diagonalize two commuting real symmetric matrices.

	Returns the complex eigenvalues and a real orthogonal matrix ``O`` with
	``(real_part + i imag_part) = O diag(values) O^T``.
	"""
	values_a, vectors = real_symmetric_eig(real_part)
	n = len(values_a)
	# Group equal eigenvalues of the real part.
	order = sorted(range(n), key=lambda i: values_a[i])
	groups: list[list[int]] = []
	for index in order:
		if groups and abs(values_a[index] - values_a[groups[-1][0]]) < 1e-9:
			groups[-1].append(index)
		else:
			groups.append([index])

	for group in groups:
		if len(group) < 2:
			continue
		# Diagonalize the imaginary part inside this degenerate subspace.
		sub = [
			[
				sum(vectors[r][group[i]] * imag_part[r][s] * vectors[s][group[j]] for r in range(n) for s in range(n))
				for j in range(len(group))
			]
			for i in range(len(group))
		]
		sub_values, sub_vectors = real_symmetric_eig(sub)
		_ = sub_values
		rotated = [
			[sum(vectors[r][group[k]] * sub_vectors[k][j].real for k in range(len(group))) for j in range(len(group))]
			for r in range(n)
		]
		for j, column in enumerate(group):
			for r in range(n):
				vectors[r][column] = rotated[r][j]

	complex_values: list[complex] = []
	for j in range(n):
		real_v = sum(vectors[r][j].real * real_part[r][s] * vectors[s][j].real for r in range(n) for s in range(n))
		imag_v = sum(vectors[r][j].real * imag_part[r][s] * vectors[s][j].real for r in range(n) for s in range(n))
		complex_values.append(complex(real_v, imag_v))
	return complex_values, vectors


def _kron_factor(matrix: Matrix) -> tuple[complex, Matrix, Matrix]:
	"""Factor a 4x4 matrix that is ``g * kron(a, b)`` into ``g``, ``a`` and ``b``."""
	best = max(((i, j) for i in range(4) for j in range(4)), key=lambda t: abs(matrix[t[0]][t[1]]))
	a, b = best
	f1 = [[0j, 0j], [0j, 0j]]
	f2 = [[0j, 0j], [0j, 0j]]
	for i in range(2):
		for j in range(2):
			f1[(a >> 1) ^ i][(b >> 1) ^ j] = matrix[a ^ (i << 1)][b ^ (j << 1)]
			f2[(a & 1) ^ i][(b & 1) ^ j] = matrix[a ^ i][b ^ j]
	for factor in (f1, f2):
		determinant = det(factor)
		if abs(determinant) > 1e-12:
			root = cmath.sqrt(determinant)
			for i in range(2):
				for j in range(2):
					factor[i][j] /= root
	global_phase = matrix[a][b] / (f1[a >> 1][b >> 1] * f2[a & 1][b & 1])
	if global_phase.real < 0:
		f1 = [[-z for z in row] for row in f1]
		global_phase = -global_phase
	return global_phase, f1, f2


def kak_decomposition(
	matrix: Sequence[Sequence[complex]],
	*,
	atol: float = 1e-8,
) -> tuple[complex, tuple[float, float, float], tuple[Matrix, Matrix], tuple[Matrix, Matrix]]:
	"""Return ``(phase, (x, y, z), (a1, a2), (b1, b2))`` for a two-qubit unitary.

	so that ``U = phase * (a1 (x) a2) exp(i(x XX + y YY + z ZZ)) (b1 (x) b2)``.
	"""
	unitary = [[complex(z) for z in row] for row in matrix]
	# The Cartan form has unit determinant, so pull the global phase of a U(4)
	# element into `unitary` first (a phase on the whole gate is irrelevant).
	phase = det(unitary) ** (-0.25)
	unitary = [[phase * z for z in row] for row in unitary]
	magic_conv = matmul(matmul(_MAGIC_DAG, unitary), _MAGIC)
	symmetric = matmul(transpose(magic_conv), magic_conv)

	real_part = [[symmetric[i][j].real for j in range(4)] for i in range(4)]
	imag_part = [[symmetric[i][j].imag for j in range(4)] for i in range(4)]
	values, vectors = _simultaneous_real_diag(real_part, imag_part)

	roots = [cmath.sqrt(value) for value in values]
	# The interaction factor must have unit determinant (it is generated by
	# XX, YY and ZZ, whose traces vanish). A single sign flip fixes the
	# "swap class" without changing the reconstructed unitary.
	if (roots[0] * roots[1] * roots[2] * roots[3]).real < 0:
		roots[0] = -roots[0]
	inverse_roots = [[1.0 / roots[i] if i == j else (0 + 0j) for j in range(4)] for i in range(4)]
	left = matmul(matmul(magic_conv, vectors), inverse_roots)
	right = transpose(vectors)

	left = [[complex(z.real, 0.0) for z in row] for row in left]
	right = [[complex(z.real, 0.0) for z in row] for row in right]

	after = matmul(matmul(_MAGIC, left), _MAGIC_DAG)
	before = matmul(matmul(_MAGIC, right), _MAGIC_DAG)

	g1, a1, a2 = _kron_factor(after)
	g2, b1, b2 = _kron_factor(before)

	# The interaction diagonal in the magic basis gives the three angles.
	x, y, z = _angles_from_diagonal(roots)
	return g1 * g2, (x, y, z), (a1, a2), (b1, b2)


def _angles_from_diagonal(roots: Sequence[complex]) -> tuple[float, float, float]:
	"""Recover the interaction angles from the magic-basis diagonal factor.

	The phases are only known modulo ``2*pi``; the four Cartan rows sum to zero,
	which pins the total shift and makes the 4x3 system consistent.
	"""
	phases = [cmath.phase(root) for root in roots]
	shift = round(-sum(phases) / (2 * math.pi))
	rows = []
	for k in range(4):
		rows.append(
			[
				_CARTAN["X"][k][k].real,
				_CARTAN["Y"][k][k].real,
				_CARTAN["Z"][k][k].real,
			]
		)
	right = list(phases)
	right[3] += 2 * math.pi * shift
	return _solve3(rows[:3], right[:3])


def _solve3(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> tuple[float, float, float]:
	a = [[float(z) for z in row] + [float(rhs[i])] for i, row in enumerate(matrix)]
	for col in range(3):
		pivot = max(range(col, 3), key=lambda r: abs(a[r][col]))
		a[col], a[pivot] = a[pivot], a[col]
		inv = 1.0 / a[col][col]
		a[col] = [z * inv for z in a[col]]
		for r in range(3):
			if r == col:
				continue
			factor = a[r][col]
			a[r] = [a[r][c] - factor * a[col][c] for c in range(4)]
	return a[0][3], a[1][3], a[2][3]


def _single_qubit_gates(matrix: Matrix, qubit: Qid) -> list[Operation]:
	"""Decompose a 2x2 unitary into ``Rz`` and ``Ry`` operations on ``qubit``."""
	from ._optimizers import single_qubit_matrix_to_gates

	return [gate.on(qubit) for gate in single_qubit_matrix_to_gates(matrix)]


def _xx_interaction(theta: float, first: Qid, second: Qid) -> list[Operation]:
	"""Return operations for ``exp(i theta XX)``.

	Conjugating an ``X`` on the control of a CNOT produces ``XX``, so the
	rotation goes on the control qubit.
	"""
	return [
		CNOTGate().on(first, second),
		Rx(-2 * theta).on(first),
		CNOTGate().on(first, second),
	]


def _yy_interaction(theta: float, first: Qid, second: Qid) -> list[Operation]:
	return [
		(S**-1).on(first),
		(S**-1).on(second),
		*_xx_interaction(theta, first, second),
		S.on(first),
		S.on(second),
	]


def _zz_interaction(theta: float, first: Qid, second: Qid) -> list[Operation]:
	"""Return operations for ``exp(i theta ZZ)`` (rotation on the CNOT target)."""
	return [
		CNOTGate().on(first, second),
		Rz(-2 * theta).on(second),
		CNOTGate().on(first, second),
	]


def two_qubit_matrix_to_operations(
	matrix: Sequence[Sequence[complex]],
	*,
	qubits: tuple[Qid, Qid],
	atol: float = 1e-8,
) -> list[Operation]:
	"""Decompose an arbitrary two-qubit unitary into single-qubit gates and CNOTs."""
	_, (x, y, z), (a1, a2), (b1, b2) = kak_decomposition(matrix, atol=atol)
	first, second = qubits
	operations: list[Operation] = []
	operations.extend(_single_qubit_gates(b1, first))
	operations.extend(_single_qubit_gates(b2, second))
	operations.extend(_xx_interaction(x, first, second))
	operations.extend(_yy_interaction(y, first, second))
	operations.extend(_zz_interaction(z, first, second))
	operations.extend(_single_qubit_gates(a1, first))
	operations.extend(_single_qubit_gates(a2, second))
	return operations


def _xx_interaction_cz(theta: float, first: Qid, second: Qid) -> list[Operation]:
	# CNOT = (H on target) CZ (H on target).
	return [
		H.on(second),
		CZGate().on(first, second),
		H.on(second),
		Rx(-2 * theta).on(first),
		H.on(second),
		CZGate().on(first, second),
		H.on(second),
	]


def _zz_interaction_cz(theta: float, first: Qid, second: Qid) -> list[Operation]:
	return [
		H.on(second),
		CZGate().on(first, second),
		H.on(second),
		Rz(-2 * theta).on(second),
		H.on(second),
		CZGate().on(first, second),
		H.on(second),
	]


def _yy_interaction_cz(theta: float, first: Qid, second: Qid) -> list[Operation]:
	return [
		(S**-1).on(first),
		(S**-1).on(second),
		*_xx_interaction_cz(theta, first, second),
		S.on(first),
		S.on(second),
	]


def two_qubit_matrix_to_cz_operations(
	matrix: Sequence[Sequence[complex]],
	*,
	qubits: tuple[Qid, Qid],
	atol: float = 1e-8,
) -> list[Operation]:
	"""Decompose an arbitrary two-qubit unitary into single-qubit gates and CZs."""
	_, (x, y, z), (a1, a2), (b1, b2) = kak_decomposition(matrix, atol=atol)
	first, second = qubits
	operations: list[Operation] = []
	operations.extend(_single_qubit_gates(b1, first))
	operations.extend(_single_qubit_gates(b2, second))
	operations.extend(_xx_interaction_cz(x, first, second))
	operations.extend(_yy_interaction_cz(y, first, second))
	operations.extend(_zz_interaction_cz(z, first, second))
	operations.extend(_single_qubit_gates(a1, first))
	operations.extend(_single_qubit_gates(a2, second))
	return operations


__all__ = [
	"CZGate",
	"kak_decomposition",
	"two_qubit_matrix_to_cz_operations",
	"two_qubit_matrix_to_operations",
]
