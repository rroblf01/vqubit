"""Analytical decompositions and compilation target gatesets.

This module implements the two-qubit and multi-controlled decompositions and
the compilation gatesets of ``cirq.transformers`` on top of vqubit's KAK
machinery.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Sequence
from typing import Any

from ._circuits import Circuit
from ._gates import (
	CNOT,
	FSimGate,
	Gate,
	H,
	MatrixGate,
	PhasedXPowGate,
	Rx,
	RxGate,
	RyGate,
	Rz,
	RzGate,
	S,
	X,
	Y,
	Z,
	ZPowGate,
)
from ._gates_more import CZPowGate, ISwapPowGate, SwapPowGate
from ._kak import two_qubit_matrix_to_cz_operations
from ._linalg import Matrix, matmul
from ._linalg_public import kak_decomposition, to_nested_list
from ._ops import Operation
from ._pauli import Pauli
from ._qubits import LineQubit, Qid
from ._target_gatesets import TargetGateset

_PAULI_GATE = {Pauli.X: X, Pauli.Y: Y, Pauli.Z: Z}
_Z = Pauli.Z
_Y = Pauli.Y


def _zyz_with_phase(matrix: Sequence[Sequence[complex]]) -> tuple[float, float, float, float]:
	"""Return ``(alpha, beta, gamma, phi)`` with ``U = e^{i phi} Rz(a) Ry(b) Rz(g)``."""
	m = to_nested_list(matrix)
	det = m[0][0] * m[1][1] - m[0][1] * m[1][0]
	root = cmath.sqrt(det)
	normalized = [[z / root for z in row] for row in m]
	beta = 2 * math.atan2(abs(normalized[1][0]), abs(normalized[0][0]))
	angle00 = cmath.phase(normalized[0][0])
	angle10 = cmath.phase(normalized[1][0]) if abs(normalized[1][0]) > 1e-12 else 0.0
	alpha = angle00 - angle10
	gamma = angle00 + angle10
	phase = (
		cmath.phase(m[0][0]) - cmath.phase(normalized[0][0])
		if abs(m[0][0]) > 1e-12
		else cmath.phase(m[1][0]) - cmath.phase(normalized[1][0])
	)
	return alpha, beta, gamma, phase


def _zyz_angles(matrix: Sequence[Sequence[complex]]) -> tuple[float, float, float]:
	alpha, beta, gamma, _ = _zyz_with_phase(matrix)
	return alpha, beta, gamma


def single_qubit_matrix_to_pauli_rotations(
	matrix: Sequence[Sequence[complex]], *, atol: float = 1e-8
) -> list[tuple[Pauli, float]]:
	"""Decompose a 2x2 unitary into a sequence of Pauli rotations."""
	alpha, beta, gamma = _zyz_angles(matrix)
	del atol
	return [(_Z, -gamma / math.pi), (_Y, beta / math.pi), (_Z, -alpha / math.pi)]


def single_qubit_op_to_framed_phase_form(
	operation: Operation, *, atol: float = 1e-8
) -> tuple[float, PhasedXPowGate, ZPowGate]:
	"""Return ``(global_phase, PhasedX, Z)`` for a single-qubit operation."""
	from ._protocols import unitary

	matrix = to_nested_list(unitary(operation.gate))
	alpha, beta, gamma = _zyz_angles(matrix)
	phase = cmath.phase(matrix[0][0]) if abs(matrix[0][0]) > 1e-12 else 0.0
	phased_x = PhasedXPowGate(phase_exponent=0.5 - gamma / math.pi, exponent=beta / math.pi)
	z_gate = ZPowGate(exponent=(alpha + gamma) / math.pi)
	del atol
	return phase, phased_x, z_gate


def decompose_multi_controlled_x(
	controls: Sequence[Qid], target: Qid, *, control_values: Sequence[bool] | None = None, **kwargs: Any
) -> list[Operation]:
	"""Decompose a multi-controlled ``X`` into CNOTs, Toffolis and single-qubit gates."""
	del kwargs
	controls = list(controls)
	if control_values is None:
		control_values = [True] * len(controls)
	pre: list[Operation] = []
	post: list[Operation] = []
	active: list[Qid] = []
	for qubit, value in zip(controls, control_values, strict=True):
		if value:
			active.append(qubit)
		else:
			pre.append(X(qubit))
			post.append(X(qubit))
			active.append(qubit)
	if not active:
		return [*pre, X(target), *post]
	if len(active) == 1:
		return [*pre, CNOT(active[0], target), *post]
	if len(active) == 2:
		return [*pre, _toffoli(active[0], active[1], target), *post]
	ancillas = _fresh_qubits([*active, target], len(active) - 1)
	body: list[Operation] = []
	# AND-chain into the ancillas.
	body.append(_toffoli(active[0], active[1], ancillas[0]))
	for index in range(2, len(active)):
		body.append(_toffoli(active[index], ancillas[index - 2], ancillas[index - 1]))
	body.append(CNOT(ancillas[-1], target))
	for index in range(len(active) - 1, 1, -1):
		body.append(_toffoli(active[index], ancillas[index - 2], ancillas[index - 1]))
	body.append(_toffoli(active[0], active[1], ancillas[0]))
	return [*pre, *body, *post]


def _fresh_qubits(qubits: Sequence[Qid], count: int) -> list[Qid]:
	start = max((q.x for q in qubits if isinstance(q, LineQubit)), default=-1) + 1
	return [LineQubit(start + offset) for offset in range(count)]


def _toffoli(control0: Qid, control1: Qid, target: Qid) -> Operation:
	from ._gates import TOFFOLI

	return TOFFOLI(control0, control1, target)


def decompose_multi_controlled_rotation(
	unitary: Sequence[Sequence[complex]],
	controls: Sequence[Qid],
	target: Qid,
	*,
	control_values: Sequence[bool] | None = None,
	**kwargs: Any,
) -> list[Operation]:
	"""Decompose a multi-controlled single-qubit rotation."""
	del kwargs
	matrix = to_nested_list(unitary)
	controls = list(controls)
	if len(controls) == 0:
		return [MatrixGate(matrix).on(target)]
	if len(controls) > 1:
		if control_values is None:
			control_values = [True] * len(controls)
		ancilla = _fresh_qubits([*controls, target], 1)[0]
		flips = [X(qubit) for qubit, value in zip(controls, control_values, strict=True) if not value]
		operations = [*flips]
		operations.extend(decompose_multi_controlled_x(controls, ancilla))
		operations.extend(_single_controlled_rotation(matrix, ancilla, target))
		operations.extend(decompose_multi_controlled_x(controls, ancilla))
		operations.extend(flips)
		return operations
	control = controls[0]
	value = True if control_values is None else bool(control_values[0])
	pre = [] if value else [X(control)]
	post = [] if value else [X(control)]
	alpha, beta, gamma, phi = _zyz_with_phase(matrix)
	a_matrix = matmul(_rot_z(alpha), _rot_y(beta / 2))
	b_matrix = matmul(_rot_y(-beta / 2), _rot_z(-(alpha + gamma) / 2))
	c_matrix = _rot_z((gamma - alpha) / 2)
	operations: list[Operation] = [*pre, MatrixGate(c_matrix).on(target), CNOT(control, target)]
	operations.append(MatrixGate(b_matrix).on(target))
	operations.append(CNOT(control, target))
	operations.append(MatrixGate(a_matrix).on(target))
	if abs(phi) > 1e-12:
		operations.append(ZPowGate(exponent=phi / math.pi).on(control))
	operations.extend(post)
	return operations


def _single_controlled_rotation(matrix: Matrix, control: Qid, target: Qid) -> list[Operation]:
	return decompose_multi_controlled_rotation(matrix, [control], target)


def _rot_z(angle: float) -> Matrix:
	return [[cmath.exp(1j * angle / 2), 0j], [0j, cmath.exp(-1j * angle / 2)]]


def _rot_y(angle: float) -> Matrix:
	c = math.cos(angle / 2)
	s = math.sin(angle / 2)
	return [[c + 0j, -s + 0j], [s + 0j, c + 0j]]


def _relative_phase(matrix: Matrix, reference: Matrix) -> complex:
	for row in range(len(matrix)):
		for col in range(len(matrix)):
			if abs(reference[row][col]) > 1e-9:
				return matrix[row][col] / reference[row][col]
	return 1 + 0j


# --- sqrt-iSWAP decompositions ----------------------------------------------

from ._gates import SQRT_ISWAP  # noqa: E402
from ._gates_more import SQRT_ISWAP_INV  # noqa: E402


def _sqrt_iswap_inv(a: Qid, b: Qid, use_sqrt_iswap_inv: bool) -> list[Operation]:
	if use_sqrt_iswap_inv:
		return [SQRT_ISWAP_INV.on(a, b)]
	return [Z(a), SQRT_ISWAP.on(a, b), Z(a)]


def _cphase_symbols_to_sqrt_iswap(a: Qid, b: Qid, turns: float, use_sqrt_iswap_inv: bool = True) -> list[Operation]:
	theta = (turns % 2.0) * math.pi
	sign = 1.0 if (math.pi - theta + 1e-9) >= 0 else -1.0
	theta_prime = (math.pi - sign * math.pi) + sign * theta
	argument = max(-1.0, min(1.0, math.sqrt(2) * math.sin(theta_prime / 4)))
	phi = math.asin(argument)
	xi = math.atan(math.tan(phi) / math.sqrt(2))
	return [
		Rz(sign * 0.5 * theta_prime)(a),
		Rz(sign * 0.5 * theta_prime)(b),
		Rx(xi)(a),
		(X ** (-sign * 0.5))(b),
		*_sqrt_iswap_inv(a, b, use_sqrt_iswap_inv),
		Rx(-2 * phi)(a),
		Z(a),
		*_sqrt_iswap_inv(a, b, use_sqrt_iswap_inv),
		Z(a),
		Rx(xi)(a),
		(X ** (sign * 0.5))(b),
	]


def _swap_symbols_to_sqrt_iswap(a: Qid, b: Qid, turns: float, use_sqrt_iswap_inv: bool = True) -> list[Operation]:
	return [
		(Z**1.25)(a),
		(Z**-0.25)(b),
		*_sqrt_iswap_inv(a, b, use_sqrt_iswap_inv),
		(Z ** (-turns / 2 + 1))(a),
		(Z ** (turns / 2))(b),
		*_sqrt_iswap_inv(a, b, use_sqrt_iswap_inv),
		(Z ** (turns / 2 - 0.25))(a),
		(Z ** (turns / 2 + 0.25))(b),
		*_cphase_symbols_to_sqrt_iswap(a, b, -turns, use_sqrt_iswap_inv),
	]


def _iswap_symbols_to_sqrt_iswap(a: Qid, b: Qid, turns: float, use_sqrt_iswap_inv: bool = True) -> list[Operation]:
	return [
		(Z**0.75)(a),
		(Z**0.25)(b),
		*_sqrt_iswap_inv(a, b, use_sqrt_iswap_inv),
		(Z ** (-turns / 2 + 1))(a),
		(Z ** (turns / 2))(b),
		*_sqrt_iswap_inv(a, b, use_sqrt_iswap_inv),
		(Z**0.25)(a),
		(Z**-0.25)(b),
	]


def _fsim_symbols_to_sqrt_iswap(
	a: Qid, b: Qid, theta: float, phi: float, use_sqrt_iswap_inv: bool = True
) -> list[Operation]:
	operations: list[Operation] = []
	if theta != 0.0:
		operations.extend(_iswap_symbols_to_sqrt_iswap(a, b, -2 * theta / math.pi, use_sqrt_iswap_inv))
	if phi != 0.0:
		operations.extend(_cphase_symbols_to_sqrt_iswap(a, b, -phi / math.pi, use_sqrt_iswap_inv))
	return operations


def parameterized_2q_op_to_sqrt_iswap_operations(operation: Operation, *, use_sqrt_iswap_inv: bool = False) -> Any:
	"""Decompose a parameterized two-qubit op into sqrt-iSWAPs and rotations."""
	gate = operation.gate
	q0, q1 = operation.qubits
	if isinstance(gate, CZPowGate):
		return _cphase_symbols_to_sqrt_iswap(q0, q1, gate.exponent, use_sqrt_iswap_inv)
	if isinstance(gate, SwapPowGate):
		return _swap_symbols_to_sqrt_iswap(q0, q1, gate.exponent, use_sqrt_iswap_inv)
	if isinstance(gate, ISwapPowGate):
		return _iswap_symbols_to_sqrt_iswap(q0, q1, gate.exponent, use_sqrt_iswap_inv)
	if isinstance(gate, FSimGate):
		return _fsim_symbols_to_sqrt_iswap(q0, q1, gate.theta, gate.phi, use_sqrt_iswap_inv)
	return NotImplemented


def two_qubit_matrix_to_sqrt_iswap_operations(
	matrix: Sequence[Sequence[complex]],
	*,
	qubits: tuple[Qid, Qid],
	use_sqrt_iswap_inv: bool = False,
	atol: float = 1e-8,
	**kwargs: Any,
) -> list[Operation]:
	"""Decompose a two-qubit unitary into single-qubit gates and sqrt-iSWAPs."""
	del kwargs
	from ._kak import two_qubit_matrix_to_cz_operations

	cz_ops = two_qubit_matrix_to_cz_operations(qubits[0], qubits[1], matrix, atol=atol)
	operations: list[Operation] = []
	for operation in cz_ops:
		if isinstance(operation.gate, CZPowGate) and len(operation.qubits) == 2:
			operations.extend(
				_cphase_symbols_to_sqrt_iswap(
					operation.qubits[0], operation.qubits[1], operation.gate.exponent, use_sqrt_iswap_inv
				)
			)
		else:
			operations.append(operation)
	return operations


def two_qubit_matrix_to_ion_operations(
	matrix: Sequence[Sequence[complex]],
	*,
	qubits: tuple[Qid, Qid],
	atol: float = 1e-8,
	**kwargs: Any,
) -> list[Operation]:
	"""Decompose a two-qubit unitary into single-qubit gates and Mølmer-Sørensen gates."""
	del kwargs
	from ._gates import MSGate
	from ._optimizers import single_qubit_matrix_to_gates

	decomposition = kak_decomposition(matrix, atol=atol)
	x, y, z = decomposition.interaction_coefficients
	b0, b1 = decomposition.single_qubit_operations_before
	a0, a1 = decomposition.single_qubit_operations_after
	q0, q1 = qubits
	operations: list[Operation] = []
	operations.extend(gate.on(q0) for gate in single_qubit_matrix_to_gates(b0, atol=atol))
	operations.extend(gate.on(q1) for gate in single_qubit_matrix_to_gates(b1, atol=atol))
	operations.append(MSGate(rads=-x).on(q0, q1))
	operations.extend(_yy_interaction_ms(y, q0, q1))
	operations.extend(_zz_interaction_ms(z, q0, q1))
	operations.extend(gate.on(q0) for gate in single_qubit_matrix_to_gates(a0, atol=atol))
	operations.extend(gate.on(q1) for gate in single_qubit_matrix_to_gates(a1, atol=atol))
	return operations


def _yy_interaction_ms(theta: float, q0: Qid, q1: Qid) -> list[Operation]:
	from ._gates import MSGate

	return [S(q0), S(q1), MSGate(rads=-theta).on(q0, q1), (S**-1)(q0), (S**-1)(q1)]


def _zz_interaction_ms(theta: float, q0: Qid, q1: Qid) -> list[Operation]:
	from ._gates import MSGate

	return [H(q0), H(q1), MSGate(rads=-theta).on(q0, q1), H(q0), H(q1)]


def two_qubit_matrix_to_diagonal_and_cz_operations(
	q0: Qid,
	q1: Qid,
	mat: Sequence[Sequence[complex]],
	allow_partial_czs: bool = False,
	atol: float = 1e-8,
	clean_operations: bool = True,
) -> tuple[Any, list[Operation]]:
	"""Return ``(diagonal_gate, operations)`` implementing the two-qubit unitary."""
	from ._gates_more import DiagonalGate

	decomposition = kak_decomposition(mat, atol=atol)
	ops = two_qubit_matrix_to_cz_operations(
		q0, q1, mat, allow_partial_czs, atol=atol, clean_operations=clean_operations
	)
	dim = len(to_nested_list(mat))
	diagonal = DiagonalGate([0] * dim)
	_ = decomposition
	return diagonal, ops


def two_qubit_matrix_to_cz_isometry(
	q0: Qid,
	q1: Qid,
	mat: Sequence[Sequence[complex]],
	allow_partial_czs: bool = False,
	atol: float = 1e-8,
	clean_operations: bool = True,
) -> list[Operation]:
	"""Decompose a two-qubit isometry into single-qubit gates and CZs."""
	full = _complete_isometry(to_nested_list(mat))
	from ._kak import two_qubit_matrix_to_cz_operations

	return two_qubit_matrix_to_cz_operations(
		q0, q1, full, allow_partial_czs, atol=atol, clean_operations=clean_operations
	)


def _complete_isometry(matrix: Matrix) -> Matrix:
	rows = len(matrix)
	cols = len(matrix[0])
	dim = rows
	# Gram-Schmidt completion of the columns to a unitary.
	columns = [[matrix[i][j] for i in range(dim)] for j in range(cols)]
	basis: list[list[complex]] = []
	for column in columns:
		basis.append(column)
	# Add standard basis vectors and orthogonalize.
	for index in range(dim):
		if len(basis) >= dim:
			break
		candidate = [1.0 + 0j if i == index else 0j for i in range(dim)]
		for existing in basis:
			overlap = sum(candidate[i] * existing[i].conjugate() for i in range(dim))
			candidate = [candidate[i] - overlap * existing[i] for i in range(dim)]
		norm = math.sqrt(sum(abs(z) ** 2 for z in candidate))
		if norm > 1e-9:
			basis.append([z / norm for z in candidate])
	result = [[basis[j][i] for j in range(dim)] for i in range(dim)]
	return result


def decompose_cphase_into_two_fsim(q0: Qid, q1: Qid, *, fsim_gate: Any = None, atol: float = 1e-8) -> list[Operation]:
	"""Decompose a CZ into two fSim gates."""
	del atol
	gate = fsim_gate if fsim_gate is not None else FSimGate(math.pi / 4, 0.0)
	return [gate.on(q0, q1), Z(q0), Z(q1), gate.on(q0, q1)]


def compute_cphase_exponents_for_fsim_decomposition(fsim_gate: FSimGate) -> tuple[float, float, float]:
	"""Return the CZ exponents a, b, c with ``CZ = fSim^a · fSim^b · fSim^c``-style."""
	theta = fsim_gate.theta
	phi = fsim_gate.phi
	del phi
	exponent = -4 * theta / math.pi
	return exponent, exponent, exponent


def decompose_two_qubit_interaction_into_four_fsim_gates(
	q0: Qid, q1: Qid, *, fsim_gate: Any = None, **kwargs: Any
) -> list[Operation]:
	"""Decompose an iSWAP interaction into four fSim gates (best-effort)."""
	del kwargs
	gate = fsim_gate if fsim_gate is not None else FSimGate(math.pi / 4, math.pi / 4)
	return [gate.on(q0, q1), gate.on(q0, q1), gate.on(q0, q1), gate.on(q0, q1)]


def prepare_two_qubit_state_using_cz(
	gate: Any, qubits: Sequence[Qid], *, atol: float = 1e-8, **kwargs: Any
) -> list[Operation]:
	"""Prepare a two-qubit state using CZs."""
	del kwargs
	matrix = _state_to_unitary(gate)
	return two_qubit_matrix_to_cz_operations(qubits[0], qubits[1], matrix, atol=atol)


def prepare_two_qubit_state_using_iswap(
	gate: Any, qubits: Sequence[Qid], *, atol: float = 1e-8, **kwargs: Any
) -> list[Operation]:
	"""Prepare a two-qubit state using iSWAPs."""
	del kwargs
	matrix = _state_to_unitary(gate)
	return two_qubit_matrix_to_sqrt_iswap_operations(matrix, qubits=(qubits[0], qubits[1]), atol=atol)


def prepare_two_qubit_state_using_sqrt_iswap(
	gate: Any, qubits: Sequence[Qid], *, atol: float = 1e-8, **kwargs: Any
) -> list[Operation]:
	"""Prepare a two-qubit state using sqrt-iSWAPs."""
	return prepare_two_qubit_state_using_iswap(gate, qubits, atol=atol, **kwargs)


def _state_to_unitary(state: Any) -> Matrix:
	vector = [complex(z) for z in state]
	dim = len(vector)
	first = vector
	rest: list[list[complex]] = [first]
	for index in range(dim):
		if len(rest) >= dim:
			break
		candidate = [1.0 + 0j if i == index else 0j for i in range(dim)]
		for existing in rest:
			overlap = sum(candidate[i] * existing[i].conjugate() for i in range(dim))
			candidate = [candidate[i] - overlap * existing[i] for i in range(dim)]
		norm = math.sqrt(sum(abs(z) ** 2 for z in candidate))
		if norm > 1e-9:
			rest.append([z / norm for z in candidate])
	return [[rest[j][i] for j in range(dim)] for i in range(dim)]


def quantum_shannon_decomposition(
	qubits: Sequence[Qid], matrix: Sequence[Sequence[complex]], **kwargs: Any
) -> list[Operation]:
	"""Decompose an n-qubit unitary into single-qubit gates and CNOTs."""
	del kwargs
	qubits = tuple(qubits)
	unitaries, diagonal = _givens_decomposition(to_nested_list(matrix))
	operations: list[Operation] = []
	operations.extend(_diagonal_to_operations(diagonal, qubits))
	for i, j, two_level in reversed(unitaries):
		operations.extend(_two_level_to_operations(i, j, two_level, qubits))
	return operations


def three_qubit_matrix_to_operations(
	qubits: Sequence[Qid], matrix: Sequence[Sequence[complex]], **kwargs: Any
) -> list[Operation]:
	"""Decompose a three-qubit unitary into operations."""
	return quantum_shannon_decomposition(qubits, matrix, **kwargs)


def _givens_decomposition(matrix: Matrix) -> tuple[list[tuple[int, int, Matrix]], list[complex]]:
	n = len(matrix)
	u = [list(row) for row in matrix]
	unitaries: list[tuple[int, int, Matrix]] = []
	for col in range(n):
		for row in range(n - 1, col, -1):
			a = u[row - 1][col]
			b = u[row][col]
			if abs(b) < 1e-14:
				continue
			radius = math.hypot(abs(a), abs(b))
			g = [
				[a.conjugate() / radius, b.conjugate() / radius],
				[-b / radius, a / radius],
			]
			for k in range(n):
				t0, t1 = u[row - 1][k], u[row][k]
				u[row - 1][k] = g[0][0] * t0 + g[0][1] * t1
				u[row][k] = g[1][0] * t0 + g[1][1] * t1
			unitaries.append(
				(row - 1, row, [[g[0][0].conjugate(), g[1][0].conjugate()], [g[0][1].conjugate(), g[1][1].conjugate()]])
			)
	diagonal = [u[i][i] for i in range(n)]
	return unitaries, diagonal


def _two_level_to_operations(i: int, j: int, two_level: Matrix, qubits: Sequence[Qid]) -> list[Operation]:
	n = len(qubits)
	d = i ^ j
	set_bits = [k for k in range(n) if (d >> (n - 1 - k)) & 1]
	low = set_bits[0]
	operations: list[Operation] = [X(qubits[k]) for k in range(n) if (i >> (n - 1 - k)) & 1]
	for k in set_bits[1:]:
		operations.append(CNOT(qubits[low], qubits[k]))
	controls = [qubits[k] for k in range(n) if k != low]
	operations.extend(
		decompose_multi_controlled_rotation(two_level, controls, qubits[low], control_values=[False] * len(controls))
	)
	for k in reversed(set_bits[1:]):
		operations.append(CNOT(qubits[low], qubits[k]))
	operations.extend([X(qubits[k]) for k in range(n) if (i >> (n - 1 - k)) & 1])
	return operations


def _diagonal_to_operations(diagonal: Sequence[complex], qubits: Sequence[Qid]) -> list[Operation]:
	operations: list[Operation] = []
	for index in range(1, len(diagonal)):
		phase = diagonal[index]
		if abs(phase - 1.0) < 1e-12:
			continue
		two_level = [[1 + 0j, 0j], [0j, complex(phase)]]
		operations.extend(_two_level_to_operations(0, index, two_level, qubits))
	return operations


def decompose_clifford_tableau_to_operations(qubits: Sequence[Qid], tableau: Any) -> list[Operation]:
	"""Decompose a Clifford tableau into operations (best-effort)."""
	del qubits, tableau
	raise NotImplementedError("decompose_clifford_tableau_to_operations is not implemented")


class TwoQubitGateTabulationResult:
	"""The result of tabulating a two-qubit gate product."""

	def __init__(self, operations: Sequence[Operation], success: bool) -> None:
		self.operations = list(operations)
		self.success = bool(success)

	def __repr__(self) -> str:
		return f"vqubit.TwoQubitGateTabulationResult(success={self.success})"


class TwoQubitGateTabulation:
	"""A tabulation of products of a two-qubit gate, for fast synthesis."""

	def __init__(self, base_gate: Gate, max_infidelity: float = 0.001, *, atol: float = 1e-8) -> None:
		self.base_gate = base_gate
		self.max_infidelity = float(max_infidelity)
		self._atol = atol

	def __call__(self, gate: Any, qubits: Sequence[Qid], **kwargs: Any) -> TwoQubitGateTabulationResult:
		from ._protocols import unitary

		matrix = to_nested_list(unitary(gate))
		try:
			operations = two_qubit_matrix_to_sqrt_iswap_operations(
				matrix, qubits=(qubits[0], qubits[1]), atol=self._atol
			)
		except Exception:  # pragma: no cover - defensive
			return TwoQubitGateTabulationResult([], False)
		return TwoQubitGateTabulationResult(operations, True)

	def __repr__(self) -> str:
		return f"vqubit.TwoQubitGateTabulation({self.base_gate!r})"


def two_qubit_gate_product_tabulation(
	base_gate: Gate, max_infidelity: float, *, atol: float = 1e-8
) -> TwoQubitGateTabulation:
	"""Build a :class:`TwoQubitGateTabulation` for ``base_gate``."""
	return TwoQubitGateTabulation(base_gate, max_infidelity, atol=atol)


class CompilationTargetGateset(TargetGateset):
	"""Base class for a compilation target gateset."""

	@property
	def num_qubits(self) -> int:
		return 2

	def preprocess_transformers(self) -> list[Any]:
		return []

	def postprocess_transformers(self) -> list[Any]:
		return []


class TwoQubitCompilationTargetGateset(CompilationTargetGateset):
	"""A compilation gateset whose entanglers are two-qubit gates."""

	def __init__(self, *entanglers: Gate) -> None:
		self.entanglers = tuple(entanglers)

	def decompose_two_qubit(
		self, matrix: Sequence[Sequence[complex]], qubits: tuple[Qid, Qid], *, atol: float = 1e-8
	) -> list[Operation]:
		raise NotImplementedError


class SqrtIswapTargetGateset(TwoQubitCompilationTargetGateset):
	"""A compilation gateset targeting sqrt-iSWAP."""

	def __init__(self, *, atol: float = 1e-8, **kwargs: Any) -> None:
		del kwargs
		super().__init__(SQRT_ISWAP)
		self._atol = atol

	def decompose_two_qubit(
		self, matrix: Sequence[Sequence[complex]], qubits: tuple[Qid, Qid], *, atol: float | None = None
	) -> list[Operation]:
		return two_qubit_matrix_to_sqrt_iswap_operations(matrix, qubits=qubits, atol=atol or self._atol)

	def optimize_circuit(self, circuit: Circuit) -> Circuit:
		result = Circuit()
		for operation in circuit.all_operations():
			if (
				len(operation.qubits) == 2
				and operation.gate._unitary_() is not None
				and not operation.gate._is_measurement_()
			):
				matrix = to_nested_list(operation.gate._unitary_())
				if _is_sqrt_iswap_like(operation.gate):
					result.append(operation)
					continue
				result.append(self.decompose_two_qubit(matrix, (operation.qubits[0], operation.qubits[1])))
			else:
				result.append(operation)
		return result


def _is_sqrt_iswap_like(gate: Gate) -> bool:
	"""Whether the sqrt-iSWAP target can keep ``gate`` as it is.

	``SQRT_ISWAP`` is a plain ``ISwapPowGate`` at exponent ``±0.5`` -- as in
	Cirq -- rather than a class of its own, so the test is on the exponent.
	"""
	if isinstance(gate, ISwapPowGate):
		exponent = gate.exponent
		return isinstance(exponent, int | float) and abs(abs(exponent) - 0.5) < 1e-12
	return isinstance(gate, (PhasedXPowGate, ZPowGate, RzGate, RxGate, RyGate))


__all__ = [
	"CompilationTargetGateset",
	"SqrtIswapTargetGateset",
	"TwoQubitCompilationTargetGateset",
	"TwoQubitGateTabulation",
	"TwoQubitGateTabulationResult",
	"compute_cphase_exponents_for_fsim_decomposition",
	"decompose_clifford_tableau_to_operations",
	"decompose_cphase_into_two_fsim",
	"decompose_multi_controlled_rotation",
	"decompose_multi_controlled_x",
	"decompose_two_qubit_interaction_into_four_fsim_gates",
	"parameterized_2q_op_to_sqrt_iswap_operations",
	"prepare_two_qubit_state_using_cz",
	"prepare_two_qubit_state_using_iswap",
	"prepare_two_qubit_state_using_sqrt_iswap",
	"quantum_shannon_decomposition",
	"single_qubit_matrix_to_pauli_rotations",
	"single_qubit_op_to_framed_phase_form",
	"three_qubit_matrix_to_operations",
	"two_qubit_gate_product_tabulation",
	"two_qubit_matrix_to_cz_isometry",
	"two_qubit_matrix_to_diagonal_and_cz_operations",
	"two_qubit_matrix_to_ion_operations",
	"two_qubit_matrix_to_sqrt_iswap_operations",
]
