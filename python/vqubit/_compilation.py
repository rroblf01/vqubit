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
	SWAP,
	FSimGate,
	Gate,
	H,
	MatrixGate,
	PhasedXPowGate,
	Rx,
	RxGate,
	Ry,
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
from ._parameter import is_symbolic
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
	matrix: Sequence[Sequence[complex]], *, atol: float = 1e-8
) -> tuple[Matrix, complex, complex]:
	"""Decompose a 2x2 unitary as ``U^-1 · diag(1, r) · U · diag(g, g)``.

	``U`` rotates the rotation axis of ``matrix`` onto the Z axis, ``r`` is the
	relative phase (the rotation about that axis) and ``g`` the global phase.
	This is ``cirq.single_qubit_op_to_framed_phase_form``, and it is what a
	controlled single-qubit rotation is decomposed with.

	Returns ``(U, r, g)``.
	"""
	from ._linalg_public import unitary_eig

	nested = to_nested_list(matrix)
	if len(nested) != 2 or any(len(row) != 2 for row in nested):
		raise ValueError(f"Expected a 2x2 matrix, got shape {len(nested)}x{len(nested[0]) if nested else 0}.")
	values, vectors = unitary_eig(nested, atol=atol)
	# ``unitary_eig`` returns the eigenvectors as columns, so ``U`` is V^dagger.
	u = [[complex(vectors[row][column]).conjugate() for row in range(2)] for column in range(2)]
	relative_phase = complex(values[1]) / complex(values[0])
	global_phase = complex(values[0])
	return u, relative_phase, global_phase


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
		operations: list[Operation] = [*flips]
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
	standard_operations: list[Operation] = [*pre, MatrixGate(c_matrix).on(target), CNOT(control, target)]
	standard_operations.append(MatrixGate(b_matrix).on(target))
	standard_operations.append(CNOT(control, target))
	standard_operations.append(MatrixGate(a_matrix).on(target))
	if abs(phi) > 1e-12:
		standard_operations.append(ZPowGate(exponent=phi / math.pi).on(control))
	standard_operations.extend(post)
	return standard_operations


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


def _diagonal_gate_of(matrix: Sequence[Sequence[complex]]) -> Any:
	"""A ``DiagonalGate`` applying the phases on the diagonal of ``matrix``."""
	from ._gates_more import DiagonalGate

	return DiagonalGate([cmath.phase(complex(matrix[i][i])) for i in range(len(matrix))])


def two_qubit_matrix_to_diagonal_and_cz_operations(
	q0: Qid,
	q1: Qid,
	mat: Sequence[Sequence[complex]],
	allow_partial_czs: bool = False,
	atol: float = 1e-8,
	clean_operations: bool = True,
) -> tuple[Any, list[Operation]]:
	"""Return ``(diagonal, operations)`` with ``mat ≈ Circuit(operations) @ diagonal``.

	The product holds up to a global phase, as in Cirq. A unitary that needs
	three CNOTs is rewritten as ``mat @ D`` (which needs two) followed by the
	diagonal ``D^-1``, so the CZ decomposition is cheaper;
	a unitary that already needs two or fewer keeps an identity diagonal, and a
	diagonal unitary needs no operations at all. This is
	``cirq.two_qubit_matrix_to_diagonal_and_cz_operations``.
	"""
	from ._linalg_public import extract_right_diag, is_diagonal, is_unitary, num_cnots_required

	nested = to_nested_list(mat)
	cols = len(nested[0]) if nested else 0
	if len(nested) != 4 or cols != 4:
		raise ValueError(f"Expected a 4x4 unitary, got shape {len(nested)}x{cols}.")
	if not is_unitary(nested, atol=max(atol, 1e-8)):
		raise ValueError("Expected a 4x4 unitary matrix.")
	if is_diagonal(nested, atol=atol):
		return _diagonal_gate_of(nested), []
	if num_cnots_required(nested, atol=atol) == 3:
		right = extract_right_diag(nested, atol=atol)
		two_cnot = matmul(nested, right)
		operations = two_qubit_matrix_to_cz_operations(
			q0, q1, two_cnot, allow_partial_czs, atol=atol, clean_operations=clean_operations
		)
		# ``mat = ops @ right``, so the diagonal returned is ``right^-1``.
		return _diagonal_gate_of(_dag(right)), operations
	operations = two_qubit_matrix_to_cz_operations(
		q0, q1, nested, allow_partial_czs, atol=atol, clean_operations=clean_operations
	)
	return _diagonal_gate_of([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), operations


def two_qubit_matrix_to_cz_isometry(
	q0: Qid,
	q1: Qid,
	mat: Sequence[Sequence[complex]],
	allow_partial_czs: bool = False,
	atol: float = 1e-8,
	clean_operations: bool = True,
) -> list[Operation]:
	"""Decompose a two-qubit isometry into single-qubit gates and CZs."""
	nested = to_nested_list(mat)
	# A two-qubit isometry is 4 rows of one to four orthonormal columns.
	cols = len(nested[0]) if nested else 0
	if len(nested) != 4 or not 1 <= cols <= 4 or any(len(row) != cols for row in nested):
		raise ValueError(f"Expected a 4xN isometry (1 <= N <= 4), got shape {len(nested)}x{cols}.")
	full = _complete_isometry(nested)
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


def _asinsin(value: float) -> float:
	"""``arcsin(sin(value))``, returned in ``[-π/2, π/2]``."""
	turns = round(value / math.pi)
	if turns % 2 == 0:
		return value - turns * math.pi
	return turns * math.pi - value


def compute_cphase_exponents_for_fsim_decomposition(fsim_gate: FSimGate) -> list[tuple[float, float]]:
	"""Return the intervals of ``CZPowGate`` exponents a two-fSim decomposition allows.

	This is ``cirq.compute_cphase_exponents_for_fsim_decomposition``: the fSim's
	two angles bound the phase angle, which gives one or two intervals (cropped
	to ``[0, 2]``). An empty list means the gate cannot decompose any CZ at all.
	"""
	theta, phi = fsim_gate.theta, fsim_gate.phi
	if is_symbolic(theta) or is_symbolic(phi):
		raise ValueError("Symbolic arguments are not supported by the fSim decomposition")
	bound_theta = abs(_asinsin(float(theta)))
	bound_phi = abs(_asinsin(float(phi) / 2))
	low = 4 * min(bound_theta, bound_phi) / math.pi
	high = 4 * max(bound_theta, bound_phi) / math.pi
	# The sine's second "hump": the same width, mirrored about one turn.
	mirror_low, mirror_high = 2 - high, 2 - low
	intervals: list[tuple[float, float]] = []
	if high < mirror_low:
		intervals = [(low, high), (mirror_low, mirror_high)]
	elif mirror_high < low:
		intervals = [(mirror_low, mirror_high), (low, high)]
	else:
		intervals = [(min(low, mirror_low), max(high, mirror_high))]
	return [(start, stop) for start, stop in intervals if start < stop and stop > 0 and start < 2]


def decompose_cphase_into_two_fsim(
	cphase_gate: Any,
	*,
	fsim_gate: FSimGate,
	qubits: Sequence[Qid] | None = None,
	atol: float = 1e-8,
) -> list[Operation]:
	"""Decompose a ``CZPowGate`` into two fSim gates and single-qubit rotations.

	The construction is section VII F of https://arxiv.org/abs/1910.11333 and
	matches ``cirq.decompose_cphase_into_two_fsim``, including the global phase.
	It is feasible only when the fSim angles and the CZ exponent satisfy

	``|sin(theta)| <= |sin(delta/4)| <= |sin(phi/2)|`` (or the reverse order),

	with ``delta = -pi * cphase_gate.exponent``; otherwise ``ValueError`` lists
	the exponents that would work.
	"""
	from ._gates import global_phase

	if is_symbolic(getattr(cphase_gate, "exponent", 1.0)):
		raise ValueError("Cannot decompose a parametrized gate.")
	if is_symbolic(fsim_gate.theta) or is_symbolic(fsim_gate.phi):
		raise ValueError("Cannot decompose into a parametrized gate.")
	if qubits is None:
		qubits = LineQubit.range(2)
	if len(qubits) != 2:
		raise ValueError(f"Expected a pair of qubits, but got {qubits!r}.")
	q0, q1 = qubits

	theta = float(fsim_gate.theta)
	phi = float(fsim_gate.phi)
	sin_half_phi = math.sin(phi / 2)
	cos_half_phi = math.cos(phi / 2)
	sin_theta = math.sin(theta)
	cos_theta = math.cos(theta)

	denominator = (sin_theta - sin_half_phi) * (sin_theta + sin_half_phi)
	if abs(denominator) < atol:
		raise ValueError(
			f"{fsim_gate} cannot decompose a CZPowGate: sin(theta)**2 is too close to "
			f"sin(phi/2)**2 (difference {denominator})."
		)

	exponent = float(getattr(cphase_gate, "exponent", 1.0))
	canonical_delta = -math.pi * (exponent % 2)
	delta = canonical_delta
	for candidate in (canonical_delta, canonical_delta + 2 * math.pi):
		sin_quarter = math.sin(candidate / 4)
		numerator = (sin_quarter - sin_half_phi) * (sin_quarter + sin_half_phi)
		sin_alpha_squared = numerator / denominator
		if 0 <= sin_alpha_squared <= 1:
			delta = candidate
			break
	else:
		intervals = compute_cphase_exponents_for_fsim_decomposition(fsim_gate)
		raise ValueError(
			f"{cphase_gate} cannot be decomposed into two {fsim_gate}. "
			f"Valid intervals for the canonical CZPowGate exponent: {intervals}."
		)
	alpha = math.asin(math.sqrt(sin_alpha_squared))

	tan_alpha = math.tan(alpha)
	xi = math.atan2(tan_alpha * cos_theta, cos_half_phi)
	eta = math.atan2(tan_alpha * sin_theta, sin_half_phi)
	if delta < 0:
		eta += math.pi

	phase = global_phase(phi / (4 * math.pi))
	return [
		Rx(xi)(q0),
		Rx(eta)(q1),
		fsim_gate.on(q0, q1),
		Rz(phi / 2)(q0),
		Rz(phi / 2)(q1),
		phase.on(),
		Rx(-2 * alpha)(q0),
		Z(q0),
		fsim_gate.on(q0, q1),
		Rz(phi / 2)(q0),
		Rz(phi / 2)(q1),
		phase.on(),
		Z(q0),
		Rx(-eta)(q1),
		Rx(xi)(q0),
		Rz(-delta / 2)(q0),
		Rz(-delta / 2)(q1),
		global_phase(-delta / (4 * math.pi)).on(),
	]


def _sticky_0_to_1(value: float, *, atol: float) -> float | None:
	"""``value`` clamped into ``[0, 1]`` when it is within ``atol`` of it."""
	if 0 <= value <= 1:
		return value
	if 1 < value <= 1 + atol:
		return 1.0
	if 0 > value >= -atol:
		return 0.0
	return None


def _dag(matrix: Sequence[Sequence[complex]]) -> Matrix:
	return [[complex(matrix[j][i]).conjugate() for j in range(len(matrix))] for i in range(len(matrix))]


class _BGate(Gate):
	"""The two-qubit gate two copies of which reach any KAK interaction."""

	def _num_qubits_(self) -> int:
		return 2

	def _decompose_(self, qubits: Sequence[Qid]) -> list[Operation]:
		from ._gates import XX, YY

		a, b = qubits
		return [(XX**-0.5)(a, b), (YY**-0.25)(a, b)]


def _decompose_xx_yy_into_two_fsims_ignoring_single_qubit_ops(
	*,
	qubits: Sequence[Qid],
	fsim_gate: FSimGate,
	canonical_x_kak_coefficient: float,
	canonical_y_kak_coefficient: float,
	atol: float = 1e-8,
) -> list[Operation]:
	x = canonical_x_kak_coefficient
	y = canonical_y_kak_coefficient
	assert 0 <= y <= x <= math.pi / 4

	eta = math.sin(x) ** 2 * math.cos(y) ** 2 + math.cos(x) ** 2 * math.sin(y) ** 2
	xi = abs(math.sin(2 * x) * math.sin(2 * y))

	t = float(fsim_gate.phi) / 2
	kappa = math.sin(float(fsim_gate.theta)) ** 2 - math.sin(t) ** 2
	s_sum = (eta - math.sin(t) ** 2) / kappa
	s_dif = 0.5 * xi / kappa

	a_dif = _sticky_0_to_1(s_sum + s_dif, atol=atol)
	a_sum = _sticky_0_to_1(s_sum - s_dif, atol=atol)
	if a_dif is None or a_sum is None:
		raise ValueError(
			f"Failed to synthesize XX^{x / math.pi}·YY^{y / math.pi} from two {fsim_gate!r} "
			"separated by single-qubit operations."
		)

	x_dif = math.asin(math.sqrt(a_dif))
	x_sum = math.asin(math.sqrt(a_sum))
	x_a, x_b = x_sum + x_dif, x_dif - x_sum
	a, b = qubits
	return [
		fsim_gate.on(a, b),
		Rz(t + math.pi)(a),
		Rz(t)(b),
		Rx(x_a)(a),
		Rx(x_b)(b),
		fsim_gate.on(a, b),
	]


def _decompose_interaction_into_two_b_gates_ignoring_single_qubit_ops(
	qubits: Sequence[Qid], kak_interaction_coefficients: Sequence[float]
) -> list[Operation]:
	"""The minimal two-B construction of a two-qubit operation (quant-ph/0312193)."""
	a, b = qubits
	x, y, z = kak_interaction_coefficients
	r = max(0.0, (math.sin(y) * math.cos(z)) ** 2)
	if r > 0.499999999999:
		rb: list[Operation] = [Ry(math.pi)(b)]
	else:
		b1 = max(0.0, min(1.0, math.cos(y * 2) * math.cos(z * 2) / (1 - 2 * r)))
		b2 = math.asin(math.sqrt(b1))
		b3 = math.acos(1 - 4 * r)
		rb = [Rz(-b2)(b), Ry(-b3)(b), Rz(-b2)(b)]
	s = 1 if z < 0 else -1
	return [_BGate().on(a, b), Ry(s * 2 * x)(a), *rb, _BGate().on(a, b)]


def _fix_single_qubit_gates_around_kak_interaction(
	*,
	desired: Any,
	operations: list[Operation],
	qubits: Sequence[Qid],
) -> list[Operation]:
	"""Complete a KAK interaction with the single-qubit gates around it."""
	from ._protocols import unitary as _unitary

	actual = kak_decomposition(_unitary(Circuit(operations)))
	result: list[Operation] = []
	for k in range(2):
		gate = MatrixGate(
			matmul(_dag(actual.single_qubit_operations_before[k]), desired.single_qubit_operations_before[k])
		)
		result.append(gate.on(qubits[k]))
	result.extend(operations)
	for k in range(2):
		gate = MatrixGate(
			matmul(desired.single_qubit_operations_after[k], _dag(actual.single_qubit_operations_after[k]))
		)
		result.append(gate.on(qubits[k]))
	phase = complex(desired.global_phase) / complex(actual.global_phase)
	result.append(global_phase_for(phase).on())
	return result


def decompose_two_qubit_interaction_into_four_fsim_gates(
	interaction: Any,
	*,
	fsim_gate: Any,
	qubits: Sequence[Qid] | None = None,
	atol: float = 1e-8,
) -> Circuit:
	"""Decompose a two-qubit interaction into exactly four fSim gates.

	This is ``cirq.decompose_two_qubit_interaction_into_four_fsim_gates``: the
	interaction becomes two ``B`` gates and each ``B`` gate becomes two copies
	of ``fsim_gate`` separated by single-qubit rotations. It only works for an
	fSim close to iSWAP strength (``3π/8 ≤ |theta| ≤ 5π/8``) and a small cphase
	(``|phi| ≤ π/4``).
	"""
	from ._gates import ISwapPowGate

	if is_symbolic(getattr(fsim_gate, "theta", None)) or is_symbolic(getattr(fsim_gate, "phi", None)):
		raise ValueError("FSimGate must not have parameterized values for angles.")
	mapped_gate = (
		FSimGate(-float(fsim_gate.exponent) * math.pi / 2, 0.0) if isinstance(fsim_gate, ISwapPowGate) else fsim_gate
	)
	theta = float(mapped_gate.theta)
	phi = float(mapped_gate.phi)
	if not 3 / 8 * math.pi <= abs(theta) <= 5 / 8 * math.pi:
		raise ValueError("Must have 3π/8 ≤ |fsim_gate.theta| ≤ 5π/8")
	if abs(phi) > math.pi / 4:
		raise ValueError("Must have abs(fsim_gate.phi) ≤ π/4")
	if qubits is None:
		from ._ops import Operation

		qubits = interaction.qubits if isinstance(interaction, Operation) else LineQubit.range(2)
	if len(qubits) != 2:
		raise ValueError(f"Expected a pair of qubits, but got {qubits!r}.")

	from ._protocols import unitary as _unitary

	try:
		interaction_matrix = _unitary(interaction)
	except TypeError:
		interaction_matrix = interaction
	kak = kak_decomposition(interaction_matrix)
	two_b = _decompose_interaction_into_two_b_gates_ignoring_single_qubit_ops(qubits, kak.interaction_coefficients)
	two_b = _fix_single_qubit_gates_around_kak_interaction(desired=kak, operations=two_b, qubits=qubits)

	b_kak = kak_decomposition(_unitary((_BGate()).on(*qubits)))
	b_as_fsims = _decompose_xx_yy_into_two_fsims_ignoring_single_qubit_ops(
		qubits=qubits,
		fsim_gate=mapped_gate,
		canonical_x_kak_coefficient=b_kak.interaction_coefficients[0],
		canonical_y_kak_coefficient=b_kak.interaction_coefficients[1],
		atol=atol,
	)
	b_as_fsims = _fix_single_qubit_gates_around_kak_interaction(desired=b_kak, operations=b_as_fsims, qubits=qubits)
	if mapped_gate is not fsim_gate:
		b_as_fsims = [
			fsim_gate.on(*operation.qubits) if operation.gate == mapped_gate else operation for operation in b_as_fsims
		]

	result: list[Operation] = []
	for operation in two_b:
		if isinstance(operation.gate, _BGate):
			result.extend(b_as_fsims)
		else:
			result.append(operation)
	return Circuit(result)


def global_phase_for(coefficient: complex) -> Any:
	"""A ``GlobalPhaseGate`` applying ``coefficient`` (an arbitrary unit phase)."""
	from ._gates import global_phase

	return global_phase(cmath.phase(coefficient) / math.pi)


def prepare_two_qubit_state_using_cz(
	gate: Any, qubits: Sequence[Qid], *, atol: float = 1e-8, **kwargs: Any
) -> list[Operation]:
	"""Prepare a two-qubit state using CZs."""
	del kwargs
	matrix = _state_to_unitary(gate, atol=atol)
	return two_qubit_matrix_to_cz_operations(qubits[0], qubits[1], matrix, atol=atol)


def prepare_two_qubit_state_using_iswap(
	gate: Any, qubits: Sequence[Qid], *, atol: float = 1e-8, **kwargs: Any
) -> list[Operation]:
	"""Prepare a two-qubit state using iSWAPs."""
	del kwargs
	matrix = _state_to_unitary(gate, atol=atol)
	return two_qubit_matrix_to_sqrt_iswap_operations(matrix, qubits=(qubits[0], qubits[1]), atol=atol)


def prepare_two_qubit_state_using_sqrt_iswap(
	gate: Any, qubits: Sequence[Qid], *, atol: float = 1e-8, **kwargs: Any
) -> list[Operation]:
	"""Prepare a two-qubit state using sqrt-iSWAPs."""
	return prepare_two_qubit_state_using_iswap(gate, qubits, atol=atol, **kwargs)


def _state_to_unitary(state: Any, *, atol: float = 1e-8) -> Matrix:
	vector = [complex(z) for z in state]
	dim = len(vector)
	if dim != 4:
		raise ValueError(f"Expected a two-qubit state with 4 amplitudes, got {dim}.")
	norm = math.sqrt(sum(abs(z) ** 2 for z in vector))
	if abs(norm - 1.0) > max(atol, 1e-6):
		raise ValueError(f"Expected a normalized state, got norm {norm}.")
	if dim != 4:
		raise ValueError
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


def _apply_s(tableau: Any, axis: int, *, dagger: bool = False) -> None:
	"""Conjugate one qubit of ``tableau`` by ``S`` (or ``S†``).

	``S`` sends ``X`` to ``Y`` and ``Y`` to ``-X``, leaving ``Z`` alone; ``S†``
	does the opposite, ``Y`` to ``X`` and ``X`` to ``-Y``. In the tableau's
	``(x, z)`` bits the ``z`` bit flips wherever the ``x`` bit is set, and the row
	picks up a sign depending on which of the two it was.
	"""
	for row in range(2 * tableau.n):
		x = tableau._xs[row][axis]
		z = tableau._zs[row][axis]
		if x and (not z if dagger else z):
			tableau._rs[row] ^= 1
		tableau._zs[row][axis] = z ^ x


def decompose_clifford_tableau_to_operations(qubits: Sequence[Qid], tableau: Any) -> list[Operation]:
	"""Decompose a Clifford tableau into one- and two-qubit operations.

	The procedure is the one in Theorem 8 of Aaronson and Gottesman, *Improved
	Simulation of Stabilizer Circuits* (https://arxiv.org/abs/quant-ph/0406196),
	row by row, as Cirq implements it. Applying the returned operations to a
	tableau that starts in ``|0...0>`` reproduces ``tableau``.
	"""
	from ._clifford import CliffordTableau

	if not isinstance(tableau, CliffordTableau):
		raise TypeError(
			f"decompose_clifford_tableau_to_operations needs a CliffordTableau, got {type(tableau).__name__}"
		)
	if len(qubits) != tableau.n:
		raise ValueError("the number of qubits must match the size of the tableau")

	state = tableau.copy()
	size = state.n
	operations: list[Operation] = []

	def hadamard(axis: int) -> None:
		state.apply_h(axis)
		operations.append(H.on(qubits[axis]))

	def pauli_x(axis: int) -> None:
		state.apply_x(axis)
		operations.append(X.on(qubits[axis]))

	def pauli_z(axis: int) -> None:
		state.apply_z(axis)
		operations.append(Z.on(qubits[axis]))

	def phase_dagger(axis: int) -> None:
		# ``S†`` goes into the tableau; ``S`` is what the inverted list needs.
		_apply_s(state, axis, dagger=True)
		operations.append(S.on(qubits[axis]))

	def cnot(control: int, target: int) -> None:
		state.apply_cx(control, target)
		operations.append(CNOT.on(qubits[control], qubits[target]))

	def swap(first: int, second: int) -> None:
		# Apply the three CNOTs that make a SWAP to the tableau, but record a
		# single SWAP: recording the CNOTs too would apply it twice.
		state.apply_cx(first, second)
		state.apply_cx(second, first)
		state.apply_cx(first, second)
		operations.append(SWAP.on(qubits[first], qubits[second]))

	xs, zs, rs = state._xs, state._zs, state._rs
	for i in range(size):
		# 1: put a 1 on the diagonal of the destabilizer block.
		if not xs[i][i] and zs[i][i]:
			hadamard(i)
		if not xs[i][i]:
			for j in range(i + 1, size):
				if xs[i][j]:
					swap(i, j)
					break
		if not xs[i][i]:
			for j in range(i + 1, size):
				if zs[i][j]:
					hadamard(j)
					swap(i, j)
					break
		# 2: clear the rest of that destabilizer row.
		for j in range(i + 1, size):
			if xs[i][j]:
				cnot(i, j)
		if any(zs[i][j] for j in range(i, size)):
			if not zs[i][i]:
				phase_dagger(i)
			for j in range(i + 1, size):
				if zs[i][j]:
					cnot(j, i)
			phase_dagger(i)
		# 3: the matching stabilizer row.
		for j in range(i + 1, size):
			if zs[i + size][j]:
				cnot(j, i)
		if any(xs[i + size][j] for j in range(i, size)):
			hadamard(i)
			for j in range(i + 1, size):
				if xs[i + size][j]:
					cnot(i, j)
			if zs[i + size][i]:
				phase_dagger(i)
			hadamard(i)
	# 4: the signs of the two blocks.
	for i in range(size):
		if rs[i]:
			pauli_z(i)
		if rs[size + i]:
			pauli_x(i)
	# 5: undo the procedure, which is (AB)^-1 = B^-1 A^-1; only ``S`` is not
	# self-inverse, and ``phase_dagger`` recorded its inverse on purpose.
	return operations[::-1]


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
		"""The number of qubits the gate set targets."""
		return 2

	def preprocess_transformers(self) -> list[Any]:
		"""The transformers to run before compiling with this gate set."""
		return []

	def postprocess_transformers(self) -> list[Any]:
		"""The transformers to run after compiling with this gate set."""
		return []


class TwoQubitCompilationTargetGateset(CompilationTargetGateset):
	"""A compilation gateset whose entanglers are two-qubit gates."""

	def __init__(self, *entanglers: Gate) -> None:
		self.entanglers = tuple(entanglers)

	def decompose_two_qubit(
		self, matrix: Sequence[Sequence[complex]], qubits: tuple[Qid, Qid], *, atol: float = 1e-8
	) -> list[Operation]:
		"""Decompose a two-qubit gate into the target gate set."""
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
		"""Decompose a two-qubit gate into the square root of iSWAP."""
		return two_qubit_matrix_to_sqrt_iswap_operations(matrix, qubits=qubits, atol=atol or self._atol)

	def optimize_circuit(self, circuit: Circuit) -> Circuit:
		"""Compile a circuit onto this gate set."""
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
