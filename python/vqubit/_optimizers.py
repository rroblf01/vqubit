"""Circuit optimizers.

An optimizer rewrites a circuit into an equivalent one that is cheaper or
simpler. Two decompositions of a single-qubit unitary live here: the Z-Y-Z
form ``Rz(a) Ry(b) Rz(c)``, and the phased-X/Z form used by most hardware,
``PhasedX`` followed by ``Z``.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Sequence

from ._gates import Gate, PhasedXPowGate, PhasedXZGate, Ry, Rz, ZPowGate
from ._ops import GateOperation
from ._qubits import Qid

# Matrices and gates are compared with this tolerance when extracting angles.
_DEFAULT_ATOL = 1e-8


def _zyz_angles(matrix: Sequence[Sequence[complex]], atol: float = _DEFAULT_ATOL) -> tuple[float, float, float]:
	"""Return ``(alpha, beta, gamma)`` with ``U = Rz(a) Ry(b) Rz(c)`` up to phase."""
	u00 = complex(matrix[0][0])
	u01 = complex(matrix[0][1])
	u10 = complex(matrix[1][0])
	u11 = complex(matrix[1][1])

	det = u00 * u11 - u01 * u10
	if abs(det) > atol:
		root = cmath.sqrt(det)
		u00, u01, u10, u11 = u00 / root, u01 / root, u10 / root, u11 / root

	if abs(u00) < atol:
		beta = math.pi
		alpha = cmath.phase(u10)
		gamma = -alpha
	elif abs(u10) < atol:
		beta = 0.0
		alpha = 2 * cmath.phase(u11)
		gamma = 0.0
	else:
		beta = 2 * math.atan2(abs(u10), abs(u00))
		alpha = cmath.phase(u11) + cmath.phase(u10)
		gamma = cmath.phase(u11) - cmath.phase(u10)
	return alpha, beta, gamma


def single_qubit_matrix_to_gates(
	matrix: Sequence[Sequence[complex]],
	*,
	atol: float = _DEFAULT_ATOL,
) -> list[Gate]:
	"""Return ``Rz``/``Ry`` gates whose product is ``matrix``, up to phase.

	The gates are returned in application order, so feeding them to
	:class:`~vqubit.Circuit` reproduces the matrix.
	"""
	alpha, beta, gamma = _zyz_angles(matrix, atol)
	return [Rz(gamma), Ry(beta), Rz(alpha)]


def single_qubit_matrix_to_phased_x_z(
	matrix: Sequence[Sequence[complex]],
	*,
	atol: float = _DEFAULT_ATOL,
) -> list[Gate]:
	"""Return ``[PhasedXPowGate, ZPowGate]`` reproducing ``matrix`` up to phase."""
	alpha, beta, gamma = _zyz_angles(matrix, atol)
	return [
		PhasedXPowGate(phase_exponent=0.5 - gamma / math.pi, exponent=beta / math.pi),
		ZPowGate(exponent=(alpha + gamma) / math.pi),
	]


def single_qubit_matrix_to_phxz(
	matrix: Sequence[Sequence[complex]],
	*,
	atol: float = _DEFAULT_ATOL,
) -> PhasedXZGate:
	"""Return a :class:`~vqubit.PhasedXZGate` equal to ``matrix`` up to phase."""
	alpha, beta, gamma = _zyz_angles(matrix, atol)
	return PhasedXZGate(
		x_exponent=beta / math.pi,
		z_exponent=(alpha + gamma) / math.pi,
		axis_phase_exponent=0.5 - gamma / math.pi,
	)


def deconstruct_single_qubit_matrix_into_angles(
	matrix: Sequence[Sequence[complex]],
	*,
	atol: float = _DEFAULT_ATOL,
) -> tuple[float, float, float]:
	"""Return ``(x, y, z)`` half-turn angles for a ``PhasedXZGate`` equal to ``matrix``.

	The gate is ``PhasedXZGate(x_exponent=x, z_exponent=z, axis_phase_exponent=y)``.
	"""
	alpha, beta, gamma = _zyz_angles(matrix, atol)
	return (beta / math.pi, 0.5 - gamma / math.pi, (alpha + gamma) / math.pi)


def single_qubit_matrix_to_operations(
	matrix: Sequence[Sequence[complex]],
	qubit: Qid,
	*,
	atol: float = _DEFAULT_ATOL,
) -> list[GateOperation]:
	"""Return the Z-Y-Z decomposition of ``matrix`` as operations on ``qubit``."""
	return [gate.on(qubit) for gate in single_qubit_matrix_to_gates(matrix, atol=atol)]


__all__ = [
	"deconstruct_single_qubit_matrix_into_angles",
	"single_qubit_matrix_to_gates",
	"single_qubit_matrix_to_operations",
	"single_qubit_matrix_to_phased_x_z",
	"single_qubit_matrix_to_phxz",
]
