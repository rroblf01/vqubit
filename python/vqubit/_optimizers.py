"""Circuit optimizers.

An optimizer rewrites a circuit into an equivalent one that is cheaper or
simpler. The first one here is the single-qubit Z-Y-Z decomposition: any
one-qubit unitary can be written as ``Rz(a) Ry(b) Rz(c)``, up to a global
phase.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Sequence

from ._gates import Gate, Ry, Rz
from ._ops import GateOperation
from ._qubits import Qid


def single_qubit_matrix_to_gates(matrix: Sequence[Sequence[complex]], *, atol: float = 1e-8) -> list[Gate]:
	"""Return ``Rz``/``Ry`` gates whose product is ``matrix``, up to phase.

	The gates are returned in application order, so feeding them to
	:class:`~vqubit.Circuit` reproduces the matrix.
	"""
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

	return [Rz(gamma), Ry(beta), Rz(alpha)]


def single_qubit_matrix_to_operations(
	matrix: Sequence[Sequence[complex]],
	qubit: Qid,
	*,
	atol: float = 1e-8,
) -> list[GateOperation]:
	"""Return the Z-Y-Z decomposition of ``matrix`` as operations on ``qubit``."""
	return [gate.on(qubit) for gate in single_qubit_matrix_to_gates(matrix, atol=atol)]


__all__ = ["single_qubit_matrix_to_gates", "single_qubit_matrix_to_operations"]
