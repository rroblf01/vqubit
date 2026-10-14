"""Tests for the Ising (XX/YY/ZZ) power gates and iSWAP."""

from __future__ import annotations

import vqubit
from tests.helpers import assert_matrix_close
from vqubit._linalg import kron, matmul

X = [[0j, 1 + 0j], [1 + 0j, 0j]]
Y = [[0j, -1j], [1j, 0j]]
Z = [[1 + 0j, 0j], [0j, -1 + 0j]]


def test_xx_is_x_tensor_x() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.XX), kron(X, X))


def test_yy_is_y_tensor_y() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.YY), kron(Y, Y))


def test_zz_is_z_tensor_z() -> None:
	assert_matrix_close(vqubit.unitary(vqubit.ZZ), kron(Z, Z))


def test_pauli_product_squares_to_identity() -> None:
	for gate in [vqubit.XX, vqubit.YY, vqubit.ZZ]:
		identity = [[1 + 0j if i == j else 0j for j in range(4)] for i in range(4)]
		assert_matrix_close(vqubit.unitary(gate**2), identity)


def test_pauli_product_power() -> None:
	# The square root squares back to the gate.
	sqrt = vqubit.unitary(vqubit.XX**0.5)
	squared = matmul(sqrt, sqrt)
	assert_matrix_close(squared, vqubit.unitary(vqubit.XX))


def test_pauli_product_can_be_parameterized() -> None:
	t = vqubit.Parameter("t")
	gate = vqubit.ZZPowGate(exponent=t)
	assert vqubit.is_parameterized(gate)
	resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({t: 1}))
	assert_matrix_close(vqubit.unitary(resolved), vqubit.unitary(vqubit.ZZ))


def test_iswap_unitary() -> None:
	expected = [
		[1 + 0j, 0j, 0j, 0j],
		[0j, 0j, 1j, 0j],
		[0j, 1j, 0j, 0j],
		[0j, 0j, 0j, 1 + 0j],
	]
	assert_matrix_close(vqubit.unitary(vqubit.ISWAP), expected)


def test_iswap_swaps_with_phase() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	# |10> becomes i|01>, so the amplitude moves from index 2 to index 1.
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.ISWAP(q0, q1))
	result = vqubit.Simulator().simulate(circuit)
	vector = result.state_vector()
	assert abs(vector[1] - 1j) < 1e-9
	assert abs(vector[2]) < 1e-9
