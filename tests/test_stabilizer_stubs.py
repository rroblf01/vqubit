"""Tests for the stabilizer helpers that used to be stubs.

Four things used to raise `NotImplementedError`: decomposing a Clifford tableau
into operations, the Kraus operators of a multi-qubit `ResetChannel`,
`CliffordTableau.kron`/`factor`/`reindex`, and applying an operation to the
CH-form simulation state.
"""

from __future__ import annotations

import random
from typing import Any

import pytest

import vqubit
from tests.helpers import assert_matrix_close


def test_a_tableau_decomposes_into_operations_that_rebuild_it() -> None:
	# A random Clifford circuit, its tableau, and the decomposition of that
	# tableau must prepare the same state.
	rng = random.Random(3)
	q = vqubit.LineQubit.range(3)
	tableau = vqubit.CliffordTableau(3)
	for _ in range(12):
		kind = rng.choice(["h", "x", "z", "cx", "swap"])
		first, second = rng.sample(range(3), 2)
		if kind == "h":
			tableau.apply_h(first)
		elif kind == "x":
			tableau.apply_x(first)
		elif kind == "z":
			tableau.apply_z(first)
		elif kind == "cx":
			tableau.apply_cx(first, second)
		else:
			tableau.apply_cx(first, second)
			tableau.apply_cx(second, first)
			tableau.apply_cx(first, second)

	operations = vqubit.decompose_clifford_tableau_to_operations(q, tableau)
	assert operations
	assert {type(operation.gate).__name__ for operation in operations} <= {
		"HPowGate",
		"XPowGate",
		"ZPowGate",
		"CXPowGate",
		"SwapPowGate",
	}
	rebuilt = vqubit.CliffordTableau(3)
	for operation in operations:
		_apply_to_tableau(rebuilt, operation)
	assert [str(stabilizer) for stabilizer in rebuilt.stabilizer_pauli_strings()] == [
		str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()
	]


def test_the_tableau_decomposition_validates_its_input() -> None:
	q = vqubit.LineQubit.range(2)
	with pytest.raises(ValueError):
		vqubit.decompose_clifford_tableau_to_operations([q[0]], vqubit.CliffordTableau(2))
	with pytest.raises(TypeError):
		vqubit.decompose_clifford_tableau_to_operations(q, object())


def _apply_to_tableau(tableau: vqubit.CliffordTableau, operation: Any) -> None:
	"""Replay one decomposed operation on a tableau."""
	gate = operation.gate
	axes = [qubit.x for qubit in operation.qubits]
	if isinstance(gate, vqubit.HPowGate):
		tableau.apply_h(axes[0])
	elif isinstance(gate, vqubit.XPowGate):
		tableau.apply_x(axes[0])
	elif isinstance(gate, vqubit.SwapPowGate):
		tableau.apply_cx(axes[0], axes[1])
		tableau.apply_cx(axes[1], axes[0])
		tableau.apply_cx(axes[0], axes[1])
	elif isinstance(gate, vqubit.CXPowGate):
		tableau.apply_cx(axes[0], axes[1])
	elif isinstance(gate, vqubit.ZPowGate):
		if abs(float(gate.exponent) - 1) < 1e-9:
			tableau.apply_z(axes[0])
		else:
			vqubit._compilation._apply_s(tableau, axes[0])
	else:  # pragma: no cover - the decomposition emits nothing else
		raise AssertionError(f"unexpected gate {type(gate).__name__}")


def test_reset_channel_kraus_operators_for_every_dimension() -> None:
	# A reset maps every basis state onto |0>, so the Kraus operators are the
	# ``dimension`` projections; their sum of squares is the identity (it is trace
	# preserving) and applying them to any state gives |0>.
	for dim in (2, 3, 4):
		operators = vqubit.kraus(vqubit.ResetChannel(dim))
		assert len(operators) == dim
		assert all(len(operator) == dim and len(operator[0]) == dim for operator in operators)
		# Each operator sends one basis state to |0...0>.
		for column, operator in enumerate(operators):
			image = [operator[row][column] for row in range(dim)]
			assert image[0] == 1
			assert all(entry == 0 for entry in image[1:])
		# Trace preservation is ``sum_k M_k^dagger M_k = I``.
		identity = [[1.0 if i == j else 0.0 for j in range(dim)] for i in range(dim)]
		total = [[0j] * dim for _ in range(dim)]
		for operator in operators:
			for i in range(dim):
				for j in range(dim):
					total[i][j] += sum(operator[row][i].conjugate() * operator[row][j] for row in range(dim))
		assert_matrix_close(total, identity)


def test_reset_channel_matches_the_single_qubit_case() -> None:
	one = vqubit.kraus(vqubit.ResetChannel())
	assert_matrix_close(one[0], [[1.0, 0.0], [0.0, 0.0]])
	assert_matrix_close(one[1], [[0.0, 1.0], [0.0, 0.0]])


def test_tableau_kron_builds_the_product_state() -> None:
	bell = vqubit.CliffordTableau(2)
	bell.apply_h(0)
	bell.apply_cx(0, 1)
	product = bell.kron(vqubit.CliffordTableau(1))
	assert product.n == 3
	assert [str(stabilizer) for stabilizer in product.stabilizer_pauli_strings()] == [
		"X(0)*X(1)",
		"Z(0)*Z(1)",
		"Z(2)",
	]
	with pytest.raises(TypeError):
		bell.kron(object())


def test_tableau_reindex_permutes_the_qubits() -> None:
	tableau = vqubit.CliffordTableau(2)
	tableau.apply_h(0)
	tableau.apply_x(1)
	assert sorted(str(stabilizer) for stabilizer in tableau.reindex([1, 0]).stabilizer_pauli_strings()) == [
		"-1*Z(0)",
		"X(1)",
	]
	with pytest.raises(ValueError):
		tableau.reindex([0, 0])
	with pytest.raises(ValueError):
		tableau.reindex([0])


def test_tableau_factor_splits_a_product_and_refuses_an_entangled_state() -> None:
	bell = vqubit.CliffordTableau(2)
	bell.apply_h(0)
	bell.apply_cx(0, 1)
	product = bell.kron(vqubit.CliffordTableau(1))
	assert product.supports_factor
	extracted, remainder = product.factor([0, 1])
	assert extracted == bell
	assert remainder.n == 1
	assert extracted.kron(remainder) == product
	with pytest.raises(ValueError):
		bell.factor([0])


def test_ch_form_state_applies_operations_and_measures() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	state = vqubit.StabilizerChFormSimulationState(qubits=[q0, q1], prng=random.Random(0))
	state.apply_operation(vqubit.H(q0))
	state.apply_operation(vqubit.CNOT(q0, q1))
	vector = state.state_vector()
	assert abs(abs(vector[0]) - 0.5**0.5) < 1e-9
	assert abs(abs(vector[3]) - 0.5**0.5) < 1e-9
	outcome = state.measure([q0], "m")
	assert outcome in (0, 1)
	assert state.log_of_measurement_results["m"] == outcome
	collapsed = state.state_vector()
	# The measurement collapsed the pair onto one of the two basis states.
	assert abs(abs(collapsed[0]) - 1) < 1e-9 or abs(abs(collapsed[3]) - 1) < 1e-9
	with pytest.raises(TypeError):
		state.apply_operation(vqubit.depolarize(0.1)(q0))
