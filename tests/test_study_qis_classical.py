"""Tests for the sweep helpers, the state helpers and the classical conditions.

These are the small utilities Cirq exports alongside the big types: building
sweeps from dictionaries, validating a state vector or a density matrix,
computing a Bloch vector or an entropy, and the `if(c == 1)` conditions of
classical control.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import vqubit
from tests.helpers import assert_matrix_close, assert_vector_close

q0, q1 = vqubit.LineQubit.range(2)


# --- sweeps ------------------------------------------------------------------


def test_unit_sweep_and_conversions() -> None:
	assert vqubit.UNIT_SWEEP.param_resolvers == [vqubit.ParamResolver({})]
	assert repr(vqubit.UNIT_SWEEP) == "vqubit.UnitSweep()"
	# A single sweep, a sequence of them and a bare dict all become a list.
	assert vqubit.to_sweeps(vqubit.UNIT_SWEEP) == [vqubit.UNIT_SWEEP]
	assert len(vqubit.to_sweeps([vqubit.UNIT_SWEEP, vqubit.Linspace("t", 0, 1, 2)])) == 2
	assert len(vqubit.flatten(vqubit.Linspace("t", 0, 1, 3))) == 3
	resolvers, params = vqubit.flatten_with_params(vqubit.Points("t", [0, 1]))
	assert len(resolvers) == 2
	assert params == {}
	resolvers, sweep = vqubit.flatten_with_sweep(vqubit.Points("t", [0, 1]))
	assert len(resolvers) == 2
	assert len(sweep.param_resolvers) == 2


def test_sweep_builders() -> None:
	product = vqubit.dict_to_product_sweep({"a": [0, 1], "b": [2]})
	assert len(product.param_resolvers) == 2
	zipped = vqubit.dict_to_zip_sweep({"a": [0, 1], "b": [2, 3]})
	assert [dict(resolver) for resolver in zipped.param_resolvers] == [
		{"a": 0, "b": 2},
		{"a": 1, "b": 3},
	]
	from_list = vqubit.list_of_dicts_to_zip([{"a": 0}, {"a": 1}])
	assert [dict(resolver) for resolver in from_list.param_resolvers] == [{"a": 0}, {"a": 1}]


def test_zip_longest_pads_the_shorter_sweep() -> None:
	zipped = vqubit.ZipLongest(vqubit.Points("a", [0, 1]), vqubit.Points("b", [5]))
	resolvers = zipped.param_resolvers
	assert len(resolvers) == 2
	assert dict(resolvers[0]) == {"a": 0, "b": 5}
	# The second point of the short sweep is an empty resolver.
	assert dict(resolvers[1]) == {"a": 1}
	with pytest.raises(ValueError, match="at least one sweep"):
		vqubit.ZipLongest()


def test_result_and_expression_maps() -> None:
	results = vqubit.ResultDict({"m": [[0]]})
	assert results["m"] == [[0]]
	expression = vqubit.ExpressionMap({"t": 0.5})
	assert expression["t"] == 0.5
	# A Parameter is looked up by its name.
	assert expression[vqubit.Parameter("t")] == 0.5


# --- states ------------------------------------------------------------------


def test_to_valid_state_vector() -> None:
	assert_vector_close(vqubit.to_valid_state_vector(0, (2,)), [1, 0])
	assert_vector_close(vqubit.to_valid_state_vector(1, (2,)), [0, 1])
	assert_vector_close(vqubit.to_valid_state_vector([0, 1], (2,)), [0, 1])
	# A non-normalized vector is rejected.
	with pytest.raises(ValueError, match="norm"):
		vqubit.to_valid_state_vector([2, 0], (2,))
	# The wrong length is rejected.
	with pytest.raises(ValueError, match="length"):
		vqubit.to_valid_state_vector([1, 0], (2, 2))
	# state_vector() infers the shape when it is not given.
	assert_vector_close(vqubit.state_vector([1, 0]), [1, 0])
	assert_vector_close(vqubit.state_vector(np.array([0, 1], dtype=complex)), [0, 1])


def test_to_valid_density_matrix() -> None:
	assert_matrix_close(vqubit.to_valid_density_matrix(0, (2,)), [[1, 0], [0, 0]])
	# A state vector becomes its projector.
	assert_matrix_close(vqubit.to_valid_density_matrix([1, 0], (2,)), [[1, 0], [0, 0]])
	mixed = [[0.5, 0], [0, 0.5]]
	assert_matrix_close(vqubit.to_valid_density_matrix(mixed, (2,)), mixed)
	with pytest.raises(ValueError, match="square"):
		vqubit.to_valid_density_matrix([[1, 0, 0], [0, 0, 0]], (2,))
	with pytest.raises(ValueError, match="trace"):
		vqubit.to_valid_density_matrix([[2, 0], [0, 0]], (2,))
	assert_matrix_close(vqubit.density_matrix(mixed), mixed)


def test_bloch_vector_of_a_known_state() -> None:
	zero: list[complex] = [1, 0]
	one: list[complex] = [0, 1]
	assert vqubit.bloch_vector_from_state_vector(zero, 0) == pytest.approx((0, 0, 1))
	assert vqubit.bloch_vector_from_state_vector(one, 0) == pytest.approx((0, 0, -1))
	# |+> points along x.
	plus: list[complex] = [1 / math.sqrt(2), 1 / math.sqrt(2)]
	assert vqubit.bloch_vector_from_state_vector(plus, 0) == pytest.approx((1, 0, 0), abs=1e-9)
	# The numpy path agrees with the list path.
	numpy_plus = np.array(plus)
	assert vqubit.bloch_vector_from_state_vector(numpy_plus, 0) == pytest.approx((1, 0, 0), abs=1e-9)
	with pytest.raises(ValueError, match="out of range"):
		vqubit.bloch_vector_from_state_vector(zero, 1)
	with pytest.raises(ValueError, match="out of range"):
		vqubit.bloch_vector_from_state_vector(np.array([1, 0]), 1)


def test_partial_trace_of_a_bell_state() -> None:
	bell = [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)]
	rho = np.outer(np.asarray(bell), np.asarray(bell).conj()).tolist()
	traced = vqubit.partial_trace(rho, [0])
	assert_matrix_close(np.asarray(traced), [[0.5, 0], [0, 0.5]])


def test_fidelity_and_entropy() -> None:
	assert vqubit.fidelity([1, 0], [1, 0]) == pytest.approx(1)
	assert vqubit.fidelity([1, 0], [0, 1]) == pytest.approx(0)
	# The mixed-state form agrees on pure states.
	assert vqubit.fidelity([[1, 0], [0, 0]], [[1, 0], [0, 0]]) == pytest.approx(1)
	# A pure state has zero entropy; the maximally mixed qubit has one bit.
	assert vqubit.von_neumann_entropy([[1, 0], [0, 0]]) == pytest.approx(0, abs=1e-9)
	assert vqubit.von_neumann_entropy([[0.5, 0], [0, 0.5]]) == pytest.approx(1)
	assert vqubit.von_neumann_entropy([[0.5, 0], [0, 0.5]], base=math.e) == pytest.approx(math.log(2))


def test_random_states_are_valid() -> None:
	state = vqubit.random_state_vector(2, seed=1)
	assert len(state) == 4
	assert sum(abs(z) ** 2 for z in state) == pytest.approx(1)
	# The same seed gives the same state.
	assert np.allclose(np.asarray(vqubit.random_state_vector(2, seed=1)), np.asarray(state))
	matrix = vqubit.random_density_matrix(2, seed=1)
	assert len(matrix) == 4
	assert np.trace(np.asarray(matrix)) == pytest.approx(1)
	assert_matrix_close(np.asarray(matrix), np.asarray(matrix).conj().T)


def test_kraus_channel_validation() -> None:
	operators = vqubit.kraus_channel([vqubit.unitary(vqubit.X)])
	assert len(operators) == 1
	assert_matrix_close(np.asarray(operators[0]), vqubit.unitary(vqubit.X))
	with pytest.raises(ValueError):
		vqubit.kraus_channel([])


# --- classical conditions ----------------------------------------------------


def test_key_conditions() -> None:
	condition = vqubit.KeyCondition("m", index=0)
	assert condition.key == "m"
	assert "m" in repr(condition)
	negated = ~condition
	assert str(negated) == "~" + str(condition)
	# The operations in a conditional circuit are classically controlled.
	circuit = vqubit.Circuit(
		vqubit.measure(q0, key="m"),
		vqubit.If(condition, vqubit.X(q1)),
	)
	assert len(list(circuit.all_operations())) == 2
	conditional = [operation for operation in circuit.all_operations() if "if" in str(operation)]
	assert conditional, str(circuit)


def test_set_variable_and_conditional_operation() -> None:
	set_variable = vqubit.SetVariable("m", 1)
	assert set_variable.key == "m"
	assert set_variable.value == 1
	assert "SetVariable" in str(set_variable)
	assert "SetVariable" in repr(set_variable)
	# The classically controlled operation carries its conditions and can be
	# remapped onto other qubits.
	condition = vqubit.KeyCondition("m", index=0)
	operation = vqubit.ClassicallyControlledOperation([condition], vqubit.X(q0))
	assert operation.conditions == (condition,) or len(operation.conditions) == 1
	assert "if" in str(operation)
	assert "ClassicallyControlledOperation" in repr(operation)
	remapped = operation.with_qubits(q1)
	assert remapped.qubits == (q1,)
	# An if/else builds two conditional operations.
	both = vqubit.If(condition, vqubit.X(q0), else_ops=[vqubit.Y(q0)])
	assert len(both) == 2
