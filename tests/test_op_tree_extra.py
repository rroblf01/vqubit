"""Additional tests for the op-tree helpers, gate families and gatesets."""

from __future__ import annotations

import pytest

import vqubit


def test_none_is_an_empty_op_tree() -> None:
	assert vqubit.flatten_to_ops(None) == []
	assert vqubit.flatten_to_ops_or_moments(None) == []
	assert vqubit.flatten_op_tree(None) == []
	assert vqubit.freeze_op_tree(None) == ()


def test_an_operation_flattens_to_itself() -> None:
	q0 = vqubit.LineQubit(0)
	operation = vqubit.H(q0)

	assert vqubit.flatten_to_ops(operation) == [operation]
	assert vqubit.flatten_to_ops_or_moments(operation) == [operation]


def test_a_moment_flattens_to_its_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	moment = vqubit.Moment([vqubit.H(q0), vqubit.X(q1)])

	assert vqubit.flatten_to_ops(moment) == [vqubit.H(q0), vqubit.X(q1)]
	assert vqubit.flatten_to_ops_or_moments(moment) == [moment]


def test_a_circuit_flattens_in_time_order() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q1), vqubit.CNOT(q0, q1))

	assert vqubit.flatten_to_ops(circuit) == [vqubit.H(q0), vqubit.X(q1), vqubit.CNOT(q0, q1)]
	assert vqubit.flatten_to_ops_or_moments(circuit) == list(circuit.moments)


def test_nested_iterables_and_none_are_flattened_depth_first() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	tree = [vqubit.H(q0), [vqubit.X(q1), None], (operation for operation in [vqubit.CNOT(q0, q1)])]

	assert vqubit.flatten_to_ops(tree) == [vqubit.H(q0), vqubit.X(q1), vqubit.CNOT(q0, q1)]


def test_strings_are_not_op_trees() -> None:
	with pytest.raises(TypeError):
		vqubit.flatten_to_ops("H(q0)")
	with pytest.raises(TypeError):
		vqubit.flatten_to_ops_or_moments(b"H(q0)")


def test_flatten_to_ops_or_moments_keeps_moments_and_operations() -> None:
	q0 = vqubit.LineQubit(0)
	moment = vqubit.Moment([vqubit.H(q0)])
	operation = vqubit.X(q0)

	assert vqubit.flatten_to_ops_or_moments([moment, operation]) == [moment, operation]


def test_flatten_to_ops_or_moments_rejects_non_op_trees() -> None:
	with pytest.raises(TypeError):
		vqubit.flatten_to_ops_or_moments(3.5)


def test_flatten_op_tree_defaults_to_operations() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	tree = [vqubit.Moment([vqubit.H(q0), vqubit.X(q1)])]

	assert vqubit.flatten_op_tree(tree) == [vqubit.H(q0), vqubit.X(q1)]
	assert isinstance(vqubit.flatten_op_tree(tree), list)


def test_flatten_op_tree_can_pair_operations_with_their_moment() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	first = vqubit.Moment([vqubit.H(q0), vqubit.X(q1)])
	second = vqubit.Moment([vqubit.CNOT(q0, q1)])

	pairs = vqubit.flatten_op_tree([first, second], preserve_moments=True)

	assert pairs == [
		(vqubit.H(q0), first),
		(vqubit.X(q1), first),
		(vqubit.CNOT(q0, q1), second),
	]


def test_flatten_op_tree_wraps_standalone_operations_in_a_new_moment() -> None:
	q0 = vqubit.LineQubit(0)
	operation = vqubit.H(q0)

	[(paired_operation, moment)] = vqubit.flatten_op_tree(operation, preserve_moments=True)

	assert paired_operation == operation
	assert moment == vqubit.Moment([operation])


def test_freeze_op_tree_is_a_flat_tuple() -> None:
	q0 = vqubit.LineQubit(0)
	frozen = vqubit.freeze_op_tree([vqubit.H(q0), [vqubit.X(q0)]])

	assert isinstance(frozen, tuple)
	assert frozen == (vqubit.H(q0), vqubit.X(q0))


def test_freeze_op_tree_can_preserve_moments() -> None:
	q0 = vqubit.LineQubit(0)
	moment = vqubit.Moment([vqubit.H(q0)])

	frozen = vqubit.freeze_op_tree(moment, preserve_moments=True)

	assert frozen == ((vqubit.H(q0), moment),)


def test_transform_op_tree_applies_the_transformation_at_every_depth() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	tree = [vqubit.X(q0), [vqubit.Y(q1)]]

	transformed = vqubit.transform_op_tree(tree, lambda op: vqubit.Z(*op.qubits))

	assert transformed == [vqubit.Z(q0), [vqubit.Z(q1)]]


def test_transform_op_tree_transforms_inside_a_moment() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	moment = vqubit.Moment([vqubit.X(q0), vqubit.Y(q1)])

	transformed = vqubit.transform_op_tree(moment, lambda op: vqubit.Z(*op.qubits))

	assert transformed == vqubit.Moment([vqubit.Z(q0), vqubit.Z(q1)])


def test_transform_op_tree_can_apply_a_moment_transformation() -> None:
	q0 = vqubit.LineQubit(0)
	moment = vqubit.Moment([vqubit.X(q0)])

	transformed = vqubit.transform_op_tree(
		moment,
		lambda op: vqubit.Z(*op.qubits),
		moment_transformation=lambda transformed_moment: ["moment", transformed_moment],
	)

	assert transformed == ["moment", vqubit.Moment([vqubit.Z(q0)])]


def test_transform_op_tree_returns_a_circuit_for_a_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.X(q0), vqubit.H(q1), vqubit.Z(q0))

	transformed = vqubit.transform_op_tree(circuit, lambda op: vqubit.Z(*op.qubits))

	assert isinstance(transformed, vqubit.Circuit)
	assert list(transformed.all_operations()) == [vqubit.Z(q0), vqubit.Z(q1), vqubit.Z(q0)]
	assert len(transformed) == 2


def test_transform_op_tree_can_transform_each_moment_of_a_circuit() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.X(q1), vqubit.CNOT(q0, q1))

	lengths = vqubit.transform_op_tree(
		circuit,
		lambda op: op,
		moment_transformation=lambda moment: len(moment.operations),
	)

	assert lengths == [2, 1]


def test_transform_op_tree_can_build_an_immutable_tree() -> None:
	q0 = vqubit.LineQubit(0)
	tree = [vqubit.X(q0), vqubit.Y(q0)]

	transformed = vqubit.transform_op_tree(tree, lambda op: op, iter_transformation=tuple)

	assert isinstance(transformed, tuple)
	assert transformed == (vqubit.X(q0), vqubit.Y(q0))


def test_transform_op_tree_returns_none_for_none() -> None:
	assert vqubit.transform_op_tree(None, lambda op: op) is None


def test_transform_op_tree_rejects_non_op_trees() -> None:
	with pytest.raises(TypeError):
		vqubit.transform_op_tree(object(), lambda op: op)


def test_gate_family_matches_instances_types_and_predicates() -> None:
	q0 = vqubit.LineQubit(0)

	instance_family = vqubit.GateFamily(vqubit.X)
	assert instance_family.name == "XPowGate(exponent=1)"
	assert vqubit.X(q0) in instance_family
	assert (vqubit.X**0.5)(q0) not in instance_family

	# A gate *class* matches any instance of that class with ``isinstance``.
	type_family = vqubit.GateFamily(vqubit.XPowGate)
	assert vqubit.X(q0) in type_family
	assert (vqubit.X**0.5)(q0) in type_family
	assert vqubit.Y(q0) not in type_family

	predicate_family = vqubit.GateFamily(predicate=lambda gate: gate.num_qubits == 1)
	assert predicate_family.name == "GateFamily"
	assert vqubit.H(q0) in predicate_family
	assert "H" not in predicate_family


def test_gate_family_requires_a_gate_or_a_predicate() -> None:
	with pytest.raises(ValueError):
		vqubit.GateFamily()


def test_gate_family_validate_describes_the_family() -> None:
	q0 = vqubit.LineQubit(0)
	family = vqubit.GateFamily(vqubit.X, name="pauli-x")

	family.validate(vqubit.X(q0))
	with pytest.raises(ValueError, match="pauli-x"):
		family.validate(vqubit.H(q0))


def test_gateset_membership_and_validation() -> None:
	q0 = vqubit.LineQubit(0)
	gateset = vqubit.Gateset([vqubit.X, vqubit.GateFamily(vqubit.H)], name="clifford")

	assert vqubit.X(q0) in gateset
	assert vqubit.H(q0) in gateset
	assert vqubit.Y(q0) not in gateset
	gateset.validate(vqubit.X(q0))
	with pytest.raises(ValueError, match="clifford"):
		gateset.validate(vqubit.Y(q0))


def test_gateset_with_added_gates_preserves_name_and_description() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gateset = vqubit.Gateset([vqubit.X], name="base", description="one-qubit gates")

	bigger = gateset.with_added_gates(vqubit.CNOT)

	assert gateset.name == "base"
	assert bigger.name == "base"
	assert bigger.description == "one-qubit gates"
	assert vqubit.CNOT(q0, q1) in bigger
	assert vqubit.CNOT(q0, q1) not in gateset
	assert bigger.with_added_gates(vqubit.Y, name="renamed").name == "renamed"


def test_gateset_repr_names_its_gates() -> None:
	gateset = vqubit.Gateset([vqubit.X], name="my-gates")

	assert repr(gateset).startswith("vqubit.Gateset(")
	assert "name='my-gates'" in repr(gateset)
	assert "XPowGate" in repr(gateset)
