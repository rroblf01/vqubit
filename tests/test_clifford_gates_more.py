"""Tests for the single-qubit Clifford group and the multi-qubit Clifford gate.

The group has 24 elements. These tests walk the constructor helpers, check that
every element decomposes into Pauli rotations that rebuild its unitary, and pin
down the algebra (merge, conjugate, commute, powers).
"""

from __future__ import annotations

import cmath
import itertools

import pytest

import vqubit
from tests.helpers import assert_matrix_close, up_to_phase

np = pytest.importorskip("numpy", exc_type=ImportError)

Gate = vqubit.SingleQubitCliffordGate
P = vqubit.Pauli


def all_cliffords() -> tuple[Gate, ...]:
	return tuple(Gate.all_single_qubit_cliffords)


def test_the_group_has_twenty_four_elements() -> None:
	cliffords = all_cliffords()
	assert len(cliffords) == 24
	assert len(set(cliffords)) == 24


def test_from_xz_map_rejects_an_impossible_map() -> None:
	# X and Z must be mapped to distinct Paulis; mapping both to X is not an
	# automorphism of the Pauli group.
	with pytest.raises(ValueError, match="No single-qubit Clifford"):
		Gate.from_xz_map((P.X, False), (P.X, False))


def test_from_single_map_infers_the_other_two_images() -> None:
	# "X goes to Z" determines the gate: the rest of the map follows from the
	# Pauli algebra, so the result is one of the 24 elements.
	gate = Gate.from_single_map(x_to=(P.Z, False))
	assert gate in all_cliffords()
	assert gate.pauli_tuple(P.X) == (P.Z, False)
	# The keyword form and the mapping form agree.
	assert Gate.from_single_map({P.X: (P.Z, False)}) == gate
	with pytest.raises(ValueError, match="exactly one transform"):
		Gate.from_single_map()
	with pytest.raises(ValueError, match="exactly one transform"):
		Gate.from_single_map(x_to=(P.Z, False), y_to=(P.Y, False))


def test_from_double_map_infers_the_third_image() -> None:
	gate = Gate.from_double_map(x_to=(P.Z, False), z_to=(P.X, False))
	assert gate == Gate.H
	with pytest.raises(ValueError, match="exactly two transforms"):
		Gate.from_double_map(x_to=(P.Z, False))
	with pytest.raises(ValueError, match="exactly two transforms"):
		Gate.from_double_map(x_to=(P.Z, False), y_to=(P.Y, False), z_to=(P.X, False))


def test_from_pauli_and_quarter_turns() -> None:
	assert Gate.from_pauli(P.X) == Gate.X
	assert Gate.from_pauli(P.Y) == Gate.Y
	assert Gate.from_pauli(P.Z) == Gate.Z
	# The square root of X is the quarter turn, and twice it is X.
	assert Gate.from_pauli(P.X, sqrt=True).merged_with(Gate.from_pauli(P.X, sqrt=True)) == Gate.X
	for pauli in (P.X, P.Y, P.Z):
		assert Gate.from_quarter_turns(pauli, 0) == Gate.I
		assert Gate.from_quarter_turns(pauli, 2) == Gate.from_pauli(pauli)
		assert Gate.from_quarter_turns(pauli, 4) == Gate.I
		assert Gate.from_quarter_turns(pauli, -1) == Gate.from_pauli(pauli, sqrt=True) ** -1
		assert Gate.from_quarter_turns(pauli, 3) == Gate.from_pauli(pauli, sqrt=True) ** -1


def test_from_unitary_finds_the_gate_or_refuses() -> None:
	for gate in all_cliffords():
		found = Gate.from_unitary(vqubit.unitary(gate))
		assert found == gate, gate
	# The T gate is not in the group, and a 3x3 matrix is not a qubit gate.
	assert Gate.from_unitary(vqubit.unitary(vqubit.T)) is None
	assert Gate.from_unitary([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) is None


def test_from_unitary_with_global_phase_returns_the_phase() -> None:
	phase = cmath.exp(0.7j)
	matrix = [[phase * entry for entry in row] for row in vqubit.unitary(Gate.H)]
	found = Gate.from_unitary_with_global_phase(matrix)
	assert found is not None
	gate, recovered = found
	assert gate == Gate.H
	assert recovered == pytest.approx(phase)
	assert Gate.from_unitary_with_global_phase(vqubit.unitary(vqubit.T)) is None


def test_pauli_tuple_of_every_element() -> None:
	for gate in all_cliffords():
		for pauli in (P.X, P.Y, P.Z):
			image, _flip = gate.pauli_tuple(pauli)
			assert image in (P.X, P.Y, P.Z)


def test_decompose_rotation_rebuilds_the_unitary() -> None:
	# The decomposition is a list of (Pauli, quarter turns); turning it into
	# operations has to give the gate back, up to the global phase the algebra
	# drops.
	for gate in all_cliffords():
		rotation = gate.decompose_rotation()
		if not rotation:
			assert up_to_phase(np.asarray(vqubit.unitary(gate)), np.eye(2))
			continue
		q = vqubit.LineQubit(0)
		operations = [
			(
				vqubit.X ** (turns / 2)
				if pauli is P.X
				else vqubit.Y ** (turns / 2)
				if pauli is P.Y
				else vqubit.Z ** (turns / 2)
			).on(q)
			for pauli, turns in rotation
		]
		circuit = vqubit.Circuit(operations)
		assert up_to_phase(np.asarray(vqubit.unitary(circuit)), np.asarray(vqubit.unitary(gate))), (gate, rotation)


def test_decompose_gate_of_hadamard_is_just_h() -> None:
	assert Gate.H.decompose_gate() == [vqubit.H]
	# Every other element decomposes into H and Pauli rotations whose product
	# is the gate.
	for gate in all_cliffords():
		if gate == Gate.H:
			continue
		operations = [gate_part.on(vqubit.LineQubit(0)) for gate_part in gate.decompose_gate()]
		if not operations:
			assert up_to_phase(np.asarray(vqubit.unitary(gate)), np.eye(2))
			continue
		assert up_to_phase(np.asarray(vqubit.unitary(vqubit.Circuit(operations))), np.asarray(vqubit.unitary(gate))), (
			gate
		)


def test_merging_two_gates_is_their_product() -> None:
	for first, second in itertools.islice(itertools.product(all_cliffords(), repeat=2), 0, 576, 37):
		merged = first.merged_with(second)
		expected = np.asarray(vqubit.unitary(second)) @ np.asarray(vqubit.unitary(first))
		assert up_to_phase(np.asarray(vqubit.unitary(merged)), expected), (first, second)


def test_equivalent_gate_before_moves_a_gate_across() -> None:
	for first in all_cliffords():
		conjugated = first.equivalent_gate_before(Gate.H)
		# first ** -1 * H * first, i.e. H conjugated by the gate.
		expected = (
			np.asarray(vqubit.unitary(first)).conj().T
			@ np.asarray(vqubit.unitary(Gate.H))
			@ np.asarray(vqubit.unitary(first))
		)
		assert up_to_phase(np.asarray(vqubit.unitary(conjugated)), expected)


def test_commutation() -> None:
	# Z commutes with Z and S, but not with H or X.
	assert Gate.Z.commutes_with_pauli(P.Z)
	assert Gate.S.commutes_with_pauli(P.Z)
	assert not Gate.H.commutes_with_pauli(P.Z)
	assert not Gate.X.commutes_with_pauli(P.Z)
	assert Gate.H.commutes_with_single_qubit_gate(Gate.H)
	assert not Gate.S.commutes_with_single_qubit_gate(Gate.H)
	assert Gate.I.commutes_with_single_qubit_gate(Gate.H)


def test_integer_powers_follow_the_order_of_the_element() -> None:
	# The group has elements of order 1, 2, 3 and 4, so the power is taken
	# modulo the element's own order (modulo twelve always works).
	orders = set()
	for gate in all_cliffords():
		order = 1
		current = gate
		while current != Gate.I:
			current = current.merged_with(gate)
			order += 1
		orders.add(order)
		assert gate**0 == Gate.I
		assert gate**12 == Gate.I
		assert gate**-12 == Gate.I
		assert gate**2 == gate.merged_with(gate)
		assert gate.merged_with(gate**-1) == Gate.I
		assert gate**3 == gate.merged_with(gate).merged_with(gate)
	assert orders == {1, 2, 3, 4}


def test_half_powers_where_they_exist() -> None:
	# Every element of the group has a square root inside the group; the
	# quarter turns expose two of them.
	assert Gate.X**0.5 in all_cliffords()
	assert (Gate.X**0.5).merged_with(Gate.X**0.5) == Gate.X
	assert Gate.I**0.5 == Gate.I
	# A twelfth of a turn is not in the group.
	with pytest.raises(NotImplementedError, match="cannot raise"):
		Gate.X**0.25


def test_unitary_and_unitary_flag() -> None:
	assert vqubit.has_unitary(Gate.H)
	assert_matrix_close(vqubit.unitary(Gate.H), np.asarray(vqubit.unitary(vqubit.H)))
	assert "H" in str(Gate.H) or "SingleQubitCliffordGate" in str(Gate.H)
	assert "SingleQubitCliffordGate" in repr(Gate.H)


def test_the_multi_qubit_clifford_gate() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	swap_like = vqubit.CliffordGate.from_unitary(vqubit.unitary(vqubit.CNOT))
	assert swap_like.num_qubits == 2
	assert vqubit.has_stabilizer_effect(swap_like)
	assert_matrix_close(vqubit.unitary(swap_like), vqubit.unitary(vqubit.CNOT))
	# Powers multiply the unitary, and the inverse undoes it.
	assert_matrix_close(np.asarray(vqubit.unitary(swap_like**-1)) @ np.asarray(vqubit.unitary(swap_like)), np.eye(4))
	assert swap_like**2 == vqubit.CliffordGate.from_unitary(np.eye(4).tolist())
	with pytest.raises(NotImplementedError, match="integer powers"):
		swap_like**0.5
	assert "2q" in str(swap_like)
	assert "CliffordGate" in repr(swap_like)
	# The same gate, composed from its operations.
	from_ops = vqubit.CliffordGate.from_op_list([vqubit.CNOT(q0, q1)], [q0, q1])
	assert from_ops == swap_like
