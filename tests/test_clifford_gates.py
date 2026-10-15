"""Tests for the Clifford gates."""

from __future__ import annotations

import vqubit
from tests.helpers import up_to_phase


def test_single_qubit_clifford_group_size() -> None:
	assert len(vqubit.SingleQubitCliffordGate.all_single_qubit_cliffords) == 24


def test_single_qubit_clifford_known_unitaries() -> None:
	assert up_to_phase(vqubit.SingleQubitCliffordGate.H._unitary_(), vqubit.unitary(vqubit.H))
	assert up_to_phase(vqubit.SingleQubitCliffordGate.S._unitary_(), vqubit.unitary(vqubit.S))
	assert up_to_phase(vqubit.SingleQubitCliffordGate.X._unitary_(), vqubit.unitary(vqubit.X))
	assert up_to_phase(vqubit.SingleQubitCliffordGate.Y._unitary_(), vqubit.unitary(vqubit.Y))
	assert up_to_phase(vqubit.SingleQubitCliffordGate.Z._unitary_(), vqubit.unitary(vqubit.Z))


def test_single_qubit_clifford_from_unitary() -> None:
	assert vqubit.SingleQubitCliffordGate.from_unitary(vqubit.unitary(vqubit.H)) == vqubit.SingleQubitCliffordGate.H
	assert vqubit.SingleQubitCliffordGate.from_unitary(vqubit.unitary(vqubit.T)) is None


def test_single_qubit_clifford_arithmetic() -> None:
	sqrt = vqubit.SingleQubitCliffordGate.X_sqrt
	assert sqrt.merged_with(sqrt) == vqubit.SingleQubitCliffordGate.X
	assert (vqubit.SingleQubitCliffordGate.S**4) == vqubit.SingleQubitCliffordGate.I
	assert (
		vqubit.SingleQubitCliffordGate.Z.merged_with(vqubit.SingleQubitCliffordGate.Z)
		== vqubit.SingleQubitCliffordGate.I
	)


def test_single_qubit_clifford_pauli_tuple() -> None:
	gate = vqubit.SingleQubitCliffordGate.H
	assert gate.pauli_tuple(vqubit.Pauli.X) == (vqubit.Pauli.Z, False)
	assert gate.pauli_tuple(vqubit.Pauli.Z) == (vqubit.Pauli.X, False)


def test_single_qubit_clifford_commutes() -> None:
	assert vqubit.SingleQubitCliffordGate.S.commutes_with_pauli(vqubit.Pauli.Z)
	assert not vqubit.SingleQubitCliffordGate.S.commutes_with_pauli(vqubit.Pauli.X)


def test_clifford_gate_swap_variants() -> None:
	for gate in (vqubit.CXSWAP, vqubit.CZSWAP):
		assert gate.num_qubits == 2
		assert vqubit.is_unitary(vqubit.unitary(gate))
		assert (gate**2) == (gate**2)
		assert (gate**-1) == (gate**-1)


def test_clifford_gate_from_op_list() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gate = vqubit.CliffordGate.from_op_list([vqubit.CNOT(q0, q1), vqubit.SWAP(q0, q1)], [q0, q1])
	assert up_to_phase(vqubit.unitary(gate), vqubit.unitary(vqubit.CXSWAP))
