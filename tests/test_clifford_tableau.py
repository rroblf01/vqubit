"""Tests for the stabilizer tableau itself (``vqubit.CliffordTableau``).

The tableau is the state the Clifford simulator runs on, so these tests walk its
API directly: the gates it can apply, the generators it reports, measuring and
sampling, and the algebra of composing tableaux.
"""

from __future__ import annotations

import pytest

import vqubit


def test_a_fresh_tableau_is_the_zero_state() -> None:
	tableau = vqubit.CliffordTableau(2)
	assert tableau.n == 2
	assert [str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()] == ["Z(0)", "Z(1)"]
	assert [str(destabilizer) for destabilizer in tableau.destabilizers()] == ["XI", "IX"]
	assert tableau._validate()
	assert tableau.matrix().shape == (4, 4)
	assert "CliffordTableau" in repr(tableau)


def test_apply_paulis_and_hadamard() -> None:
	tableau = vqubit.CliffordTableau(1)
	# X flips |0> to |1>, whose stabilizer is -Z; the gates that follow map it on.
	tableau.apply_x(0)
	assert [str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()] == ["-1*Z(0)"]
	tableau.apply_y(0)
	assert [str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()] == ["Z(0)"]
	tableau.apply_z(0)
	assert [str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()] == ["Z(0)"]
	tableau.apply_h(0)
	assert [str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()] == ["X(0)"]


def test_apply_entangling_gates() -> None:
	# A Bell pair is stabilized by XX and ZZ.
	tableau = vqubit.CliffordTableau(2)
	tableau.apply_h(0)
	tableau.apply_cx(0, 1)
	assert sorted(str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()) == ["X(0)*X(1)", "Z(0)*Z(1)"]
	# H, CZ, H on the second qubit is the Bell preparation again.
	other = vqubit.CliffordTableau(2)
	other.apply_h(0)
	other.apply_cz(0, 1)
	other.apply_h(1)
	assert sorted(str(stabilizer) for stabilizer in other.stabilizer_pauli_strings()) == ["X(0)*X(1)", "X(1)"]


def test_apply_cz_validates_its_exponent() -> None:
	tableau = vqubit.CliffordTableau(2)
	with pytest.raises(ValueError):
		tableau.apply_cz(0, 1, exponent=0.5)


def test_apply_global_phase_changes_nothing() -> None:
	tableau = vqubit.CliffordTableau(1)
	before = tableau.copy()
	tableau.apply_global_phase(1j)
	assert tableau == before


def test_the_tableau_matches_the_clifford_simulator() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	tableau = vqubit.CliffordTableau(2)
	tableau.apply_h(0)
	tableau.apply_cx(0, 1)
	state = vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))).final_state
	expected = [str(stabilizer) for stabilizer in state.stabilizers()]
	assert sorted(str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()) == sorted(expected)
	assert tableau._validate()


def test_measuring_collapses_the_tableau() -> None:
	tableau = vqubit.CliffordTableau(1)
	tableau.apply_h(0)
	outcome = tableau.measure([0], seed=0)
	assert outcome[0] in (0, 1)
	# The state is one of the two basis states now, so Z is a stabilizer again.
	assert [str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()] in (["Z(0)"], ["-1*Z(0)"])


def test_sampling_does_not_change_the_tableau() -> None:
	tableau = vqubit.CliffordTableau(1)
	tableau.apply_h(0)
	before = tableau.copy()
	samples = tableau.sample([0], repetitions=10, seed=1)
	assert len(samples) == 10
	assert all(value in (0, 1) for value in samples)
	assert tableau == before


def test_copy_is_independent() -> None:
	tableau = vqubit.CliffordTableau(1)
	clone = tableau.copy()
	clone.apply_x(0)
	assert clone != tableau
	assert [str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()] == ["Z(0)"]


def test_inverse_undoes_the_tableau() -> None:
	tableau = vqubit.CliffordTableau(2)
	tableau.apply_h(0)
	tableau.apply_cx(0, 1)
	tableau.apply_x(1)
	assert tableau.then(tableau.inverse()).matrix().tolist() == vqubit.CliffordTableau(2).matrix().tolist()


def test_then_composes_tableaux() -> None:
	first = vqubit.CliffordTableau(2)
	first.apply_h(0)
	second = vqubit.CliffordTableau(2)
	second.apply_cx(0, 1)
	composed = first.then(second)
	# The composition prepares the Bell state, the same one the two gates do.
	assert sorted(str(stabilizer) for stabilizer in composed.stabilizer_pauli_strings()) == [
		"X(0)*X(1)",
		"Z(0)*Z(1)",
	]
	with pytest.raises(ValueError):
		first.then(vqubit.CliffordTableau(3))
	with pytest.raises(TypeError):
		first.then("not a tableau")  # ty: ignore[invalid-argument-type]


def test_from_flat_reads_the_native_layout() -> None:
	import vqubit_native

	# A fresh native tableau is the |00> state, and ``from_flat`` reads that.
	native = vqubit_native.NativeTableau(2)
	rebuilt = vqubit.CliffordTableau.from_flat(native.stabilizer_tableau(), 2)
	assert [str(stabilizer) for stabilizer in rebuilt.stabilizer_pauli_strings()] == ["Z(0)", "Z(1)"]
	assert rebuilt._validate()


def test_tableau_equality_and_hash() -> None:
	first = vqubit.CliffordTableau(1)
	second = vqubit.CliffordTableau(1)
	assert first == second
	assert hash(first) == hash(second)
	assert first != vqubit.CliffordTableau(2)
	assert first != "not a tableau"
	first.apply_h(0)
	assert first != second
