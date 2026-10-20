"""Tests for the half-integer powers and the remaining tableau plumbing.

The tableau can apply any Clifford power of X, Y, Z and H (``S = Z ** 0.5`` and
``S ** 3``, the quarter turns of Cirq), measure and sample without mutating
itself, split into factors, and report itself as a stabilizer state. These are
the paths a circuit with ``S``, ``S ** 3`` or ``sqrt(X)`` goes through.
"""

from __future__ import annotations

import numpy as np
import pytest

import vqubit


def stabilizers(tableau: vqubit.CliffordTableau) -> list[str]:
	return [str(stabilizer) for stabilizer in tableau.stabilizer_pauli_strings()]


def test_sqrt_x_twice_is_x() -> None:
	# Applying the square root twice has to land on the gate itself, whatever
	# the internal representation does in between.
	twice = vqubit.CliffordTableau(1)
	twice.apply_x(0, 0.5)
	twice.apply_x(0, 0.5)
	once = vqubit.CliffordTableau(1)
	once.apply_x(0, 1)
	assert stabilizers(twice) == stabilizers(once) == ["-1*Z(0)"]


def test_sqrt_z_is_the_s_gate() -> None:
	# On |+>, S rotates X into Y and S ** 3 rotates it into -Y.
	s = vqubit.CliffordTableau(1)
	s.apply_h(0)
	s.apply_z(0, 0.5)
	assert stabilizers(s) == ["Y(0)"]
	s_dagger = vqubit.CliffordTableau(1)
	s_dagger.apply_h(0)
	s_dagger.apply_z(0, 1.5)
	assert stabilizers(s_dagger) == ["-1*Y(0)"]


def test_quarter_turns_of_y_and_h() -> None:
	# The same e^{iπ/4 P} construction for the other Paulis: two quarter turns
	# make the full gate.
	y_quarter = vqubit.CliffordTableau(1)
	y_quarter.apply_y(0, 0.5)
	y_quarter.apply_y(0, 0.5)
	y_full = vqubit.CliffordTableau(1)
	y_full.apply_y(0, 1)
	assert stabilizers(y_quarter) == stabilizers(y_full)
	# H is an involution, and ``apply_h`` only accepts integer powers.
	h_twice = vqubit.CliffordTableau(1)
	h_twice.apply_h(0)
	h_twice.apply_h(0)
	assert stabilizers(h_twice) == ["Z(0)"]


def test_even_powers_are_the_identity() -> None:
	for axis in (0,):
		tableau = vqubit.CliffordTableau(1)
		tableau.apply_h(axis)
		before = stabilizers(tableau)
		tableau.apply_x(axis, 2)
		tableau.apply_y(axis, 4)
		tableau.apply_z(axis, -2)
		assert stabilizers(tableau) == before


def test_global_shift_is_accepted_and_ignored() -> None:
	# A global phase does not change a stabilizer state.
	shifted = vqubit.CliffordTableau(1)
	shifted.apply_z(0, 1, global_shift=0.5)
	plain = vqubit.CliffordTableau(1)
	plain.apply_z(0, 1)
	assert shifted == plain


def test_non_clifford_powers_are_refused() -> None:
	for gate, exponent in (("x", 0.25), ("y", 0.3), ("z", 0.25)):
		tableau = vqubit.CliffordTableau(1)
		with pytest.raises(ValueError, match=f"{gate.upper()} exponent must be half integer"):
			getattr(tableau, f"apply_{gate}")(0, exponent)
	with pytest.raises(ValueError, match="H exponent must be integer"):
		vqubit.CliffordTableau(1).apply_h(0, 0.5)
	with pytest.raises(ValueError, match="CX exponent must be integer"):
		vqubit.CliffordTableau(2).apply_cx(0, 1, 0.5)
	with pytest.raises(ValueError, match="CZ exponent must be integer"):
		vqubit.CliffordTableau(2).apply_cz(0, 1, 0.5)


def test_cx_and_cz_with_an_even_exponent_do_nothing() -> None:
	tableau = vqubit.CliffordTableau(2)
	before = tableau.copy()
	tableau.apply_cx(0, 1, 2)
	tableau.apply_cz(0, 1, 4)
	assert tableau == before


def test_str_draws_the_rows() -> None:
	tableau = vqubit.CliffordTableau(2)
	tableau.apply_h(0)
	tableau.apply_cz(0, 1)
	text = str(tableau)
	assert isinstance(text, str)
	# One line per stabilizer, each starting with the sign.
	lines = [line for line in text.splitlines() if line.strip()]
	assert len(lines) == 2
	assert all(line[0] in "+-" for line in lines)
	assert {"X", "Z"} <= {character for line in lines for character in line}


def test_matrix_product_composes_in_circuit_order() -> None:
	first = vqubit.CliffordTableau(2)
	first.apply_h(0)
	second = vqubit.CliffordTableau(2)
	second.apply_cx(0, 1)
	# ``second @ first`` is "first, then second", the same as ``first.then``.
	assert second @ first == first.then(second)
	with pytest.raises(TypeError):
		first @ "not a tableau"  # ty: ignore[unsupported-operator]


def test_measuring_an_out_of_range_axis_fails() -> None:
	with pytest.raises(IndexError):
		vqubit.CliffordTableau(1).measure([3])


def test_measuring_a_deterministic_axis_ignores_the_seed() -> None:
	tableau = vqubit.CliffordTableau(2)
	outcomes = [tableau.measure([0, 1], seed=seed) for seed in range(5)]
	assert outcomes == [[0, 0]] * 5


def test_sampling_several_axes_returns_a_matrix() -> None:
	tableau = vqubit.CliffordTableau(2)
	tableau.apply_h(0)
	samples = tableau.sample([0, 1], repetitions=20, seed=2)
	assert isinstance(samples, np.ndarray)
	assert samples.shape == (20, 2)
	# The second axis is deterministic, the first is not.
	assert set(samples[:, 1].tolist()) == {0}
	assert set(samples[:, 0].tolist()) <= {0, 1}


def test_factor_validates_its_axes() -> None:
	tableau = vqubit.CliffordTableau(2)
	with pytest.raises(ValueError, match="not distinct"):
		tableau.factor([0, 0])
	with pytest.raises(ValueError, match="not distinct"):
		vqubit.CliffordTableau(2).factor([0, 2])
	extracted, remainder = tableau.factor([0])
	assert extracted.n == 1
	assert remainder.n == 1


def test_kron_rejects_anything_but_a_tableau() -> None:
	with pytest.raises(TypeError, match="kron needs a CliffordTableau"):
		vqubit.CliffordTableau(1).kron(3)


def test_from_flat_validates_the_layout() -> None:
	with pytest.raises(ValueError, match="needs 20 entries"):
		vqubit.CliffordTableau.from_flat([0, 1, 2], 2)


def test_a_clifford_state_reports_itself() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	state = vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))).final_state
	assert state.qubits == (q0, q1)
	assert "2 qubits" in repr(state)
	assert [str(stabilizer) for stabilizer in state.stabilizers()] == ["X(0)*X(1)", "Z(0)*Z(1)"]
	assert state.dirac_notation() == "0.707|00⟩ + 0.707|11⟩"


def test_the_clifford_simulator_refuses_what_it_cannot_run() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	with pytest.raises(ValueError, match="cannot run a circuit with measurements"):
		vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.measure(q0)))
	with pytest.raises(TypeError, match="doesn't support"):
		vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.T(q0)))
	# A qubit order that leaves out a used qubit is an error, not a surprise.
	with pytest.raises(ValueError, match="qubit_order is missing"):
		vqubit.CliffordSimulator().simulate(vqubit.Circuit(vqubit.CNOT(q0, q1)), qubit_order=[q0])
