"""Tests for the arithmetic of Pauli strings and the gate base classes.

The Pauli types support the algebra Cirq gives them (addition builds a sum,
multiplication composes, division scales), and the gate base classes define a
contract every concrete gate follows: an abstract piece raises
`NotImplementedError`, the symbolic forms have no unitary, and the diagram
suffixes carry the exponent.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import vqubit
from tests.helpers import assert_matrix_close

P = vqubit.Pauli


# --- Pauli algebra -----------------------------------------------------------


def test_pauli_string_arithmetic() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	first = vqubit.PauliString({q0: P.X})
	second = vqubit.PauliString({q1: P.Z})
	# X + Z is a sum of two terms, and so is the difference.
	assert isinstance(first + second, vqubit.PauliSum)
	assert len(first + second) == 2
	assert len(first - second) == 2
	# Adding a scalar makes a sum with the identity term.
	summed = first + 0.5
	assert isinstance(summed, vqubit.PauliSum)
	assert len(summed) == 2
	assert summed[vqubit.PauliString()] == 0.5
	# Division and multiplication scale the coefficient.
	assert (first / 2).coefficient == 0.5
	assert (2 * first).coefficient == 2
	assert (-first).coefficient == -1
	assert (first * 2j).coefficient == 2j
	assert first * second == vqubit.PauliString({q0: P.X, q1: P.Z})
	# The product of X with itself is the identity.
	assert first * first == vqubit.PauliString()


def test_pauli_string_expectation_values() -> None:
	q0 = vqubit.LineQubit(0)
	string = vqubit.PauliString({q0: P.Z})
	assert string.expectation_from_state_vector([1, 0], qubit_map={q0: 0}) == pytest.approx(1)
	assert string.expectation_from_state_vector([0, 1], qubit_map={q0: 0}) == pytest.approx(-1)
	# The density-matrix form agrees, and a scaling coefficient multiplies it.
	zero_state = [[1, 0], [0, 0]]
	assert string.expectation_from_density_matrix(zero_state, qubit_map={q0: 0}) == pytest.approx(1)
	assert (2 * string).expectation_from_state_vector([1, 0], qubit_map={q0: 0}) == pytest.approx(2)
	with pytest.raises(ValueError, match="power of two"):
		string.expectation_from_state_vector([1, 0, 0], qubit_map={q0: 0})


def test_pauli_sum_arithmetic_and_matrix() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	x = vqubit.PauliString({q0: P.X})
	z = vqubit.PauliString({q1: P.Z})
	total = x + z
	assert total[x] == 1
	assert len(total) == 2
	assert (total / 2).items() is not None
	assert (2 * total)[x] == 2
	assert (-total)[x] == -1
	assert (total - z)[x] == 1
	assert isinstance(total * total, vqubit.PauliSum)
	assert isinstance(total * x, vqubit.PauliSum)
	# X*X and Z*Z are identities and the cross terms cancel, so two terms remain.
	assert len(total * total) == 2
	# (X + Z) ** 2 has the four products, and the diagonal ones are identities.
	combined = total * total
	assert combined[vqubit.PauliString()] == 2
	assert len(combined) == 2  # the identity and X*Z; the cross terms cancel
	# wrap and from_pauli_strings build the same sum.
	assert vqubit.PauliSum.wrap(x) == vqubit.PauliSum([(x, 1)])
	assert vqubit.PauliSum.from_pauli_strings([x, z]) == total


def test_pauli_sum_expectation_values() -> None:
	q0 = vqubit.LineQubit(0)
	total = 0.5 * vqubit.PauliString({q0: P.Z}) + 0.25
	# <0| 0.5 Z + 0.25 |0> = 0.75.
	assert total.expectation_from_state_vector([1, 0], qubit_map={q0: 0}) == pytest.approx(0.75)
	assert total.expectation_from_density_matrix([[1, 0], [0, 0]], qubit_map={q0: 0}) == pytest.approx(0.75)


def test_dense_pauli_string_arithmetic() -> None:
	first = vqubit.DensePauliString("XY")
	second = vqubit.DensePauliString("ZY")
	# Elementwise: X*Z = -iY and Y*Y = I, so the coefficient collects the phase.
	assert (first * second) == vqubit.DensePauliString("YI", -1j)
	assert (2 * first) == vqubit.DensePauliString("XY", 2)
	assert (first * 2) == vqubit.DensePauliString("XY", 2)
	assert len(first) == 2
	assert first[0] is P.X
	assert first.to_pauli_string(vqubit.LineQubit.range(2)) == vqubit.PauliString(
		{vqubit.LineQubit(0): P.X, vqubit.LineQubit(1): P.Y}
	)
	assert_matrix_close(first.matrix(), np.kron(vqubit.unitary(vqubit.X), vqubit.unitary(vqubit.Y)))


def test_pauli_enum() -> None:
	q0 = vqubit.LineQubit(0)
	assert str(P.X) == "X"
	assert_matrix_close(P.Z.matrix(), [[1, 0], [0, -1]])
	assert P.Y.on(q0) == vqubit.Y(q0)
	# All four Paulis are distinct values with distinct matrices.
	matrices = {P.I, P.X, P.Y, P.Z}
	assert len(matrices) == 4


# --- gate base contracts -----------------------------------------------------


class _DummyEigenGate(vqubit.EigenGate):
	"""A minimal gate that only provides its eigencomponents."""

	def _num_qubits_(self) -> int:
		return 1

	def _eigen_components(self) -> list[tuple[float, list[list[complex]]]]:
		return [(1.0, [[1, 0], [0, 0]]), (-1.0, [[0, 0], [0, 1]])]


def test_eigen_gate_base_class() -> None:
	# The base class cannot supply eigencomponents, so asking for them fails...
	with pytest.raises(NotImplementedError):
		vqubit.EigenGate()._eigen_components()
	# ... and a subclass only provides them: the base builds the unitary from
	# them, with the exponent and the global phase.
	gate = _DummyEigenGate(exponent=0.5)
	assert_matrix_close(vqubit.unitary(gate), vqubit.unitary(vqubit.Z**0.5))
	shifted = _DummyEigenGate(exponent=1, global_shift=0.5)
	assert_matrix_close(vqubit.unitary(shifted), [[1j, 0], [0, -1j]])
	# A parameterized exponent has no unitary.
	parameterized = _DummyEigenGate(exponent=vqubit.Parameter("t"))
	assert vqubit.is_parameterized(parameterized)
	with pytest.raises(TypeError, match="has no unitary"):
		vqubit.unitary(parameterized)


def test_rotation_gate_base_class_has_no_matrix_of_its_own() -> None:
	# A rotation supplies its matrix; the base class asks for it.
	from vqubit._gates import _RotationGate

	with pytest.raises(NotImplementedError):
		_RotationGate(0.5)._rotation_matrix(0.5)


def test_fixed_gate_base_class_contract() -> None:
	# CSWAP is the fixed gate that ships; its base class asks for a matrix and
	# diagram symbols, which the concrete gate supplies.
	assert vqubit.CSWAP.num_qubits == 3
	assert str(vqubit.CSWAP) == "CSWAPGate"
	assert repr(vqubit.CSWAP) == "vqubit.CSWAPGate()"
	assert vqubit.CSWAP._inverse_() == vqubit.CSWAP


def test_power_gate_strings_and_symbolic_suffixes() -> None:
	t = vqubit.Parameter("t")
	q0 = vqubit.LineQubit(0)
	# A symbolic exponent is drawn as `^t` rather than dropped.
	diagram = str(vqubit.Circuit((vqubit.X**t).on(q0)))  # ty: ignore[unsupported-operator]
	assert "X^t" in diagram
	assert repr(vqubit.CZPowGate(exponent=0.5)) == "vqubit.CZPowGate(exponent=0.5, global_shift=0.0)"
	assert str(vqubit.SwapPowGate(exponent=0.5)) == "SwapPowGate(exponent=0.5)"
	assert_matrix_close(vqubit.unitary(vqubit.H**0.5), vqubit.unitary(vqubit.HPowGate(exponent=0.5)))
	assert_matrix_close(vqubit.unitary(vqubit.H**-0.5), vqubit.unitary(vqubit.HPowGate(exponent=-0.5)))


def test_identity_gate_forms() -> None:
	assert vqubit.IdentityGate(2).num_qubits == 2
	assert vqubit.I**0.5 == vqubit.I
	assert repr(vqubit.IdentityGate(3)) == "vqubit.IdentityGate(3)"
	assert str(vqubit.IdentityGate(3)) == "IdentityGate"
	assert_matrix_close(vqubit.unitary(vqubit.IdentityGate(2)), np.eye(4))


def test_measurement_gate_forms() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	gate = vqubit.MeasurementGate("m", num_qubits=2, invert_mask=(True, False))
	assert gate.full_invert_mask == (True, False)
	assert "m" in repr(gate)
	with pytest.raises(ValueError, match="invert_mask"):
		vqubit.MeasurementGate("m", num_qubits=1, invert_mask=(True, True))
	# measure() flattens an iterable of qubits, and needs at least one.
	operation = vqubit.measure([q0, q1], key="both")
	assert operation.qubits == (q0, q1)
	with pytest.raises(ValueError, match="at least one qubit"):
		vqubit.measure()
	# measure_each keys by str by default and by the given function otherwise.
	assert [vqubit.measurement_key_name(op) for op in vqubit.measure_each(q0, q1)] == [str(q0), str(q1)]
	assert [vqubit.measurement_key_name(op) for op in vqubit.measure_each(q0, q1, key_func=lambda q: "k")] == ["k", "k"]


def test_sqrt_iswap_can_be_built_explicitly() -> None:
	gate = vqubit.SqrtISwapGate()
	assert isinstance(gate, vqubit.ISwapPowGate)
	assert_matrix_close(vqubit.unitary(gate), vqubit.unitary(vqubit.SQRT_ISWAP))
	assert vqubit.SQRT_ISWAP == vqubit.ISWAP**0.5


def test_phased_gates_report_their_symbols_and_resolve() -> None:
	t = vqubit.Parameter("t")
	q0 = vqubit.LineQubit(0)
	assert vqubit.PhasedXZGate()._diagram_symbols_([q0]) == ("PhXZ",)
	assert "PhXZ" in str(vqubit.Circuit(vqubit.PhasedXZGate(x_exponent=0.5).on(q0)))
	# A concrete gate resolves to itself; a symbolic one takes the value.
	concrete = vqubit.PhasedXPowGate(phase_exponent=0.25)
	assert vqubit.resolve_parameters(concrete, vqubit.ParamResolver({})) is concrete
	symbolic = vqubit.PhasedXPowGate(phase_exponent=t)
	assert vqubit.resolve_parameters(symbolic, vqubit.ParamResolver({"t": 0.25})) == concrete
	assert vqubit.is_parameterized(vqubit.PhasedXZGate(z_exponent=t))
	assert (
		vqubit.resolve_parameters(vqubit.PhasedXZGate(z_exponent=t), vqubit.ParamResolver({"t": 0}))
		== vqubit.PhasedXZGate()
	)


def test_pauli_product_powers() -> None:
	t = vqubit.Parameter("t")
	assert vqubit.XXPowGate(exponent=0.5).exponent == 0.5
	assert vqubit.resolve_parameters(
		vqubit.XXPowGate(exponent=t), vqubit.ParamResolver({"t": 0.5})
	) == vqubit.XXPowGate(exponent=0.5)
	assert vqubit.XXPowGate(exponent=0.5)._inverse_() == vqubit.XXPowGate(exponent=-0.5)
	# The symbolic form has no unitary, and the concrete one is XX ** e.
	assert vqubit.is_parameterized(vqubit.YYPowGate(exponent=t))
	assert_matrix_close(vqubit.unitary(vqubit.ZZPowGate(exponent=0.5)), vqubit.unitary(vqubit.ZZ**0.5))


def test_controlled_gate_forms() -> None:
	t = vqubit.Parameter("t")
	symbolic = vqubit.ControlledGate(vqubit.ZPowGate(exponent=t), num_controls=1)
	with pytest.raises(TypeError, match="has no unitary"):
		vqubit.unitary(symbolic)
	assert "ControlledGate(ZPowGate" in str(symbolic)
	concrete = vqubit.ControlledGate(vqubit.Z, num_controls=1)
	assert vqubit.resolve_parameters(concrete, vqubit.ParamResolver({})) is concrete
	# A gate with no inverse makes the controlled gate have none either.
	matrix_gate = vqubit.MatrixGate([[1, 0], [0, 1]])
	assert vqubit.ControlledGate(matrix_gate, num_controls=1)._inverse_() is None


def test_global_phase_gate_forms() -> None:
	t = vqubit.Parameter("t")
	assert vqubit.GlobalPhaseGate(0.5)._diagram_symbols_([]) == ()
	assert str(vqubit.global_phase(0.5)) == "GlobalPhaseGate(exponent=0.5)"
	assert repr(vqubit.global_phase(0.5)) == "vqubit.GlobalPhaseGate(0.5)"
	concrete = vqubit.global_phase(0.5)
	assert vqubit.resolve_parameters(concrete, vqubit.ParamResolver({})) is concrete
	assert vqubit.resolve_parameters(vqubit.global_phase(t), vqubit.ParamResolver({"t": 0.5})) == concrete


def test_wait_gate_forms() -> None:
	gate = vqubit.WaitGate(5.0, 2)
	assert gate.num_qubits == 2
	assert repr(gate) == "vqubit.WaitGate(5.0, num_qubits=2)"
	assert gate == vqubit.WaitGate(5.0, 2)
	assert gate != vqubit.WaitGate(6.0, 2)


def test_ms_and_fsim_gates() -> None:
	t = vqubit.Parameter("t")
	ms = vqubit.MSGate(rads=0.5)
	assert ms._diagram_symbols_([vqubit.LineQubit(0), vqubit.LineQubit(1)]) == ("MS", "MS")
	assert vqubit.resolve_parameters(ms, vqubit.ParamResolver({})) is ms
	assert vqubit.is_parameterized(vqubit.MSGate(rads=t))
	assert_matrix_close(vqubit.unitary(vqubit.MSGate(rads=0)), np.eye(4))
	fsim = vqubit.FSimGate(theta=0.1, phi=0.2)
	assert fsim._diagram_symbols_([vqubit.LineQubit(0), vqubit.LineQubit(1)]) == ("fSim", "fSim")
	assert vqubit.resolve_parameters(fsim, vqubit.ParamResolver({})) is fsim
	assert vqubit.is_parameterized(vqubit.FSimGate(theta=t, phi=0))
	assert_matrix_close(vqubit.unitary(vqubit.FSimGate(theta=0, phi=0)), np.eye(4))
	assert math.isclose(abs(vqubit.unitary(fsim)[0][0]), 1.0)
