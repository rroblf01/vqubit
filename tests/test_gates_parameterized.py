"""Tests for the gates that carry parameters, matrices or controls.

A gate is allowed to be symbolic: it then reports that it is parameterized, has
no unitary, and turns into the concrete gate once a resolver supplies a value.
These tests pin that contract down for every gate family that supports it, and
the validation the hand-written gates do on their arguments.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

import vqubit
from tests.helpers import assert_matrix_close


def symbolic_gates() -> list[tuple[str, Any, Any]]:
	"""(name, gate built on ``Parameter('t')``, gate with ``t = 0.5``)."""
	t = vqubit.Parameter("t")
	return [
		("x", vqubit.XPowGate(exponent=t), vqubit.XPowGate(exponent=0.5)),
		("cz", vqubit.CZPowGate(exponent=t), vqubit.CZPowGate(exponent=0.5)),
		("ccx", vqubit.CCXPowGate(exponent=t), vqubit.CCXPowGate(exponent=0.5)),
		("swap", vqubit.SwapPowGate(exponent=t), vqubit.SwapPowGate(exponent=0.5)),
		("iswap", vqubit.ISwapPowGate(exponent=t), vqubit.ISwapPowGate(exponent=0.5)),
		("xx", vqubit.XXPowGate(exponent=t), vqubit.XXPowGate(exponent=0.5)),
		("phased_x", vqubit.PhasedXPowGate(phase_exponent=t), vqubit.PhasedXPowGate(phase_exponent=0.5)),
		("phased_xz", vqubit.PhasedXZGate(x_exponent=t), vqubit.PhasedXZGate(x_exponent=0.5)),
		("ms", vqubit.MSGate(rads=t), vqubit.MSGate(rads=0.5)),
		("fsim", vqubit.FSimGate(theta=t, phi=0.1), vqubit.FSimGate(theta=0.5, phi=0.1)),
		("global_phase", vqubit.GlobalPhaseGate(t), vqubit.GlobalPhaseGate(0.5)),
		("rz", vqubit.rz(t / math.pi), vqubit.rz(0.5 / math.pi)),
	]


@pytest.mark.parametrize(("name", "symbolic", "resolved"), symbolic_gates(), ids=lambda value: str(value)[:12])
def test_a_symbolic_gate_has_no_unitary(name: str, symbolic: Any, resolved: Any) -> None:
	assert vqubit.is_parameterized(symbolic), name
	assert vqubit.unitary(resolved) is not None, name
	# Some gates fail when the fallback decomposition needs a concrete
	# exponent (the Toffoli, for one), so both errors are acceptable.
	with pytest.raises((TypeError, ValueError)):
		vqubit.unitary(symbolic)


@pytest.mark.parametrize(("name", "symbolic", "resolved"), symbolic_gates(), ids=lambda value: str(value)[:12])
def test_resolving_a_symbolic_gate_gives_the_concrete_gate(name: str, symbolic: Any, resolved: Any) -> None:
	resolver = vqubit.ParamResolver({"t": 0.5})
	resolved_gate = vqubit.resolve_parameters(symbolic, resolver)
	assert not vqubit.is_parameterized(resolved_gate), name
	assert resolved_gate == resolved, name
	assert_matrix_close(vqubit.unitary(resolved_gate), vqubit.unitary(resolved))


def test_a_concrete_gate_resolves_to_itself() -> None:
	gate = vqubit.XPowGate(exponent=0.25)
	assert vqubit.resolve_parameters(gate, vqubit.ParamResolver({"t": 1})) is gate


def test_a_partially_resolved_gate_stays_symbolic() -> None:
	s = vqubit.Parameter("s")
	gate = vqubit.PhasedXZGate(x_exponent=0.5, z_exponent=s)
	resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({"t": 1}))
	assert vqubit.is_parameterized(resolved)
	assert resolved.x_exponent == 0.5


def test_powers_of_a_symbolic_gate_keep_the_symbol() -> None:
	t = vqubit.Parameter("t")
	gate = vqubit.ZPowGate(exponent=t) ** 2
	assert vqubit.is_parameterized(gate)
	assert vqubit.resolve_parameters(gate, vqubit.ParamResolver({"t": 0.5})) == vqubit.ZPowGate(exponent=1.0)
	# ... and the inverse of a symbolic gate is still symbolic.
	assert vqubit.is_parameterized(vqubit.ZPowGate(exponent=t)._inverse_())


def test_unresolved_parameters_keep_symbolic_equality() -> None:
	# Two gates written with the same parameters are equal; a different
	# expression is not.
	t = vqubit.Parameter("t")
	assert vqubit.XPowGate(exponent=t) == vqubit.XPowGate(exponent=t)
	assert vqubit.XPowGate(exponent=t) != vqubit.XPowGate(exponent=t * 2)
	# The symbolic gate only compares equal to itself, not to the concrete one.
	assert vqubit.XPowGate(exponent=t) != vqubit.XPowGate(exponent=0.5)


@pytest.mark.parametrize(("name", "symbolic", "resolved"), symbolic_gates(), ids=lambda value: str(value)[:12])
def test_reprs_rebuild_the_gate(name: str, symbolic: Any, resolved: Any) -> None:
	namespace = {"vqubit": vqubit}
	for gate in (symbolic, resolved):
		assert eval(repr(gate), namespace) == gate, name


def test_diagram_symbols_of_the_odd_gates() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	# A phased X shows up as PhX, a matrix gate as Matrix, a wait as Wait and
	# the Mølmer-Sørensen gate as MS.
	assert "PhX" in str(vqubit.Circuit(vqubit.PhasedXPowGate(phase_exponent=0.25).on(q0)))
	matrix_gate = vqubit.MatrixGate([[0, 1], [1, 0]])
	assert "Matrix" in str(vqubit.Circuit(matrix_gate.on(q0)))
	assert "Wait" in str(vqubit.Circuit(vqubit.WaitGate(10.0).on(q0)))
	assert "MS" in str(vqubit.Circuit(vqubit.MSGate(rads=0.5).on(q0, q1)))


def test_str_of_the_plain_gates() -> None:
	assert str(vqubit.PhasedXPowGate(phase_exponent=0.25)) == "PhasedXPowGate(phase_exponent=0.25, exponent=1)"
	assert str(vqubit.WaitGate(5.0, 2)) == "WaitGate(5.0, num_qubits=2)"
	assert str(vqubit.MatrixGate([[1, 0], [0, 1]])) == "MatrixGate(1q)"
	assert str(vqubit.MSGate(rads=0.5)) == "MSGate(rads=0.5)"
	assert str(vqubit.CSWAP) == "CSWAPGate"


def test_matrix_gate_validates_its_shape() -> None:
	with pytest.raises(ValueError, match="power of two"):
		vqubit.MatrixGate([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
	with pytest.raises(ValueError, match="square"):
		vqubit.MatrixGate([[1, 0], [0, 1], [0, 0], [0, 0]])


def test_matrix_gate_is_its_matrix() -> None:
	matrix: list[list[complex]] = [[0, 1j], [1, 0]]
	gate = vqubit.MatrixGate(matrix)
	assert gate.num_qubits == 1
	assert_matrix_close(gate.matrix, matrix)
	assert_matrix_close(vqubit.unitary(gate.on(vqubit.LineQubit(0))), [[0, 1j], [1, 0]])
	# The four by four matrix is a two-qubit gate.
	assert vqubit.MatrixGate([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]).num_qubits == 2


def test_controlled_gate_validates_its_arguments() -> None:
	with pytest.raises(ValueError, match="at least one control"):
		vqubit.ControlledGate(vqubit.X, num_controls=0)
	with pytest.raises(TypeError, match="needs a Gate"):
		vqubit.ControlledGate("X", num_controls=1)  # ty: ignore[invalid-argument-type]
	with pytest.raises(ValueError, match="control_values"):
		vqubit.ControlledGate(vqubit.X, num_controls=2, control_values=[True])


def test_controlled_gate_acts_only_on_its_control_value() -> None:
	# The X only fires on the |11> branch, which is the CNOT.
	gate = vqubit.ControlledGate(vqubit.X, num_controls=1)
	assert_matrix_close(vqubit.unitary(gate), vqubit.unitary(vqubit.CNOT))
	# ... and the open circle means the control is active on |0>.
	inverted = vqubit.ControlledGate(vqubit.X, num_controls=1, control_values=[False])
	q0, q1 = vqubit.LineQubit.range(2)
	state = vqubit.Simulator().simulate(vqubit.Circuit(vqubit.X(q1), inverted.on(q0, q1))).state_vector()
	# |01> -- the inverted control fires and flips the target back to |00>.
	assert abs(state[0]) == pytest.approx(1.0)
	assert "○" in str(vqubit.Circuit(inverted.on(q0, q1)))
	assert "@" in str(vqubit.Circuit(gate.on(q0, q1)))


def test_controlled_gate_of_a_parameterized_gate_is_parameterized() -> None:
	t = vqubit.Parameter("t")
	gate = vqubit.ControlledGate(vqubit.ZPowGate(exponent=t), num_controls=1)
	assert vqubit.is_parameterized(gate)
	resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({"t": 1}))
	assert resolved == vqubit.ControlledGate(vqubit.Z, num_controls=1)


def test_wait_gate_validates_its_duration() -> None:
	with pytest.raises(ValueError, match="non-negative"):
		vqubit.WaitGate(-1.0)
	assert vqubit.WaitGate(0.0).duration_nanos == 0.0


def test_wait_gate_is_the_identity() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	assert_matrix_close(
		vqubit.unitary(vqubit.Circuit(vqubit.WaitGate(5.0, 2).on(q0, q1))),
		[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
	)


def test_the_toffoli_decomposition_only_exists_at_exponent_one() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	with pytest.raises(ValueError, match="only applies at exponent 1"):
		vqubit.decompose_once((vqubit.CCX**0.5).on(q0, q1, q2))
	expanded = vqubit.decompose_once(vqubit.CCX.on(q0, q1, q2))
	assert_matrix_close(
		vqubit.unitary(vqubit.Circuit(*expanded)),
		vqubit.unitary(vqubit.Circuit(vqubit.CCX.on(q0, q1, q2))),
	)


def test_measurement_gate_requires_a_string_key() -> None:
	q0 = vqubit.LineQubit(0)
	with pytest.raises(TypeError, match="key must be a str"):
		vqubit.MeasurementGate(1, 3)  # ty: ignore[invalid-argument-type]
	assert repr(vqubit.measure(q0, key="m").gate) == "vqubit.MeasurementGate('m', num_qubits=1, invert_mask=())"


def test_global_phase_gate_acts_on_no_qubits() -> None:
	gate = vqubit.global_phase(0.5)
	assert gate.num_qubits == 0
	assert_matrix_close(vqubit.unitary(gate), [[1j]])
	# The circuit is a single moment with no qubits, and the diagram is empty.
	circuit = vqubit.Circuit(gate.on())
	assert circuit.all_qubits() == ()
	assert str(circuit).strip() == ""
