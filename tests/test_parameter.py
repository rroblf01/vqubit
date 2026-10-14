"""Tests for symbolic parameters and their resolution."""

from __future__ import annotations

import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close
from vqubit._parameter import resolve_value


def test_parameter_basics() -> None:
	t = vqubit.Parameter("t")
	assert t.name == "t"
	assert str(t) == "t"
	assert repr(t) == "vqubit.Parameter('t')"
	assert vqubit.Parameter("t") == vqubit.Parameter("t")
	assert vqubit.Parameter("t") != vqubit.Parameter("u")
	assert hash(vqubit.Parameter("t")) == hash(vqubit.Parameter("t"))


def test_parameter_validates_name() -> None:
	with pytest.raises(ValueError):
		vqubit.Parameter("")


def test_parameter_arithmetic_and_names() -> None:
	t = vqubit.Parameter("t")
	expression = (t + 1) * 2 - t / 3
	assert vqubit.is_parameterized(expression)
	assert expression.names() == {"t"}


def test_param_resolver() -> None:
	t = vqubit.Parameter("t")
	resolver = vqubit.ParamResolver({t: 3})
	assert resolver[t] == 3
	assert resolver["t"] == 3
	assert resolver.value_of(t) == 3
	assert resolver.value_of(5) == 5
	assert len(resolver) == 1
	assert list(resolver) == ["t"]


def test_param_resolver_supplies_kwargs() -> None:
	assert vqubit.ParamResolver(theta=0.5).value_of(vqubit.Parameter("theta")) == 0.5


def test_expression_substitution() -> None:
	t = vqubit.Parameter("t")
	resolver = vqubit.ParamResolver({t: 4})
	assert resolve_value((t + 1) * 2, resolver) == 10
	assert resolve_value(t**2, resolver) == 16
	assert resolve_value(-t, resolver) == -4


def test_unresolved_parameter_stays_symbolic() -> None:
	t = vqubit.Parameter("t")
	assert resolve_value(t + 1, vqubit.ParamResolver()) == t + 1


def test_gate_is_parameterized() -> None:
	t = vqubit.Parameter("t")
	assert vqubit.is_parameterized(vqubit.Rx(t))
	assert not vqubit.is_parameterized(vqubit.Rx(0.5))
	assert vqubit.is_parameterized(vqubit.XPowGate(exponent=t))
	assert not vqubit.is_parameterized(vqubit.X)


def test_parameterized_gate_unitary_raises() -> None:
	t = vqubit.Parameter("t")
	with pytest.raises(TypeError):
		vqubit.unitary(vqubit.Rx(t))


def test_resolve_parameters_of_gate() -> None:
	t = vqubit.Parameter("t")
	gate = vqubit.Rx(t / 2)
	resolved = vqubit.resolve_parameters(gate, vqubit.ParamResolver({t: math.pi / 2}))
	assert resolved == vqubit.Rx(math.pi / 4)
	assert_matrix_close(vqubit.unitary(resolved), vqubit.unitary(vqubit.Rx(math.pi / 4)))


def test_resolve_parameters_of_non_parameterized_gate_is_identity() -> None:
	gate = vqubit.H
	assert vqubit.resolve_parameters(gate, vqubit.ParamResolver()) is gate


def test_operation_and_circuit_parameterization() -> None:
	t = vqubit.Parameter("t")
	q0 = vqubit.LineQubit(0)
	operation = vqubit.Rz(t)(q0)
	assert vqubit.is_parameterized(operation)
	resolved = vqubit.resolve_parameters(operation, vqubit.ParamResolver({t: 0.5}))
	assert resolved == vqubit.Rz(0.5)(q0)

	circuit = vqubit.Circuit(vqubit.Rz(t)(q0), vqubit.H(q0))
	assert vqubit.is_parameterized(circuit)
	resolved_circuit = vqubit.resolve_parameters(circuit, vqubit.ParamResolver({t: 0.5}))
	assert not vqubit.is_parameterized(resolved_circuit)
	assert list(resolved_circuit.all_operations()) == [vqubit.Rz(0.5)(q0), vqubit.H(q0)]
