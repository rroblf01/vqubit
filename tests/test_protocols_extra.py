"""Tests for the extra protocol helpers in ``vqubit._protocols_extra``."""

from __future__ import annotations

import dataclasses
import math

import pytest

import vqubit
from tests.helpers import assert_matrix_close, assert_vector_close
from vqubit import _protocols_extra as pe


class _Action:
	def __init__(self) -> None:
		self.seen: tuple[object, object] | None = None

	def _act_on_(self, args: object, qubits: object) -> str:
		self.seen = (args, qubits)
		return "acted"


class _ApplyOperationArgs:
	def __init__(self) -> None:
		self.operation: object = None

	def apply_operation(self, operation: object) -> None:
		self.operation = operation


class _Channel:
	def _apply_channel_(self, args: object) -> tuple[str, object]:
		return "channel", args


class _Mixture:
	def _apply_mixture_(self, args: object) -> tuple[str, object]:
		return "mixture", args


class _Symbols:
	num_qubits = 2

	def _diagram_symbols_(self, qubits: object) -> tuple[str, str]:
		return "A", "B"


class _DiagramInfo:
	num_qubits = 1

	def _circuit_diagram_info_(self, args: object) -> pe.CircuitDiagramInfo:
		return pe.CircuitDiagramInfo(("Q",), exponent=0.5)


class _TupleDiagram:
	num_qubits = 1

	def _circuit_diagram_info_(self, args: object) -> tuple[str, ...]:
		return ("Z",)


class _Plain:
	num_qubits = 2


class _Decomposer:
	def _decompose_(self, qubits: list[vqubit.Qid]) -> list[vqubit.Operation]:
		return [vqubit.X(qubits[0]), vqubit.Y(qubits[0])]


class _EmptyDecomposer:
	def _decompose_(self, qubits: list[vqubit.Qid]) -> None:
		return None


@dataclasses.dataclass
class _Point:
	x: int = 1
	y: int = 2


def test_apply_unitary_args_subspace_index() -> None:
	args = pe.ApplyUnitaryArgs([0] * 8, axes=[1, 2])
	assert args.subspace_index() == (0, 0)
	assert args.subspace_index(1) == (0, 1)
	assert args.subspace_index(3) == (1, 1)
	assert pe.ApplyUnitaryArgs([0] * 2).subspace_index() == ()


def test_apply_unitary_args_defaults_and_axes_transposed() -> None:
	tensor = [1, 0]
	buffer = [0, 0]
	args = pe.ApplyUnitaryArgs(tensor, buffer, [0])
	assert args.target_tensor is tensor
	assert args.available_buffer is buffer
	assert args.axes == [0]
	# Without an explicit buffer the target is reused.
	assert pe.ApplyUnitaryArgs(tensor).available_buffer is tensor
	assert args.with_axes_transposed_to_start() is args


def test_apply_channel_and_mixture_args_defaults() -> None:
	tensor = [1, 0]
	channel_args = pe.ApplyChannelArgs(tensor, left_axes=[0])
	assert channel_args.target_tensor is tensor
	assert channel_args.out_buffer is tensor
	assert channel_args.auxiliary_buffer0 is None
	assert channel_args.auxiliary_buffer1 is None
	assert channel_args.left_axes == [0]
	assert channel_args.right_axes == []
	mixture_args = pe.ApplyMixtureArgs(tensor, right_axes=[1])
	assert mixture_args.out_buffer is tensor
	assert mixture_args.auxiliary_buffer is None
	assert mixture_args.left_axes == []
	assert mixture_args.right_axes == [1]


def test_act_on_prefers_the_act_on_protocol() -> None:
	action = _Action()
	q0 = vqubit.LineQubit(0)
	assert pe.act_on(action, "state", (q0,)) == "acted"
	assert action.seen == ("state", (q0,))


def test_act_on_falls_back_to_apply_operation() -> None:
	args = _ApplyOperationArgs()
	assert pe.act_on(vqubit.X(q0 := vqubit.LineQubit(0)), args) is True
	assert args.operation == vqubit.X(q0)


def test_act_on_returns_not_implemented_when_unsupported() -> None:
	assert pe.act_on(object(), object()) is NotImplemented


def test_apply_unitaries_composes_in_order() -> None:
	args = pe.ApplyUnitaryArgs([1, 0, 0, 0], axes=[0])
	assert_vector_close(pe.apply_unitaries([vqubit.X], args), [0, 0, 1, 0])
	# X then Z: Z X |0> = -|1>
	assert_vector_close(pe.apply_unitaries([vqubit.X, vqubit.Z], args), [0, 0, -1, 0])


def test_apply_unitaries_without_axes_returns_flat_tensor() -> None:
	args = pe.ApplyUnitaryArgs([1, 0, 0, 0])
	assert_vector_close(pe.apply_unitaries([vqubit.X], args), [1, 0, 0, 0])
	# Nested tensors are flattened when no axes are given.
	nested = pe.ApplyUnitaryArgs([[1, 0], [0, 1]])
	assert_vector_close(pe.apply_unitaries([vqubit.X], nested), [1, 0, 0, 1])


def test_apply_channel_and_mixture_dispatch() -> None:
	channel_args = pe.ApplyChannelArgs([1, 0])
	assert pe.apply_channel(_Channel(), channel_args) == ("channel", channel_args)
	mixture_args = pe.ApplyMixtureArgs([1, 0])
	assert pe.apply_mixture(_Mixture(), mixture_args) == ("mixture", mixture_args)


def test_apply_channel_and_mixture_defaults() -> None:
	channel_args = pe.ApplyChannelArgs([1, 0])
	mixture_args = pe.ApplyMixtureArgs([1, 0])
	assert pe.apply_channel(vqubit.X, channel_args) is None
	assert pe.apply_mixture(vqubit.X, mixture_args) is None
	assert pe.apply_channel(vqubit.X, channel_args, default="none") == "none"
	assert pe.apply_mixture(vqubit.X, mixture_args, default="none") == "none"


def test_mixture_and_has_mixture_for_channels_and_operations() -> None:
	identity = [[1, 0], [0, 1]]
	x = [[0, 1], [1, 0]]
	channel = vqubit.MixedUnitaryChannel([(0.25, identity), (0.75, x)])
	assert pe.has_mixture(channel)
	entries = pe.mixture(channel)
	assert entries is not None
	assert [probability for probability, _ in entries] == pytest.approx([0.25, 0.75])
	assert_matrix_close(entries[1][1], x)
	operation = channel.on(vqubit.LineQubit(0))
	assert pe.has_mixture(operation)
	assert pe.mixture(operation) is not None


def test_mixture_absent_returns_none() -> None:
	assert not pe.has_mixture(vqubit.X)
	assert pe.mixture(vqubit.X) is None


def test_validate_mixture() -> None:
	pe.validate_mixture([(0.25, None), (0.75, None)])
	with pytest.raises(ValueError, match="sum to 1"):
		pe.validate_mixture([(0.2, None), (0.2, None)])
	with pytest.raises(ValueError, match="sum to 1"):
		pe.validate_mixture([])


def test_has_stabilizer_effect() -> None:
	assert pe.has_stabilizer_effect(vqubit.X)
	assert pe.has_stabilizer_effect(vqubit.CNOT)
	assert not pe.has_stabilizer_effect(vqubit.T)
	assert pe.has_stabilizer_effect(vqubit.measure(vqubit.LineQubit(0), key="m"))

	class Explicit:
		def _has_stabilizer_effect_(self) -> bool:
			return False

	assert not pe.has_stabilizer_effect(Explicit())


def test_circuit_diagram_info_for_gates_and_operations() -> None:
	q0 = vqubit.LineQubit(0)
	assert pe.circuit_diagram_info(vqubit.X) == pe.CircuitDiagramInfo(("X",))
	assert pe.circuit_diagram_info(vqubit.CNOT).wire_symbols == ("@", "X")
	assert pe.circuit_diagram_info(vqubit.measure(q0, key="m")).wire_symbols == ("M",)
	assert pe.circuit_diagram_info(vqubit.X(q0)).wire_symbols == ("X",)


def test_circuit_diagram_info_custom_protocols_and_default() -> None:
	assert pe.circuit_diagram_info(_Symbols()).wire_symbols == ("A", "B")
	info = pe.circuit_diagram_info(_DiagramInfo())
	assert info.wire_symbols == ("Q",)
	assert info.exponent == 0.5
	assert pe.circuit_diagram_info(_TupleDiagram()).wire_symbols == ("Z",)
	assert len(pe.circuit_diagram_info(_Plain()).wire_symbols) == 2


def test_circuit_diagram_info_dataclasses() -> None:
	info = pe.CircuitDiagramInfo(("X",), exponent=0.5, connected=False)
	renamed = info.with_wire_symbols(["Y", "Z"])
	assert renamed.wire_symbols == ("Y", "Z")
	assert renamed.exponent == 0.5
	assert renamed.connected is False
	args = pe.CircuitDiagramInfoArgs(precision=3)
	assert args.format_real(1.5) == "1.5"
	context = pe.DecompositionContext()
	assert context.qubits == ()
	assert context.atol == 1e-8
	assert pe.DEFAULT_RESOLVERS == ()


def test_decompose_once_with_qubits() -> None:
	q0, q1, q2 = vqubit.LineQubit.range(3)
	operations = pe.decompose_once_with_qubits(vqubit.CCX, [q0, q1, q2])
	assert operations[0] == vqubit.H(q2)
	assert vqubit.equal_up_to_global_phase(vqubit.Circuit(operations), vqubit.CCX(q0, q1, q2))
	assert len(pe.decompose_once_with_qubits(_Decomposer(), [q0])) == 2
	assert pe.decompose_once_with_qubits(_EmptyDecomposer(), [q0]) == []
	assert pe.decompose_once_with_qubits(vqubit.X, [q0]) == []
	assert pe.decompose_once_with_qubits(object(), [q0]) == []


def test_parameter_names_in_nested_equality_values() -> None:
	parameter = vqubit.Parameter("t")

	class Nested:
		def _value_equality_values_(self) -> tuple[object, ...]:
			return ((parameter, 1),)

	assert pe.parameter_names(Nested()) == {"t"}
	assert pe.parameter_symbols(Nested()) == {vqubit.Parameter("t")}


def test_parameter_names_and_symbols() -> None:
	parameter = vqubit.Parameter("t")
	gate = vqubit.XPowGate(exponent=parameter)
	assert pe.parameter_names(gate) == {"t"}
	assert pe.parameter_symbols(gate) == {vqubit.Parameter("t")}
	assert pe.parameter_names(vqubit.X) == set()


def test_resolve_parameters_once() -> None:
	parameter = vqubit.Parameter("t")
	gate = vqubit.XPowGate(exponent=parameter)
	resolved = pe.resolve_parameters_once(gate, vqubit.ParamResolver({"t": 0.5}))
	assert resolved == vqubit.X**0.5
	assert not vqubit.is_parameterized(resolved)
	assert pe.resolve_parameters_once(vqubit.X, vqubit.ParamResolver()) is vqubit.X
	plain = object()
	assert pe.resolve_parameters_once(plain, vqubit.ParamResolver()) is plain


def test_mul_handles_symbolic_and_unsupported_values() -> None:
	assert pe.mul(2, 3) == 6
	parameter = vqubit.Parameter("t")
	assert "t * 2" in str(pe.mul(parameter, 2))
	assert "2 * t" in str(pe.mul(2, parameter))
	assert pe.mul(object(), 2, default="default") == "default"


def test_measurement_key_name() -> None:
	q0 = vqubit.LineQubit(0)
	operation = vqubit.measure(q0, key="alpha")
	assert pe.measurement_key_name(operation) == "alpha"
	assert pe.measurement_key_name(operation.gate) == "alpha"
	assert pe.measurement_key_name(vqubit.X) is None
	assert pe.measurement_key_name(q0) is None


def test_with_key_path_helpers() -> None:
	key = pe.with_key_path("k", ["a", "b"])
	assert isinstance(key, vqubit.MeasurementKey)
	assert str(key) == "a.b.k"
	assert key.name == "a.b.k"
	assert key.path == ()
	assert str(pe.with_key_path_prefix(vqubit.MeasurementKey("k"), ["a"])) == "a.k"
	assert str(pe.with_key_path("k", [])) == "k"
	assert str(pe.with_rescoped_keys("k", ["a", "b"])) == "a.b.k"


def test_json_namespace() -> None:
	assert pe.json_namespace(vqubit.X) == ""
	assert pe.json_namespace(pe.HasJSONNamespace()) == ""

	class Namespaced:
		def _json_namespace_(self) -> str:
			return "custom"

	assert pe.json_namespace(Namespaced()) == "custom"


def test_dataclass_json_dict() -> None:
	assert pe.dataclass_json_dict(_Point(3, 4)) == {"x": 3, "y": 4}
	assert pe.dataclass_json_dict(_Point()) == {"x": 1, "y": 2}
	assert pe.dataclass_json_dict(_Point) == {}
	assert pe.dataclass_json_dict(vqubit.X) == {}


def test_obj_to_dict_helper() -> None:
	entity = pe.LabelEntity(name="n", display_name="d")
	assert pe.obj_to_dict_helper(entity, ["name", "display_name"]) == {"name": "n", "display_name": "d"}
	with_namespace = pe.obj_to_dict_helper(entity, ["name"], namespace="vqubit.LabelEntity")
	assert with_namespace == {"name": "n", "cirq_type": "vqubit.LabelEntity"}


def test_json_resolver() -> None:
	resolver = pe.JsonResolver({"Thing": int})
	assert "Thing" in resolver
	assert "Other" not in resolver
	assert resolver["Thing"] is int
	with pytest.raises(KeyError):
		resolver["Missing"]
	assert "Thing" not in pe.JsonResolver()


def test_serializable_by_key() -> None:
	with pytest.raises(NotImplementedError):
		pe.SerializableByKey()._serialization_key_()

	class Keyed(pe.SerializableByKey):
		def _serialization_key_(self) -> str:
			return "keyed"

	assert Keyed()._serialization_key_() == "keyed"


def test_label_entity_defaults() -> None:
	entity = pe.LabelEntity(display_name="Display")
	assert entity.name == "Display"
	assert entity.display_name == "Display"
	assert entity.description == ""
	assert str(entity) == "Display"
	assert "LabelEntity" in repr(entity)
	assert pe.LabelEntity(name="n").display_name == "n"


def test_trace_distance_from_angle_list() -> None:
	assert pe.trace_distance_from_angle_list([]) == 0.0
	assert pe.trace_distance_from_angle_list([0, 0]) == pytest.approx(0)
	assert pe.trace_distance_from_angle_list([math.pi, math.pi]) == pytest.approx(0)
	# A single angle always has |mean phase| = 1, so the bound is zero.
	assert pe.trace_distance_from_angle_list([math.pi / 2]) == pytest.approx(0)
	assert pe.trace_distance_from_angle_list([0, math.pi]) == pytest.approx(1)
	assert pe.trace_distance_from_angle_list([0, math.pi / 2]) == pytest.approx(math.sqrt(0.5))


def test_supports_protocols_runtime_checks() -> None:
	class Unit:
		def _unitary_(self) -> list[list[complex]]:
			return [[1]]

	class Kraus:
		def _kraus_(self) -> list[list[list[complex]]]:
			return []

	class Apply:
		def _apply_unitary_(self, args: object) -> object:
			return args

	class Act:
		def _act_on_(self, args: object, qubits: object = None) -> bool:
			return True

	class Expansion:
		def _pauli_expansion_(self) -> object:
			return None

	assert isinstance(Unit(), pe.SupportsUnitary)
	assert isinstance(Kraus(), pe.SupportsKraus)
	assert isinstance(Apply(), pe.SupportsApplyUnitary)
	assert isinstance(Act(), pe.SupportsActOn)
	assert isinstance(Expansion(), pe.SupportsPauliExpansion)
	assert not isinstance(Unit(), pe.SupportsKraus)
	assert not isinstance(object(), pe.SupportsUnitary)
	assert isinstance(vqubit.X, pe.SupportsUnitary)
	assert isinstance(vqubit.X, pe.SupportsParameterization)
	assert not isinstance(vqubit.X, pe.SupportsMixture)
	assert pe.SupportsConsistentApplyUnitary is pe.SupportsApplyUnitary
	assert pe.SupportsActOnQubits is pe.SupportsActOn
