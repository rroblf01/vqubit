"""Extra ``cirq.protocols`` helpers: apply args, diagram info, mixtures,
JSON helpers, key paths and the ``Supports*`` protocol classes.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Iterable, Mapping, Sequence
from typing import Any, Protocol, runtime_checkable

from ._linalg import Matrix
from ._measurement_key import MeasurementKey
from ._parameter import Parameter, ParamResolver, is_symbolic
from ._qubits import LineQubit, Qid

# --- apply arguments ---------------------------------------------------------


class ApplyUnitaryArgs:
	"""The arguments passed to a gate's ``_apply_unitary_`` method."""

	def __init__(self, target_tensor: Any, available_buffer: Any = None, axes: Sequence[int] | None = None) -> None:
		self.target_tensor = target_tensor
		self.available_buffer = target_tensor if available_buffer is None else available_buffer
		self.axes = list(axes) if axes is not None else []

	def subspace_index(self, big_endian_bits_int: int = 0) -> Any:
		"""Return the index into ``target_tensor`` for a subspace basis state."""
		index: list[Any] = []
		for axis in range(len(self.axes)):
			index.append((big_endian_bits_int >> (len(self.axes) - 1 - axis)) & 1)
		return tuple(index)

	def with_axes_transposed_to_start(self) -> ApplyUnitaryArgs:
		return self


class ApplyChannelArgs:
	"""The arguments passed to a channel's ``_apply_channel_`` method."""

	def __init__(
		self,
		target_tensor: Any,
		out_buffer: Any = None,
		auxiliary_buffer0: Any = None,
		auxiliary_buffer1: Any = None,
		left_axes: Sequence[int] | None = None,
		right_axes: Sequence[int] | None = None,
	) -> None:
		self.target_tensor = target_tensor
		self.out_buffer = target_tensor if out_buffer is None else out_buffer
		self.auxiliary_buffer0 = auxiliary_buffer0
		self.auxiliary_buffer1 = auxiliary_buffer1
		self.left_axes = list(left_axes) if left_axes is not None else []
		self.right_axes = list(right_axes) if right_axes is not None else []


class ApplyMixtureArgs:
	"""The arguments passed to a channel's ``_apply_mixture_`` method."""

	def __init__(
		self,
		target_tensor: Any,
		out_buffer: Any = None,
		auxiliary_buffer: Any = None,
		left_axes: Sequence[int] | None = None,
		right_axes: Sequence[int] | None = None,
	) -> None:
		self.target_tensor = target_tensor
		self.out_buffer = target_tensor if out_buffer is None else out_buffer
		self.auxiliary_buffer = auxiliary_buffer
		self.left_axes = list(left_axes) if left_axes is not None else []
		self.right_axes = list(right_axes) if right_axes is not None else []


def act_on(action: Any, args: Any, qubits: Sequence[Qid] | None = None, *, allow_decompose: bool = True) -> Any:
	"""Apply ``action`` to a simulation state ``args``."""
	del allow_decompose
	if hasattr(action, "_act_on_"):
		return action._act_on_(args, qubits)
	if hasattr(args, "apply_operation"):
		args.apply_operation(action)
		return True
	return NotImplemented


def apply_unitaries(unitaries: Iterable[Any], args: ApplyUnitaryArgs, *, default: Any = None) -> Any:
	"""Apply a sequence of unitaries to the target tensor of ``args``."""
	result = args.target_tensor
	for value in unitaries:
		from ._protocols import unitary

		matrix = unitary(value)
		result = _left_multiply(matrix, result, args.axes)
	if result is None:
		return default
	return result


def _left_multiply(matrix: Matrix, tensor: Any, axes: Sequence[int]) -> Any:
	flat = _flatten(tensor)
	from ._linalg import apply_gate

	return apply_gate(flat, matrix, list(axes)) if axes else flat


def _flatten(tensor: Any) -> list[complex]:
	if isinstance(tensor, Sequence) and tensor and not isinstance(tensor[0], Sequence):
		return [complex(z) for z in tensor]
	if isinstance(tensor, Sequence):
		return [complex(z) for row in tensor for z in row]
	return [complex(tensor)]


def apply_channel(channel: Any, args: ApplyChannelArgs, *, default: Any = None) -> Any:
	"""Apply a channel to the target tensor of ``args``."""
	if hasattr(channel, "_apply_channel_"):
		return channel._apply_channel_(args)
	return default


def apply_mixture(value: Any, args: ApplyMixtureArgs, *, default: Any = None) -> Any:
	"""Apply a mixture to the target tensor of ``args``."""
	if hasattr(value, "_apply_mixture_"):
		return value._apply_mixture_(args)
	return default


# --- mixtures ----------------------------------------------------------------


def has_mixture(value: Any) -> bool:
	"""Report whether ``value`` has a mixture representation."""
	return hasattr(value, "_mixture_") or hasattr(getattr(value, "gate", None), "_mixture_")


def mixture(value: Any) -> list[tuple[float, Any]] | None:
	"""Return the mixture of ``value``, or ``None``."""
	subject = getattr(value, "gate", value)
	if hasattr(subject, "_mixture_"):
		return list(subject._mixture_())
	return None


def validate_mixture(mixture: Sequence[tuple[float, Any]]) -> None:
	"""Raise ``ValueError`` if ``mixture`` is not a valid mixture."""
	total = sum(probability for probability, _ in mixture)
	if abs(total - 1) > 1e-6:
		raise ValueError(f"mixture probabilities must sum to 1, got {total}")


def has_stabilizer_effect(value: Any) -> bool:
	"""Report whether ``value`` acts as a Clifford (stabilizer) operation."""
	subject = getattr(value, "gate", value)
	if hasattr(subject, "_has_stabilizer_effect_"):
		return bool(subject._has_stabilizer_effect_())
	from ._clifford import _is_clifford_gate
	from ._gates import MeasurementGate

	return isinstance(subject, MeasurementGate) or _is_clifford_gate(subject)


# --- circuit diagram ---------------------------------------------------------


@dataclasses.dataclass
class CircuitDiagramInfo:
	"""The wire symbols and exponent used to draw a gate."""

	wire_symbols: tuple[str, ...]
	exponent: Any = 1
	connected: bool = True

	def with_wire_symbols(self, wire_symbols: Sequence[str]) -> CircuitDiagramInfo:
		return CircuitDiagramInfo(tuple(wire_symbols), self.exponent, self.connected)


@dataclasses.dataclass
class CircuitDiagramInfoArgs:
	"""The context passed to a gate's ``_circuit_diagram_info_`` method."""

	known_qubits: tuple[Qid, ...] | None = None
	known_qubit_count: int | None = None
	qubit_map: Mapping[Qid, int] | None = None
	precision: int | None = None
	qubit_map_entries: Any = None
	use_unicode_characters: bool = True

	def format_real(self, value: Any) -> str:
		return str(value)


def circuit_diagram_info(value: Any, args: CircuitDiagramInfoArgs | None = None) -> CircuitDiagramInfo:
	"""Return the :class:`CircuitDiagramInfo` for ``value``."""
	subject = getattr(value, "gate", value)
	qubits = getattr(value, "qubits", None)
	if qubits is None:
		qubits = tuple(LineQubit.range(getattr(subject, "num_qubits", 1)))
	if hasattr(subject, "_diagram_symbols_"):
		return CircuitDiagramInfo(tuple(subject._diagram_symbols_(qubits)))
	if hasattr(subject, "_circuit_diagram_info_"):
		result = subject._circuit_diagram_info_(args)
		return result if isinstance(result, CircuitDiagramInfo) else CircuitDiagramInfo(tuple(result))
	return CircuitDiagramInfo((str(subject),) * len(qubits))


# --- decomposition -----------------------------------------------------------


@dataclasses.dataclass
class DecompositionContext:
	"""The context passed to ``_decompose_`` methods."""

	qubits: tuple[Qid, ...] = ()
	preserve_measurements: bool = False
	preserve_tags: bool = False
	atol: float = 1e-8


def decompose_once_with_qubits(
	value: Any, qubits: Sequence[Qid], *, context: DecompositionContext | None = None
) -> list[Any]:
	"""Decompose ``value`` one layer, on the given ``qubits``."""
	del context
	subject = getattr(value, "gate", value)
	if hasattr(subject, "_decompose_"):
		result = subject._decompose_(tuple(qubits))
		return list(result) if result is not None else []
	return []


# --- parameters --------------------------------------------------------------


DEFAULT_RESOLVERS: tuple[ParamResolver, ...] = ()


def _parameters_in(value: Any, depth: int = 0) -> Iterable[Parameter]:
	"""The parameters named by a value, looking inside expressions and children."""
	if depth > 8:
		return
	if isinstance(value, Parameter):
		yield value
		return
	names = getattr(value, "names", None)
	if callable(names) and not isinstance(value, (list, tuple, set, dict)):
		for name in sorted(names()):
			yield Parameter(name)
		return
	if isinstance(value, (tuple, list, set, frozenset)):
		for item in value:
			yield from _parameters_in(item, depth + 1)
		return
	inner = getattr(value, "_value_equality_values_", None)
	if callable(inner):
		try:
			values = inner()
		except Exception:
			values = None
		if values is not None:
			yield from _parameters_in(values, depth + 1)
	# Gates keep their exponent, angle or parameter in a private attribute, so
	# walking the object's state is what finds the symbol inside a circuit.
	attributes = getattr(value, "__dict__", None)
	if isinstance(attributes, dict):
		for item in attributes.values():
			yield from _parameters_in(item, depth + 1)


def _iter_parameters(value: Any) -> Iterable[Parameter]:
	return _parameters_in(value)


def parameter_names(value: Any) -> set[str]:
	"""Return the set of parameter names used by ``value``."""
	return {parameter.name for parameter in _iter_parameters(value)}


def parameter_symbols(value: Any) -> set[Parameter]:
	"""Return the set of parameters used by ``value``."""
	return set(_iter_parameters(value))


def resolve_parameters_once(value: Any, resolver: ParamResolver) -> Any:
	"""Resolve one layer of parameters in ``value``."""
	if hasattr(value, "_resolve_parameters_"):
		return value._resolve_parameters_(resolver)
	return value


def mul(a: Any, b: Any, default: Any = NotImplemented) -> Any:
	"""Multiply two values, returning ``default`` when unsupported."""
	if is_symbolic(a) or is_symbolic(b):
		try:
			return a * b
		except TypeError:
			return default
	try:
		return a * b
	except TypeError:
		return default


# --- keys --------------------------------------------------------------------


def measurement_key_name(value: Any) -> str | None:
	"""Return the measurement key name of ``value``, or ``None``."""
	subject = getattr(value, "gate", value)
	if hasattr(subject, "_measurement_key_name_"):
		return subject._measurement_key_name_()
	return None


def with_key_path(key: Any, path: Sequence[str]) -> MeasurementKey:
	"""Return ``key`` prefixed by ``path``."""
	name = str(key)
	for prefix in reversed(path):
		name = f"{prefix}.{name}"
	return MeasurementKey(name)


def with_key_path_prefix(key: Any, prefix: Sequence[str]) -> MeasurementKey:
	"""Return ``key`` prefixed by ``prefix``."""
	return with_key_path(key, prefix)


def with_rescoped_keys(key: Any, path: Sequence[str], bindable_keys: Iterable[Any] = ()) -> MeasurementKey:
	"""Return ``key`` rescoped by ``path``, unless it is bindable."""
	del bindable_keys
	return with_key_path(key, path)


# --- JSON --------------------------------------------------------------------


class HasJSONNamespace:
	"""Mixin for values that declare a JSON namespace."""

	def _json_namespace_(self) -> str:
		return ""


def json_namespace(value: Any) -> str:
	"""Return the JSON namespace of ``value``."""
	if hasattr(value, "_json_namespace_"):
		return value._json_namespace_()
	return ""


def dataclass_json_dict(obj: Any) -> dict[str, Any]:
	"""Return a JSON dictionary for a dataclass instance."""
	result: dict[str, Any] = {}
	if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
		for field in dataclasses.fields(obj):
			result[field.name] = getattr(obj, field.name)
	return result


def obj_to_dict_helper(obj: Any, attribute_names: Iterable[str], *, namespace: str | None = None) -> dict[str, Any]:
	"""Build a JSON dictionary from the named attributes of ``obj``."""
	result: dict[str, Any] = {name: getattr(obj, name) for name in attribute_names}
	if namespace is not None:
		result["cirq_type"] = namespace
	return result


class JsonResolver:
	"""Resolves ``cirq_type`` names to classes."""

	def __init__(self, resolvers: Mapping[str, type] | None = None) -> None:
		self._resolvers = dict(resolvers or {})

	def __getitem__(self, name: str) -> type:
		return self._resolvers[name]

	def __contains__(self, name: str) -> bool:
		return name in self._resolvers


class SerializableByKey:
	"""Base class for values serialized by a stable key."""

	def _serialization_key_(self) -> str:
		raise NotImplementedError


class LabelEntity:
	"""A labelled entity, used to tag diagram elements."""

	def __init__(self, *, name: str | None = None, display_name: str | None = None, description: str = "") -> None:
		self.name = name if name is not None else (display_name or "entity")
		self.display_name = display_name if display_name is not None else self.name
		self.description = description

	def __str__(self) -> str:
		return self.display_name

	def __repr__(self) -> str:
		return f"vqubit.LabelEntity(name={self.name!r})"


# --- trace distance ----------------------------------------------------------


def trace_distance_from_angle_list(angle_list: Sequence[float]) -> float:
	"""Return a trace-distance bound from a list of rotation angles."""
	if not angle_list:
		return 0.0
	real = sum(math.cos(angle) for angle in angle_list) / len(angle_list)
	imag = sum(math.sin(angle) for angle in angle_list) / len(angle_list)
	return math.sqrt(max(0.0, 1 - (real * real + imag * imag)))


# --- Supports* protocols -----------------------------------------------------


@runtime_checkable
class SupportsUnitary(Protocol):
	def _unitary_(self) -> Any: ...


@runtime_checkable
class SupportsExplicitHasUnitary(Protocol):
	def _has_unitary_(self) -> bool: ...


@runtime_checkable
class SupportsKraus(Protocol):
	def _kraus_(self) -> Any: ...


@runtime_checkable
class SupportsMixture(Protocol):
	def _mixture_(self) -> Any: ...


@runtime_checkable
class SupportsApplyUnitary(Protocol):
	def _apply_unitary_(self, args: ApplyUnitaryArgs) -> Any: ...


SupportsConsistentApplyUnitary = SupportsApplyUnitary


@runtime_checkable
class SupportsApplyChannel(Protocol):
	def _apply_channel_(self, args: ApplyChannelArgs) -> Any: ...


@runtime_checkable
class SupportsApplyMixture(Protocol):
	def _apply_mixture_(self, args: ApplyMixtureArgs) -> Any: ...


@runtime_checkable
class SupportsActOn(Protocol):
	def _act_on_(self, args: Any, qubits: Sequence[Qid] | None = None) -> Any: ...


SupportsActOnQubits = SupportsActOn


@runtime_checkable
class SupportsApproximateEquality(Protocol):
	def _approx_eq_(self, other: Any, atol: float) -> bool: ...


@runtime_checkable
class SupportsCircuitDiagramInfo(Protocol):
	def _circuit_diagram_info_(self, args: CircuitDiagramInfoArgs) -> Any: ...


@runtime_checkable
class SupportsCommutes(Protocol):
	def _commutes_(self, other: Any, *, atol: float) -> bool: ...


@runtime_checkable
class SupportsControlKey(Protocol):
	def _control_keys_(self) -> Any: ...


@runtime_checkable
class SupportsDecompose(Protocol):
	def _decompose_(self) -> Any: ...


@runtime_checkable
class SupportsDecomposeWithQubits(Protocol):
	def _decompose_(self, qubits: Sequence[Qid]) -> Any: ...


@runtime_checkable
class SupportsEqualUpToGlobalPhase(Protocol):
	def _equal_up_to_global_phase_(self, other: Any, *, atol: float) -> bool: ...


@runtime_checkable
class SupportsExplicitQidShape(Protocol):
	def _qid_shape_(self) -> tuple[int, ...]: ...


@runtime_checkable
class SupportsExplicitNumQubits(Protocol):
	def _num_qubits_(self) -> int: ...


@runtime_checkable
class SupportsJSON(Protocol):
	def _json_dict_(self) -> Any: ...


@runtime_checkable
class SupportsMeasurementKey(Protocol):
	def _measurement_key_obj_(self) -> Any: ...


@runtime_checkable
class SupportsParameterization(Protocol):
	def _is_parameterized_(self) -> bool: ...


@runtime_checkable
class SupportsPauliExpansion(Protocol):
	def _pauli_expansion_(self) -> Any: ...


@runtime_checkable
class SupportsPhase(Protocol):
	def _phase_by_(self, phase_turns: float, qubit_index: int) -> Any: ...


@runtime_checkable
class SupportsQasm(Protocol):
	def _qasm_(self) -> Any: ...


@runtime_checkable
class SupportsQasmWithArgs(Protocol):
	def _qasm_(self, args: Any) -> Any: ...


@runtime_checkable
class SupportsQasmWithArgsAndQubits(Protocol):
	def _qasm_(self, args: Any, qubits: Sequence[Qid]) -> Any: ...


@runtime_checkable
class SupportsTraceDistanceBound(Protocol):
	def _trace_distance_bound_(self) -> float: ...


__all__ = [
	"DEFAULT_RESOLVERS",
	"ApplyChannelArgs",
	"ApplyMixtureArgs",
	"ApplyUnitaryArgs",
	"CircuitDiagramInfo",
	"CircuitDiagramInfoArgs",
	"DecompositionContext",
	"HasJSONNamespace",
	"JsonResolver",
	"LabelEntity",
	"SerializableByKey",
	"SupportsActOn",
	"SupportsActOnQubits",
	"SupportsApplyChannel",
	"SupportsApplyMixture",
	"SupportsApproximateEquality",
	"SupportsCircuitDiagramInfo",
	"SupportsCommutes",
	"SupportsConsistentApplyUnitary",
	"SupportsControlKey",
	"SupportsDecompose",
	"SupportsDecomposeWithQubits",
	"SupportsEqualUpToGlobalPhase",
	"SupportsExplicitHasUnitary",
	"SupportsExplicitNumQubits",
	"SupportsExplicitQidShape",
	"SupportsJSON",
	"SupportsKraus",
	"SupportsMeasurementKey",
	"SupportsMixture",
	"SupportsParameterization",
	"SupportsPauliExpansion",
	"SupportsPhase",
	"SupportsQasm",
	"SupportsQasmWithArgs",
	"SupportsQasmWithArgsAndQubits",
	"SupportsTraceDistanceBound",
	"SupportsUnitary",
	"act_on",
	"apply_channel",
	"apply_mixture",
	"apply_unitaries",
	"circuit_diagram_info",
	"dataclass_json_dict",
	"decompose_once_with_qubits",
	"has_mixture",
	"has_stabilizer_effect",
	"json_namespace",
	"measurement_key_name",
	"mixture",
	"mul",
	"obj_to_dict_helper",
	"parameter_names",
	"parameter_symbols",
	"resolve_parameters_once",
	"trace_distance_from_angle_list",
	"validate_mixture",
	"with_key_path",
	"with_key_path_prefix",
	"with_rescoped_keys",
]
