"""JSON serialization, compatible with Cirq's ``cirq_type`` convention.

Every serializable value becomes a JSON object carrying a ``cirq_type`` field
naming its class, plus its defining fields. :func:`read_json` rebuilds the
value from that field through a registry, so a circuit survives a round trip
through text.
"""

from __future__ import annotations

import gzip
import json
from typing import Any

from ._circuits import Circuit, FrozenCircuit
from ._gates import (
	CCXGate,
	CNOTGate,
	ControlledGate,
	CSWAPGate,
	CZGate,
	HPowGate,
	IdentityGate,
	ISWAPGate,
	MeasurementGate,
	RxGate,
	RyGate,
	RzGate,
	SWAPGate,
	XPowGate,
	XXPowGate,
	YPowGate,
	YYPowGate,
	ZPowGate,
	ZZPowGate,
)
from ._ops import GateOperation, Moment
from ._parameter import Parameter, _BinaryOp, _UnaryOp
from ._pauli import Pauli, PauliString
from ._qubits import GridQubit, LineQubit, NamedQubit

_PAULI_BY_VALUE = {p.value: p for p in Pauli}

# Decoders keyed by `cirq_type`. Each takes the raw JSON dict (with the
# `cirq_type` field already removed) and returns the value.
_DECODERS: dict[str, Any] = {}


def _decoder(cirq_type: str) -> Any:
	def register(fn: Any) -> Any:
		_DECODERS[cirq_type] = fn
		return fn

	return register


def _encode(value: Any) -> Any:
	if value is None or isinstance(value, (bool, int, float, str)):
		return value
	if isinstance(value, complex):
		return {"cirq_type": "complex", "real": value.real, "imag": value.imag}
	if isinstance(value, (list, tuple)):
		return [_encode(item) for item in value]
	if isinstance(value, dict):
		return {key: _encode(item) for key, item in value.items()}
	if isinstance(value, Parameter):
		return {"cirq_type": "Parameter", "name": value.name}
	if isinstance(value, _BinaryOp):
		return {"cirq_type": "BinaryOp", "op": value._op, "left": _encode(value._left), "right": _encode(value._right)}
	if isinstance(value, _UnaryOp):
		return {"cirq_type": "UnaryOp", "op": value._op, "operand": _encode(value._operand)}
	if isinstance(value, LineQubit):
		return {"cirq_type": "LineQubit", "x": value.x}
	if isinstance(value, GridQubit):
		return {"cirq_type": "GridQubit", "row": value.row, "col": value.col}
	if isinstance(value, NamedQubit):
		return {"cirq_type": "NamedQubit", "name": value.name}
	if isinstance(value, (XPowGate, YPowGate, ZPowGate, HPowGate, XXPowGate, YYPowGate, ZZPowGate)):
		return {
			"cirq_type": type(value).__name__,
			"exponent": _encode(value.exponent),
			"global_shift": _encode(getattr(value, "global_shift", 0.0)),
		}
	if isinstance(value, (RxGate, RyGate, RzGate)):
		return {"cirq_type": type(value).__name__, "radians": _encode(value.radians)}
	if isinstance(value, IdentityGate):
		return {"cirq_type": "IdentityGate", "num_qubits": value.num_qubits}
	if isinstance(value, MeasurementGate):
		return {
			"cirq_type": "MeasurementGate",
			"key": value.key,
			"num_qubits": value.num_qubits,
			"invert_mask": list(value.invert_mask),
		}
	if isinstance(value, ControlledGate):
		return {
			"cirq_type": "ControlledGate",
			"sub_gate": _encode(value.sub_gate),
			"num_controls": value.num_controls,
			"control_values": list(value.control_values),
		}
	if isinstance(value, (CNOTGate, CZGate, SWAPGate, CCXGate, CSWAPGate, ISWAPGate)):
		return {"cirq_type": type(value).__name__}
	if isinstance(value, GateOperation):
		return {"cirq_type": "GateOperation", "gate": _encode(value.gate), "qubits": _encode(list(value.qubits))}
	if isinstance(value, Moment):
		return {"cirq_type": "Moment", "operations": _encode(list(value.operations))}
	if isinstance(value, (Circuit, FrozenCircuit)):
		cirq_type = "FrozenCircuit" if isinstance(value, FrozenCircuit) else "Circuit"
		return {"cirq_type": cirq_type, "moments": [_encode(moment) for moment in value.moments]}
	if isinstance(value, PauliString):
		return {
			"cirq_type": "PauliString",
			"coefficient": _encode(value.coefficient),
			"qubit_pauli_map": [[_encode(qubit), pauli.value] for qubit, pauli in value.items()],
		}
	raise TypeError(f"cannot serialize {value!r} to JSON")


_TYPES: dict[str, Any] = {
	"LineQubit": LineQubit,
	"GridQubit": GridQubit,
	"NamedQubit": NamedQubit,
	"Parameter": Parameter,
	"XPowGate": XPowGate,
	"YPowGate": YPowGate,
	"ZPowGate": ZPowGate,
	"HPowGate": HPowGate,
	"XXPowGate": XXPowGate,
	"YYPowGate": YYPowGate,
	"ZZPowGate": ZZPowGate,
	"RxGate": RxGate,
	"RyGate": RyGate,
	"RzGate": RzGate,
	"IdentityGate": IdentityGate,
	"MeasurementGate": MeasurementGate,
	"ControlledGate": ControlledGate,
	"CNOTGate": CNOTGate,
	"CZGate": CZGate,
	"SWAPGate": SWAPGate,
	"CCXGate": CCXGate,
	"CSWAPGate": CSWAPGate,
	"ISWAPGate": ISWAPGate,
	"GateOperation": GateOperation,
	"Moment": Moment,
	"Circuit": Circuit,
	"FrozenCircuit": FrozenCircuit,
	"PauliString": PauliString,
	"complex": complex,
}


def cirq_type_from_json(cirq_type: str) -> Any:
	"""Return the class registered for ``cirq_type``, or ``None``."""
	return _TYPES.get(cirq_type)


def json_cirq_type(value: Any) -> str:
	"""Return the ``cirq_type`` name for a value or class."""
	return value.__name__ if isinstance(value, type) else type(value).__name__


def to_json(value: Any, *, indent: int | None = None) -> str:
	"""Return ``value`` as a JSON string."""
	return json.dumps(_encode(value), indent=indent)


def to_json_gzip(value: Any) -> bytes:
	"""Return ``value`` as gzip-compressed JSON bytes."""
	return gzip.compress(to_json(value).encode("utf-8"))


def _decode(data: Any) -> Any:
	if isinstance(data, list):
		return [_decode(item) for item in data]
	if not isinstance(data, dict):
		return data
	cirq_type = data.get("cirq_type")
	if cirq_type is None:
		return {key: _decode(item) for key, item in data.items()}
	if cirq_type == "complex":
		return complex(data["real"], data["imag"])
	decoder = _DECODERS.get(cirq_type)
	if decoder is None:
		raise ValueError(f"no decoder registered for {cirq_type!r}")
	return decoder(data)


def read_json(text: str) -> Any:
	"""Rebuild a value from its JSON representation."""
	return _decode(json.loads(text))


def read_json_gzip(data: bytes) -> Any:
	"""Rebuild a value from gzip-compressed JSON bytes."""
	return read_json(gzip.decompress(data).decode("utf-8"))


def _decode_pauli_string(data: dict[str, Any]) -> PauliString:
	qubit_pauli_map = {_decode(qubit): _PAULI_BY_VALUE[value] for qubit, value in data["qubit_pauli_map"]}
	return PauliString(qubit_pauli_map, coefficient=_decode(data["coefficient"]))


def _register_all() -> None:
	_DECODERS["LineQubit"] = lambda d: LineQubit(d["x"])
	_DECODERS["GridQubit"] = lambda d: GridQubit(d["row"], d["col"])
	_DECODERS["NamedQubit"] = lambda d: NamedQubit(d["name"])
	_DECODERS["Parameter"] = lambda d: Parameter(d["name"])
	_DECODERS["BinaryOp"] = lambda d: _BinaryOp(d["op"], _decode(d["left"]), _decode(d["right"]))
	_DECODERS["UnaryOp"] = lambda d: _UnaryOp(d["op"], _decode(d["operand"]))
	_DECODERS["XPowGate"] = lambda d: XPowGate(
		exponent=_decode(d["exponent"]), global_shift=_decode(d.get("global_shift", 0.0))
	)
	_DECODERS["YPowGate"] = lambda d: YPowGate(
		exponent=_decode(d["exponent"]), global_shift=_decode(d.get("global_shift", 0.0))
	)
	_DECODERS["ZPowGate"] = lambda d: ZPowGate(
		exponent=_decode(d["exponent"]), global_shift=_decode(d.get("global_shift", 0.0))
	)
	_DECODERS["HPowGate"] = lambda d: HPowGate(
		exponent=_decode(d["exponent"]), global_shift=_decode(d.get("global_shift", 0.0))
	)
	_DECODERS["XXPowGate"] = lambda d: XXPowGate(exponent=_decode(d["exponent"]))
	_DECODERS["YYPowGate"] = lambda d: YYPowGate(exponent=_decode(d["exponent"]))
	_DECODERS["ZZPowGate"] = lambda d: ZZPowGate(exponent=_decode(d["exponent"]))
	_DECODERS["RxGate"] = lambda d: RxGate(_decode(d["radians"]))
	_DECODERS["RyGate"] = lambda d: RyGate(_decode(d["radians"]))
	_DECODERS["RzGate"] = lambda d: RzGate(_decode(d["radians"]))
	_DECODERS["IdentityGate"] = lambda d: IdentityGate(d["num_qubits"])
	_DECODERS["MeasurementGate"] = lambda d: MeasurementGate(
		d["key"], num_qubits=d["num_qubits"], invert_mask=tuple(d.get("invert_mask", ()))
	)
	_DECODERS["ControlledGate"] = lambda d: ControlledGate(
		_decode(d["sub_gate"]), num_controls=d["num_controls"], control_values=tuple(d["control_values"])
	)
	_DECODERS["CNOTGate"] = lambda d: CNOTGate()
	_DECODERS["CZGate"] = lambda d: CZGate()
	_DECODERS["SWAPGate"] = lambda d: SWAPGate()
	_DECODERS["CCXGate"] = lambda d: CCXGate()
	_DECODERS["CSWAPGate"] = lambda d: CSWAPGate()
	_DECODERS["ISWAPGate"] = lambda d: ISWAPGate()
	_DECODERS["GateOperation"] = lambda d: GateOperation(_decode(d["gate"]), [_decode(q) for q in d["qubits"]])
	_DECODERS["Moment"] = lambda d: Moment([_decode(op) for op in d["operations"]])
	_DECODERS["Circuit"] = lambda d: Circuit.from_moments([_decode(m) for m in d["moments"]])
	_DECODERS["FrozenCircuit"] = lambda d: FrozenCircuit.from_moments([_decode(m) for m in d["moments"]])
	_DECODERS["PauliString"] = _decode_pauli_string


_register_all()
