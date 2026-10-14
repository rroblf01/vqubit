"""vqubit: a dependency-free quantum computing framework for Python.

vqubit lets you build quantum circuits and run them on local simulators, with
no third-party dependencies at all. The numerically heavy parts -- linear
algebra, gate matrices and simulation -- are written in V and compiled into the
native extension ``vqubit_native``; this package is the public, dependency-free
Python API on top of it.

The public API mirrors Cirq closely::

    import vqubit as vq

    q0, q1 = vq.LineQubit.range(2)
    circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1))
    print(circuit)

    result = vq.Simulator().run(circuit, repetitions=100)
    print(result.histogram(key="result"))
"""

from __future__ import annotations

from vqubit_native import version as _native_version

from ._circuits import Circuit, FrozenCircuit, InsertStrategy
from ._gates import (
	CCX,
	CNOT,
	CSWAP,
	CX,
	CZ,
	FREDKIN,
	ISWAP,
	SWAP,
	TOFFOLI,
	XX,
	YY,
	ZZ,
	CCXGate,
	CNOTGate,
	ControlledGate,
	CSWAPGate,
	CZGate,
	EigenGate,
	Gate,
	H,
	HPowGate,
	I,
	IdentityGate,
	ISWAPGate,
	MeasurementGate,
	Rx,
	RxGate,
	Ry,
	RyGate,
	Rz,
	RzGate,
	S,
	SWAPGate,
	T,
	X,
	XPowGate,
	XXPowGate,
	Y,
	YPowGate,
	YYPowGate,
	Z,
	ZPowGate,
	ZZPowGate,
	measure,
	measure_each,
)
from ._ops import GateOperation, Moment, Operation
from ._parameter import Parameter, ParamResolver
from ._pauli import Pauli, PauliString
from ._protocols import (
	commutes,
	decompose,
	definitely_commutes,
	has_unitary,
	inverse,
	is_measurement,
	is_parameterized,
	measurement_key_names,
	measurement_keys_touched,
	num_qubits,
	pow,
	qid_shape,
	resolve_parameters,
	unitary,
)
from ._qubits import GridQubit, LineQubit, NamedQubit, Qid, Qubit
from ._sim import Result, SimulationResult, Simulator, StateVector
from ._value import value_equality

__version__ = _native_version()

__all__ = [
	"CCX",
	"CNOT",
	"CSWAP",
	"CX",
	"CZ",
	"FREDKIN",
	"ISWAP",
	"SWAP",
	"TOFFOLI",
	"XX",
	"YY",
	"ZZ",
	"CCXGate",
	"CNOTGate",
	"CSWAPGate",
	"CZGate",
	"Circuit",
	"ControlledGate",
	"EigenGate",
	"FrozenCircuit",
	"Gate",
	"GateOperation",
	"GridQubit",
	"H",
	"HPowGate",
	"I",
	"ISWAPGate",
	"IdentityGate",
	"InsertStrategy",
	"LineQubit",
	"MeasurementGate",
	"Moment",
	"NamedQubit",
	"Operation",
	"ParamResolver",
	"Parameter",
	"Pauli",
	"PauliString",
	"Qid",
	"Qubit",
	"Result",
	"Rx",
	"RxGate",
	"Ry",
	"RyGate",
	"Rz",
	"RzGate",
	"S",
	"SWAPGate",
	"SimulationResult",
	"Simulator",
	"StateVector",
	"T",
	"X",
	"XPowGate",
	"XXPowGate",
	"Y",
	"YPowGate",
	"YYPowGate",
	"Z",
	"ZPowGate",
	"ZZPowGate",
	"__version__",
	"commutes",
	"decompose",
	"definitely_commutes",
	"has_unitary",
	"inverse",
	"is_measurement",
	"is_parameterized",
	"measure",
	"measure_each",
	"measurement_key_names",
	"measurement_keys_touched",
	"num_qubits",
	"pow",
	"qid_shape",
	"resolve_parameters",
	"unitary",
	"value_equality",
]
