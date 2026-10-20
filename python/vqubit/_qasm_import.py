"""OpenQASM 2.0 import: a small, dependency-free parser.

:func:`circuit_from_qasm` reads OpenQASM 2.0 programs — the ones :func:`vqubit.qasm`
writes, so the two round-trip, and a good part of what other tools emit: the
``qelib1`` gate set (plus Cirq's ``sx``/``sxdg``), custom ``gate`` definitions
(expanded inline), measurements into classical registers, ``reset``, ``barrier``
and classically controlled ``if`` statements.

Cirq's importer lives in ``cirq.contrib`` and needs the third-party ``ply``; this
one is a hand-written recursive-descent parser with no dependencies.

A few conventions, so that the round trip is exact:

* every declared qubit becomes a :class:`~vqubit.LineQubit`, in declaration
  order, so ``qreg q[3];`` is ``LineQubit.range(3)``;
* a classical register is a measurement key: ``creg c[2];`` means the key
  ``"c"``, and the measurements written into it are gathered into one
  ``MeasurementGate`` over their qubits, in register order;
* ``u3``/``u2``/``u1`` are written as ``rz``/``ry`` rotations (they are equal up
  to a global phase, which OpenQASM does not see), while their *controlled*
  forms (``cu3``, ``cu1``, ``crz``) use exact matrices, because there the phase
  is relative and does matter.
"""

from __future__ import annotations

import cmath
import math
import re
from collections.abc import Sequence
from typing import Any

from ._channels import reset
from ._circuits import Circuit, InsertStrategy
from ._classical import ClassicallyControlledOperation, KeyCondition
from ._gates import (
	CCX,
	CNOT,
	CSWAP,
	CZ,
	SWAP,
	ControlledGate,
	Gate,
	H,
	I,
	MatrixGate,
	RzGate,
	X,
	Y,
	Z,
	measure,
)
from ._ops import Operation
from ._parameter import Parameter
from ._qubits import LineQubit, Qid


class QasmParseError(ValueError):
	"""A malformed or unsupported OpenQASM 2.0 program, with its line number."""

	def __init__(self, message: str, line: int | None = None) -> None:
		self.line = line
		super().__init__(f"line {line}: {message}" if line is not None else message)


# --- tokens ------------------------------------------------------------------

_TOKEN = re.compile(
	r"""
	(?P<space>[ \t\r]+)
	|(?P<comment>//[^\n]*|/\*.*?\*/)
	|(?P<newline>\n)
	|(?P<number>(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?)
	|(?P<name>[A-Za-z_][A-Za-z0-9_]*)
	|(?P<string>"[^"]*")
	|(?P<arrow>->)
	|(?P<op>==|&&|[-+*/^(),\[\]{};=])
	""",
	re.VERBOSE | re.DOTALL,
)

Token = tuple[str, str, int]


def _tokenize(source: str) -> list[Token]:
	"""Return the ``(kind, text, line)`` tokens, without blanks or comments."""
	tokens: list[Token] = []
	position = 0
	line = 1
	while position < len(source):
		match = _TOKEN.match(source, position)
		if match is None:
			raise QasmParseError(f"unexpected character {source[position]!r}", line)
		kind = match.lastgroup or ""
		text = match.group()
		position = match.end()
		if kind == "newline":
			line += 1
			continue
		if kind in ("space", "comment"):
			line += text.count("\n")
			continue
		tokens.append((kind, text, line))
	return tokens


# --- gates -------------------------------------------------------------------


def _u3_matrix(theta: Any, phi: Any, lam: Any) -> list[list[complex]]:
	"""The exact ``u3`` matrix, phase included (its controlled form needs it)."""
	cosine, sine = math.cos(float(theta) / 2), math.sin(float(theta) / 2)
	return [
		[complex(cosine), -cmath.exp(1j * float(lam)) * sine],
		[cmath.exp(1j * float(phi)) * sine, cmath.exp(1j * (float(phi) + float(lam))) * cosine],
	]


def _u1_matrix(lam: Any) -> list[list[complex]]:
	return [[1 + 0j, 0j], [0j, cmath.exp(1j * float(lam))]]


def _plain_gates() -> dict[str, tuple[Gate, int]]:
	"""Gates whose name says it all, with the number of qubits they take."""
	return {
		"id": (I, 1),
		"x": (X, 1),
		"y": (Y, 1),
		"z": (Z, 1),
		"h": (H, 1),
		"s": (Z**0.5, 1),
		"sdg": (Z**-0.5, 1),
		"t": (Z**0.25, 1),
		"tdg": (Z**-0.25, 1),
		"sx": (X**0.5, 1),
		"sxdg": (X**-0.5, 1),
		"cx": (CNOT, 2),
		"CX": (CNOT, 2),
		"cz": (CZ, 2),
		"swap": (SWAP, 2),
		"ccx": (CCX, 3),
		"cswap": (CSWAP, 3),
		"cy": (ControlledGate(Y, num_controls=1), 2),
		"ch": (ControlledGate(H, num_controls=1), 2),
	}


# --- parser ------------------------------------------------------------------


class _Parser:
	"""A recursive-descent parser for the part of OpenQASM 2.0 we support."""

	def __init__(self, tokens: Sequence[Token]) -> None:
		self._tokens = list(tokens)
		self._index = 0
		self._qubits: list[Qid] = []
		self._qregs: dict[str, list[Qid]] = {}
		self._cregs: dict[str, int] = {}
		self._gates: dict[str, tuple[list[str], list[str], list[Token]]] = {}
		self._opaque: set[str] = set()
		self._pending_measurements: list[tuple[str, Qid]] = []
		self._operations: list[Operation] = []
		self._scope_qubits: dict[str, Any] = {}
		self._scope_params: dict[str, Any] = {}
		self._plain = _plain_gates()

	# -- token helpers --
	def _peek(self, offset: int = 0) -> Token:
		index = self._index + offset
		if index >= len(self._tokens):
			line = self._tokens[-1][2] if self._tokens else 1
			return ("eof", "", line)
		return self._tokens[index]

	def _next(self) -> Token:
		token = self._peek()
		self._index += 1
		return token

	def _accept(self, text: str) -> bool:
		if self._peek()[1] == text:
			self._index += 1
			return True
		return False

	def _expect(self, text: str) -> Token:
		token = self._peek()
		if token[1] != text:
			raise QasmParseError(f"expected {text!r}, found {token[1]!r}", token[2])
		self._index += 1
		return token

	def _expect_name(self) -> tuple[str, int]:
		kind, text, line = self._peek()
		if kind != "name":
			raise QasmParseError(f"expected a name, found {text!r}", line)
		self._index += 1
		return text, line

	def _expect_integer(self) -> int:
		kind, text, line = self._next()
		if kind != "number":
			raise QasmParseError(f"expected a number, found {text!r}", line)
		value = float(text)
		if value != int(value):
			raise QasmParseError(f"expected an integer, found {text!r}", line)
		return int(value)

	# -- program --
	def parse(self) -> Circuit:
		while self._peek()[0] != "eof":
			self._statement()
		self._flush_measurements()
		return Circuit(*self._operations, strategy=InsertStrategy.NEW)

	def _statement(self) -> None:
		kind, text, line = self._peek()
		if kind != "name":
			raise QasmParseError(f"unexpected token {text!r}", line)
		handlers = {
			"OPENQASM": self._header,
			"include": self._include,
			"qreg": self._declare_qreg,
			"creg": self._declare_creg,
			"gate": self._define_gate,
			"opaque": self._declare_opaque,
			"barrier": self._barrier,
			"measure": self._measure,
			"reset": self._reset,
			"if": self._conditional,
		}
		handler = handlers.get(text)
		if handler is not None:
			self._next()
			handler()
			return
		if self._peek(1)[1] == "=" or (self._peek(1)[1] == "[" and self._peek(4)[1] == "="):
			self._classical_measure()
			return
		self._gate_application()

	def _header(self) -> None:
		self._expect_name_is_version()
		self._expect(";")

	def _expect_name_is_version(self) -> None:
		kind, text, line = self._next()
		if kind != "number" or text != "2.0":
			raise QasmParseError(f"expected OpenQASM version 2.0, found {text!r}", line)

	def _include(self) -> None:
		kind, text, line = self._next()
		if kind != "string":
			raise QasmParseError(f"expected a file name after include, found {text!r}", line)
		self._expect(";")

	def _barrier(self) -> None:
		while self._peek()[1] != ";":
			if self._peek()[0] == "eof":
				raise QasmParseError("unterminated barrier")
			self._next()
		self._expect(";")

	# -- declarations --
	def _declare_qreg(self) -> None:
		name, line = self._expect_name()
		self._expect("[")
		size = self._expect_integer()
		self._expect("]")
		self._expect(";")
		if name in self._qregs:
			raise QasmParseError(f"quantum register {name!r} is already declared", line)
		register: list[Qid] = [LineQubit(len(self._qubits) + index) for index in range(size)]
		self._qubits.extend(register)
		self._qregs[name] = register

	def _declare_creg(self) -> None:
		name, line = self._expect_name()
		self._expect("[")
		size = self._expect_integer()
		self._expect("]")
		self._expect(";")
		if name in self._cregs:
			raise QasmParseError(f"classical register {name!r} is already declared", line)
		self._cregs[name] = size

	def _declare_opaque(self) -> None:
		name, _line = self._expect_name()
		if self._peek()[1] == "(":
			self._skip_balanced()
		while self._peek()[1] != ";":
			self._next()
		self._expect(";")
		self._opaque.add(name)

	def _skip_balanced(self) -> None:
		depth = 0
		while True:
			kind, text, line = self._next()
			if kind == "eof":
				raise QasmParseError("unterminated parenthesis", line)
			if text == "(":
				depth += 1
			elif text == ")":
				depth -= 1
				if depth == 0:
					return

	def _define_gate(self) -> None:
		name, line = self._expect_name()
		params: list[str] = []
		if self._peek()[1] == "(":
			params = self._name_list("(", ")")
		qargs = self._name_list(None, None, terminator="{")
		self._expect("{")
		start = self._index
		depth = 1
		while depth:
			kind, text, _body_line = self._next()
			if kind == "eof":
				raise QasmParseError(f"unterminated definition of gate {name!r}", line)
			if text == "{":
				depth += 1
			elif text == "}":
				depth -= 1
		self._gates[name] = (params, qargs, self._tokens[start : self._index - 1])

	def _name_list(self, opening: str | None, closing: str | None, *, terminator: str | None = None) -> list[str]:
		"""Read ``a, b, c`` (optionally inside brackets, optionally up to ``{``)."""
		if opening is not None:
			self._expect(opening)
		names: list[str] = []
		while True:
			name, _line = self._expect_name()
			names.append(name)
			if not self._accept(","):
				break
		if closing is not None:
			self._expect(closing)
		if terminator is not None:
			while self._peek()[1] != terminator:
				if self._peek()[0] == "eof":
					raise QasmParseError(f"expected {terminator!r}")
				self._next()
		return names

	# -- arguments --
	def _expression_list(self) -> list[Any]:
		self._expect("(")
		values = [self._expression()]
		while self._accept(","):
			values.append(self._expression())
		self._expect(")")
		return values

	def _qubit_args(self) -> list[list[Any]]:
		"""Read the qubit arguments of a call, one list per argument."""
		args: list[list[Any]] = []
		while True:
			name, line = self._expect_name()
			args.append(self._resolve_qubit_name(name, line))
			if not self._accept(","):
				break
		return args

	def _qubit_arg(self) -> list[Any]:
		name, line = self._expect_name()
		return self._resolve_qubit_name(name, line)

	def _resolve_qubit_name(self, name: str, line: int) -> list[Any]:
		if name in self._scope_qubits:
			return [self._scope_qubits[name]]
		if name not in self._qregs:
			raise QasmParseError(f"undefined quantum register {name!r}", line)
		register = self._qregs[name]
		if self._accept("["):
			index = self._expect_integer()
			self._expect("]")
			if not 0 <= index < len(register):
				raise QasmParseError(f"{name}[{index}] is out of range", line)
			return [register[index]]
		return list(register)

	def _creg_arg(self) -> tuple[str, int | None]:
		name, line = self._expect_name()
		if name not in self._cregs:
			raise QasmParseError(f"undefined classical register {name!r}", line)
		if self._accept("["):
			index = self._expect_integer()
			self._expect("]")
			if not 0 <= index < self._cregs[name]:
				raise QasmParseError(f"{name}[{index}] is out of range", line)
			return name, index
		return name, None

	# -- expressions --
	def _expression(self) -> Any:
		value = self._term()
		while self._peek()[1] in ("+", "-"):
			operator = self._next()[1]
			right = self._term()
			value = value + right if operator == "+" else value - right
		return value

	def _term(self) -> Any:
		value = self._unary()
		while self._peek()[1] in ("*", "/"):
			operator = self._next()[1]
			right = self._unary()
			value = value * right if operator == "*" else value / right
		return value

	def _unary(self) -> Any:
		if self._accept("-"):
			return -self._unary()
		if self._accept("+"):
			return self._unary()
		return self._power()

	def _power(self) -> Any:
		value = self._atom()
		if self._accept("^"):
			value = value ** self._unary()
		return value

	def _atom(self) -> Any:
		kind, text, line = self._next()
		if kind == "number":
			number = float(text)
			return int(number) if number == int(number) else number
		if kind == "name":
			if text == "pi":
				return math.pi
			if text in self._scope_params:
				return self._scope_params[text]
			return Parameter(text)
		if text == "(":
			value = self._expression()
			self._expect(")")
			return value
		raise QasmParseError(f"unexpected token {text!r} in an expression", line)

	# -- statements --
	def _gate_application(self) -> None:
		name, line = self._expect_name()
		params = self._expression_list() if self._peek()[1] == "(" else []
		args = self._qubit_args()
		self._expect(";")
		self._append(self._build_operations(name, params, args, line))

	def _build_operations(self, name: str, params: list[Any], args: list[list[Any]], line: int) -> list[Operation]:
		"""Apply ``name(params)`` to ``args``, expanding registers elementwise."""
		if name in self._opaque:
			raise QasmParseError(f"gate {name!r} is opaque and cannot be expanded", line)
		if name in self._gates:
			return self._expand_custom_gate(name, params, args, line)
		plain = self._plain.get(name)
		if plain is not None:
			gate, num_qubits = plain
			if params:
				raise QasmParseError(f"gate {name!r} takes no parameters", line)
			return self._apply(gate, args, num_qubits, line)
		builders: dict[str, tuple[int, Any]] = {
			"u3": (3, lambda p: [_rz(p[2]), _ry(p[0]), _rz(p[1])]),
			"U": (3, lambda p: [_rz(p[2]), _ry(p[0]), _rz(p[1])]),
			"u2": (2, lambda p: [_rz(p[1]), _ry(math.pi / 2), _rz(p[0])]),
			"u1": (1, lambda p: [_rz(p[0])]),
			"u0": (1, lambda p: [I]),
			"rx": (1, lambda p: [_rx(p[0])]),
			"ry": (1, lambda p: [_ry(p[0])]),
			"rz": (1, lambda p: [_rz(p[0])]),
			"crz": (1, lambda p: [ControlledGate(RzGate(p[0]), num_controls=1)]),
			"cu1": (1, lambda p: [ControlledGate(MatrixGate(_u1_matrix(p[0])), num_controls=1)]),
			"cu3": (3, lambda p: [ControlledGate(MatrixGate(_u3_matrix(p[0], p[1], p[2])), num_controls=1)]),
		}
		if name not in builders:
			raise QasmParseError(f"unknown gate {name!r}", line)
		arity, builder = builders[name]
		if len(params) != arity:
			raise QasmParseError(f"gate {name!r} takes {arity} parameter(s), got {len(params)}", line)
		gates = builder(params)
		if len(gates) == 1:
			# A builder may return a controlled gate: ask the gate itself how
			# many qubits it acts on instead of assuming one.
			return self._apply(gates[0], args, gates[0].num_qubits, line)
		# A multi-gate expansion (``u3``) applies to one qubit at a time.
		operations: list[Operation] = []
		for (qubit,) in self._pairs(args, 1, line):
			operations.extend(gate.on(qubit) for gate in gates)
		return operations

	def _pairs(self, args: list[list[Any]], num_qubits: int, line: int) -> list[tuple[Any, ...]]:
		if len(args) != num_qubits:
			raise QasmParseError(f"expected {num_qubits} qubit argument(s), got {len(args)}", line)
		count = len(args[0])
		if any(len(arg) != count for arg in args):
			raise QasmParseError("the qubit arguments have different register sizes", line)
		return [tuple(arg[index] for arg in args) for index in range(count)]

	def _apply(self, gate: Gate, args: list[list[Any]], num_qubits: int, line: int) -> list[Operation]:
		return [gate.on(*qubits) for qubits in self._pairs(args, num_qubits, line)]

	def _expand_custom_gate(self, name: str, params: list[Any], args: list[list[Any]], line: int) -> list[Operation]:
		param_names, qubit_names, body = self._gates[name]
		if len(params) != len(param_names):
			raise QasmParseError(f"gate {name!r} takes {len(param_names)} parameter(s), got {len(params)}", line)
		pairs = self._pairs(args, len(qubit_names), line)
		bindings = dict(zip(param_names, params, strict=True))
		operations: list[Operation] = []
		for pair in pairs:
			sub_parser = _Parser([*list(body), ("eof", "", line)])
			sub_parser._qregs, sub_parser._cregs = self._qregs, self._cregs
			sub_parser._gates, sub_parser._opaque = self._gates, self._opaque
			sub_parser._qubits = self._qubits
			sub_parser._scope_qubits = dict(zip(qubit_names, pair, strict=True))
			sub_parser._scope_params = bindings
			while sub_parser._peek()[0] != "eof":
				sub_parser._statement()
			sub_parser._flush_measurements()
			operations.extend(sub_parser._operations)
		return operations

	def _measure(self) -> None:
		qubits = self._qubit_arg()
		self._expect("->")
		key, _index = self._creg_arg()
		self._expect(";")
		for qubit in qubits:
			self._pending_measurements.append((key, qubit))

	def _classical_measure(self) -> None:
		key, _index = self._creg_arg()
		self._expect("=")
		self._expect("measure")
		qubits = self._qubit_arg()
		self._expect(";")
		for qubit in qubits:
			self._pending_measurements.append((key, qubit))

	def _reset(self) -> None:
		qubits = self._qubit_arg()
		self._expect(";")
		self._append([reset(qubit) for qubit in qubits])

	def _conditional(self) -> None:
		self._expect("(")
		conditions: list[Any] = []
		while True:
			key, index = self._creg_arg()
			self._expect("==")
			value = self._expect_integer()
			positions = [index] if index is not None else list(range(self._cregs[key]))
			for position in positions:
				bit = value if index is not None else (value >> position) & 1
				condition = KeyCondition(key, position)
				conditions.append(condition if bit else ~condition)
			if not self._accept("&&"):
				break
		self._expect(")")
		self._flush_measurements()
		before = len(self._operations)
		self._gate_application()
		controlled = self._operations[before:]
		del self._operations[before:]
		self._append([ClassicallyControlledOperation(conditions, operation) for operation in controlled])

	# -- state --
	def _append(self, operations: Sequence[Operation]) -> None:
		self._flush_measurements()
		self._operations.extend(operations)

	def _flush_measurements(self) -> None:
		"""Turn the pending measurements into one gate per classical register."""
		if not self._pending_measurements:
			return
		pending, self._pending_measurements = self._pending_measurements, []
		keys = {key for key, _qubit in pending}
		if len(keys) != 1:
			raise QasmParseError("measurements into different registers need a gate between them")
		key = keys.pop()
		self._operations.append(measure(*(qubit for _key, qubit in pending), key=key))


def _rx(angle: Any) -> Gate:
	from ._gates import RxGate

	return RxGate(angle)


def _ry(angle: Any) -> Gate:
	from ._gates import RyGate

	return RyGate(angle)


def _rz(angle: Any) -> Gate:
	return RzGate(angle)


def circuit_from_qasm(source: str) -> Circuit:
	"""Parse an OpenQASM 2.0 program into a :class:`~vqubit.Circuit`.

	See the module docstring for the conventions (qubits, measurement keys and
	the ``u3`` family). Raises :class:`QasmParseError`, with the line number, on
	a malformed or unsupported program.
	"""
	if not isinstance(source, str):
		raise TypeError(f"circuit_from_qasm needs the QASM text, got {type(source).__name__}")
	return _Parser(_tokenize(source)).parse()


__all__ = ["QasmParseError", "circuit_from_qasm"]
