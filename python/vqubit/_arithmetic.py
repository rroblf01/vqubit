"""Reversible classical arithmetic and Boolean-Hamiltonian gates."""

from __future__ import annotations

import ast
import math
import operator
from collections.abc import Iterable, Sequence
from itertools import product
from typing import Any

from ._gates import Gate
from ._linalg import Matrix
from ._value import value_equality

_BINOP: dict[type, Any] = {ast.BitAnd: operator.and_, ast.BitOr: operator.or_, ast.BitXor: operator.xor}


class ArithmeticGate(Gate):
	"""Base class for reversible classical arithmetic on registers of qubits.

	Subclasses override :meth:`registers`, :meth:`with_registers` and
	:meth:`apply`; the base class turns ``apply`` into a permutation matrix.
	Registers are big-endian, and a register may be a classical constant (an
	``int``) that ``apply`` must leave unchanged.
	"""

	def registers(self) -> Sequence[int | Sequence[int]]:
		"""The registers the gate acts on (constants or qudit-dimension lists)."""
		raise NotImplementedError

	def with_registers(self, *new_registers: int | Sequence[int]) -> ArithmeticGate:
		"""Return the same gate targeting different registers."""
		raise NotImplementedError

	def apply(self, *register_values: int) -> int | Iterable[int]:
		"""Return the new register values (see Cirq's contract for omissions)."""
		raise NotImplementedError

	def _num_qubits_(self) -> int:
		return sum(len(r) for r in self.registers() if not isinstance(r, int))

	def _unitary_(self) -> Matrix:
		registers = self.registers()
		sizes: list[int] = []
		ranges: list[Any] = []
		for register in registers:
			if isinstance(register, int):
				sizes.append(1)
				ranges.append([register])
			else:
				size = math.prod(register)
				sizes.append(size)
				ranges.append(range(size))
		dim = math.prod(sizes)

		def flat(indices: Sequence[int]) -> int:
			value = 0
			for index, size in zip(indices, sizes, strict=True):
				value = value * size + index
			return value

		matrix = [[0j] * dim for _ in range(dim)]
		for combo in product(*ranges):
			output = self.apply(*combo)
			outputs = [output] if isinstance(output, int) else list(output)
			if len(outputs) < len(combo):
				outputs += list(combo[len(outputs) :])
			out_indices: list[int] = []
			in_indices: list[int] = []
			for value, size, register in zip(combo, sizes, registers, strict=True):
				if isinstance(register, int):
					out_indices.append(0)
					in_indices.append(0)
				else:
					out_indices.append(int(outputs[len(out_indices)]) % size)
					in_indices.append(value)
			matrix[flat(out_indices)][flat(in_indices)] = 1 + 0j
		return matrix

	def _has_unitary_(self) -> bool:
		return True

	def __str__(self) -> str:
		return f"{type(self).__name__}({self.registers()})"

	def __repr__(self) -> str:
		return f"vqubit.{type(self).__name__}({self.registers()!r})"


def _parse_boolean(expression: str, names: Sequence[str]) -> ast.expr:
	tree = ast.parse(expression, mode="eval")
	_validate_boolean(tree.body, set(names))
	return tree.body


def _validate_boolean(node: ast.AST, names: set[str]) -> None:
	if isinstance(node, ast.Name):
		if node.id not in names:
			raise ValueError(f"Unknown symbol {node.id!r} in boolean expression")
		return
	if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Invert):
		_validate_boolean(node.operand, names)
		return
	if isinstance(node, ast.BinOp) and type(node.op) in _BINOP:
		_validate_boolean(node.left, names)
		_validate_boolean(node.right, names)
		return
	raise ValueError("Boolean expressions may only use ~, &, | and ^ over the parameter names")


def _eval_boolean(node: ast.AST, env: dict[str, bool]) -> bool:
	if isinstance(node, ast.Name):
		return env[node.id]
	if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Invert):
		return not _eval_boolean(node.operand, env)
	if isinstance(node, ast.BinOp) and type(node.op) in _BINOP:
		return bool(_BINOP[type(node.op)](_eval_boolean(node.left, env), _eval_boolean(node.right, env)))
	raise ValueError("invalid boolean expression node")


@value_equality
class BooleanHamiltonianGate(Gate):
	"""The diagonal gate ``sum_x exp(i theta/2 sum_k f_k(x)) |x><x|``."""

	def __init__(self, parameter_names: Sequence[str], boolean_strs: Sequence[str], theta: float) -> None:
		self._parameter_names = tuple(parameter_names)
		self._boolean_strs = tuple(boolean_strs)
		self._theta = float(theta)
		self._expressions = tuple(_parse_boolean(expr, self._parameter_names) for expr in self._boolean_strs)

	def _num_qubits_(self) -> int:
		return len(self._parameter_names)

	def _unitary_(self) -> Matrix:
		n = len(self._parameter_names)
		dim = 1 << n
		matrix = [[0j] * dim for _ in range(dim)]
		for x in range(dim):
			env = {name: bool((x >> (n - 1 - i)) & 1) for i, name in enumerate(self._parameter_names)}
			total = sum(1 for expression in self._expressions if _eval_boolean(expression, env))
			matrix[x][x] = complex(math.cos(self._theta / 2 * total), math.sin(self._theta / 2 * total))
		return matrix

	def _has_unitary_(self) -> bool:
		return True

	def _value_equality_values_(self) -> Any:
		return (self._parameter_names, self._boolean_strs, self._theta)

	def __str__(self) -> str:
		return f"BooleanHamiltonianGate({self._parameter_names}, {self._boolean_strs}, {self._theta})"

	def __repr__(self) -> str:
		return (
			f"vqubit.BooleanHamiltonianGate(parameter_names={self._parameter_names!r}, "
			f"boolean_strs={self._boolean_strs!r}, theta={self._theta!r})"
		)


__all__ = ["ArithmeticGate", "BooleanHamiltonianGate"]
