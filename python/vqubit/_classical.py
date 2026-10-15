"""Classical control: conditions, classically-controlled operations and
classical variable assignment.

This is a minimal classical-data model (enough for :func:`If`,
:class:`ClassicallyControlledOperation` and :class:`SetVariable`); it will be
completed together with the rest of ``cirq.value``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from ._ops import Operation
from ._qubits import Qid
from ._value import value_equality


class Condition:
	"""A condition on classical data (e.g. a measurement result)."""

	def __invert__(self) -> Condition:
		return _NegatedCondition(self)

	def __bool__(self) -> bool:
		raise NotImplementedError


@value_equality
class KeyCondition(Condition):
	"""A condition on the bit at ``index`` of the measurement with ``key``."""

	def __init__(self, key: Any, index: int = 0) -> None:
		self.key = str(key)
		self.index = int(index)

	def _value_equality_values_(self) -> Any:
		return (self.key, self.index)

	def __repr__(self) -> str:
		return f"vqubit.KeyCondition({self.key!r}, index={self.index})"


@value_equality
class _NegatedCondition(Condition):
	"""The negation of another condition."""

	def __init__(self, condition: Condition) -> None:
		self.condition = condition

	def __invert__(self) -> Condition:
		return self.condition

	def _value_equality_values_(self) -> Any:
		return (self.condition,)

	def __repr__(self) -> str:
		return f"~{self.condition!r}"


@value_equality
class SetVariable(Operation):
	"""Stores a classical ``value`` under the measurement key ``key``."""

	def __init__(self, key: Any, value: Any) -> None:
		self.key = str(key)
		self.value = value

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return ()

	def with_qubits(self, *new_qubits: Qid) -> SetVariable:
		return self

	def _value_equality_values_(self) -> Any:
		return (self.key, self.value)

	def __str__(self) -> str:
		return f"SetVariable({self.key}, {self.value})"

	def __repr__(self) -> str:
		return f"vqubit.SetVariable({self.key!r}, {self.value!r})"


@value_equality
class ClassicallyControlledOperation(Operation):
	"""Applies ``sub_operation`` only when its conditions hold."""

	def __init__(self, conditions: Sequence[Condition], sub_operation: Operation) -> None:
		self._conditions = tuple(conditions)
		self._sub_operation = sub_operation

	@property
	def conditions(self) -> tuple[Condition, ...]:
		return self._conditions

	@property
	def sub_operation(self) -> Operation:
		return self._sub_operation

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._sub_operation.qubits

	@property
	def gate(self) -> Any:
		raise TypeError("a classically controlled operation has no unitary")

	def with_qubits(self, *new_qubits: Qid) -> ClassicallyControlledOperation:
		return ClassicallyControlledOperation(self._conditions, self._sub_operation.with_qubits(*new_qubits))

	def _value_equality_values_(self) -> Any:
		return (self._conditions, self._sub_operation)

	def __str__(self) -> str:
		return f"if({', '.join(map(str, self._conditions))}): {self._sub_operation}"

	def __repr__(self) -> str:
		return f"vqubit.ClassicallyControlledOperation({self._conditions!r}, {self._sub_operation!r})"


def If(condition: Condition, *then_ops: Operation, else_ops: Sequence[Operation] = ()) -> Any:
	"""Return classically controlled operations for an if/else on ``condition``."""
	result: list[Operation] = [ClassicallyControlledOperation([condition], op) for op in then_ops]
	result += [ClassicallyControlledOperation([~condition], op) for op in else_ops]
	if len(result) == 1:
		return result[0]
	return tuple(result)


__all__ = [
	"ClassicallyControlledOperation",
	"Condition",
	"If",
	"KeyCondition",
	"SetVariable",
]
