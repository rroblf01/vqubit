"""An operation applied only when its control qubits have the given values."""

from __future__ import annotations

from collections.abc import Sequence
from itertools import product
from typing import Any

from ._control_values import AbstractControlValues
from ._gates import ControlledGate, Gate
from ._ops import Operation
from ._qubits import Qid
from ._value import value_equality


@value_equality
class ControlledOperation(Operation):
	"""An operation with additional control qubits."""

	def __init__(
		self,
		controls: Sequence[Qid],
		sub_operation: Operation,
		control_values: Sequence[bool | int] | AbstractControlValues | None = None,
	) -> None:
		self._controls = tuple(controls)
		self._sub_operation = sub_operation
		if control_values is None:
			self._control_values: Sequence[bool | int] | AbstractControlValues = (True,) * len(self._controls)
		elif isinstance(control_values, AbstractControlValues):
			self._control_values = control_values
		else:
			self._control_values = tuple(bool(v) for v in control_values)

	@property
	def controls(self) -> tuple[Qid, ...]:
		return self._controls

	@property
	def sub_operation(self) -> Operation:
		return self._sub_operation

	@property
	def control_values(self) -> Sequence[bool | int] | AbstractControlValues:
		return self._control_values

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._controls + self._sub_operation.qubits

	def _bool_control_values(self) -> tuple[bool, ...]:
		values = self._control_values
		if isinstance(values, AbstractControlValues):
			for combo in product([False, True], repeat=len(self._controls)):
				if combo in values:
					return combo
			raise ValueError("no satisfying control assignment")
		return tuple(bool(v) for v in values)

	@property
	def gate(self) -> Gate:
		return ControlledGate(
			self._sub_operation.gate,
			num_controls=len(self._controls),
			control_values=self._bool_control_values(),
		)

	def _unitary_(self) -> Any:
		return self.gate._unitary_()

	def with_qubits(self, *new_qubits: Qid) -> ControlledOperation:
		n = len(self._controls)
		return ControlledOperation(
			new_qubits[:n], self._sub_operation.with_qubits(*new_qubits[n:]), self._control_values
		)

	def _value_equality_values_(self) -> Any:
		return (self._controls, self._sub_operation, self._control_values)

	def __str__(self) -> str:
		return f"ControlledOperation({self._controls}, {self._sub_operation})"

	def __repr__(self) -> str:
		return f"vqubit.ControlledOperation({self._controls!r}, {self._sub_operation!r}, control_values={self._control_values!r})"


__all__ = ["ControlledOperation"]
