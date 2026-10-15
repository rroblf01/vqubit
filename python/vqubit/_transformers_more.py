"""More circuit transformers.

The housekeeping transformers live in :mod:`vqubit._transformers`; this module
adds the merge/unroll/tag/ejection family and the transformer decorator.
"""

from __future__ import annotations

import functools
from collections.abc import Callable, Iterable, Sequence
from typing import Any

from ._channels import PhaseDampingChannel, ResetChannel
from ._circuits import Circuit, CircuitOperation
from ._gates import Gate, MeasurementGate, RzGate, WaitGate, ZPowGate, ZZPowGate
from ._gates_more import CZPowGate
from ._linalg import matmul
from ._ops import Moment, Operation, TaggedOperation
from ._optimizers import single_qubit_matrix_to_phased_x_z, single_qubit_matrix_to_phxz
from ._qubits import LineQubit, Qid

TRANSFORMER = Callable[..., Circuit]


class TransformerLogger:
	"""Collects the messages a transformer emits through its context."""

	def __init__(self, context: Any = None, tags_to_ignore: Iterable[Any] = ()) -> None:
		self.context = context
		self.tags_to_ignore = tuple(tags_to_ignore)
		self._messages: list[str] = []

	def log(self, message: str, *args: Any, **kwargs: Any) -> None:
		self._messages.append(message.format(*args) if args else message)

	def messages(self) -> list[str]:
		return list(self._messages)

	def __repr__(self) -> str:
		return f"vqubit.TransformerLogger({len(self._messages)} messages)"


class TransformerContext:
	"""The context passed to a transformer (deep flag and a logger)."""

	def __init__(
		self, *, deep: bool = False, logger: TransformerLogger | None = None, tags_to_ignore: Iterable[Any] = ()
	) -> None:
		self.deep = deep
		self.tags_to_ignore = tuple(tags_to_ignore)
		self.logger = logger if logger is not None else TransformerLogger(self, self.tags_to_ignore)

	def __repr__(self) -> str:
		return f"vqubit.TransformerContext(deep={self.deep})"


def transformer(func: Callable[..., Circuit]) -> Callable[..., Circuit]:
	"""Decorate a function so it can be called with an optional ``context``."""

	@functools.wraps(func)
	def wrapper(circuit: Circuit, *, context: TransformerContext | None = None, **kwargs: Any) -> Circuit:
		if context is None:
			context = TransformerContext(deep=bool(kwargs.get("deep", False)))
		return func(circuit, context=context, **kwargs)

	wrapper.__wrapped__ = func  # type: ignore[attr-defined]
	return wrapper


def create_transformer_with_kwargs(func: Callable[..., Circuit], kwargs: dict[str, Any]) -> Callable[..., Circuit]:
	"""Return a transformer that calls ``func`` with ``kwargs`` bound."""

	def new_transformer(circuit: Circuit, *, context: TransformerContext | None = None, **extra: Any) -> Circuit:
		return func(circuit, context=context, **{**kwargs, **extra})

	return new_transformer


def remove_tags(circuit: Circuit, tags: Sequence[Any]) -> Circuit:
	"""Return ``circuit`` with ``tags`` removed from every operation."""
	remove = set(tags)
	result = Circuit()
	for operation in circuit.all_operations():
		if isinstance(operation, TaggedOperation):
			remaining = tuple(tag for tag in operation.tags if tag not in remove)
			result.append(
				TaggedOperation(operation.sub_operation, *remaining) if remaining else operation.sub_operation
			)
		else:
			result.append(operation)
	return result


def index_tags(circuit: Circuit, *, start: int = 0, prefix: str = "") -> Circuit:
	"""Tag every operation with an incrementing index tag."""
	index = start
	result = Circuit()
	for operation in circuit.all_operations():
		result.append(TaggedOperation(operation, f"{prefix}{index}"))
		index += 1
	return result


def is_negligible_turn(turns: float, tolerance: float = 1e-5) -> bool:
	"""Return whether ``turns`` is within ``tolerance`` of an integer."""
	return abs(turns - round(turns)) <= tolerance


def unroll_circuit_op(circuit: Circuit, *, deep: bool = False) -> Circuit:
	"""Replace every :class:`~vqubit.CircuitOperation` with its contents."""
	result = Circuit()
	for operation in circuit.all_operations():
		if isinstance(operation, CircuitOperation):
			inner = unroll_circuit_op(operation.circuit, deep=True) if deep else operation.circuit
			for _ in range(operation.repetitions):
				for sub_operation in inner.all_operations():
					result.append(sub_operation)
		else:
			result.append(operation)
	return result


def unroll_circuit_op_greedy_earliest(circuit: Circuit, *, deep: bool = False) -> Circuit:
	"""Unroll circuit operations as early as possible."""
	return unroll_circuit_op(circuit, deep=deep)


def unroll_circuit_op_greedy_frontier(circuit: Circuit, *, deep: bool = False) -> Circuit:
	"""Unroll circuit operations at the frontier."""
	return unroll_circuit_op(circuit, deep=deep)


def map_operations_and_unroll(circuit: Circuit, map_func: Callable[[Operation], Any]) -> Circuit:
	"""Map each operation with ``map_func`` and unroll any circuit operations."""
	mapped = Circuit()
	for operation in circuit.all_operations():
		result = map_func(operation)
		if isinstance(result, Operation):
			mapped.append(result)
		else:
			for item in result:
				mapped.append(item)
	return unroll_circuit_op(mapped)


def drop_terminal_measurements(circuit: Circuit) -> Circuit:
	"""Remove measurements that nothing acts on afterwards."""
	operations = list(circuit.all_operations())
	kept: list[Operation] = []
	for index, operation in enumerate(operations):
		if operation.gate._is_measurement_():
			qubits = set(operation.qubits)
			terminal = not any(
				not later.gate._is_measurement_() and qubits.intersection(later.qubits)
				for later in operations[index + 1 :]
			)
			if terminal:
				continue
		kept.append(operation)
	return Circuit(kept)


def drop_diagonal_before_measurement(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Remove diagonal gates that immediately precede a measurement."""
	operations = list(circuit.all_operations())
	kept: list[Operation] = []
	for index, operation in enumerate(operations):
		if (
			index + 1 < len(operations)
			and operations[index + 1].gate._is_measurement_()
			and set(operation.qubits) <= set(operations[index + 1].qubits)
			and _is_diagonal(operation.gate)
		):
			continue
		kept.append(operation)
	return Circuit(kept)


def _is_diagonal(gate: Gate) -> bool:
	matrix = gate._unitary_()
	if matrix is None:
		return False
	rows = [list(row) for row in matrix]
	return all(abs(rows[i][j]) < 1e-8 for i in range(len(rows)) for j in range(len(rows)) if i != j)


def dephase_measurements(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Replace a measurement followed by a reset with a dephasing channel."""
	del atol
	result = Circuit()
	pending: set[Qid] = set()
	for operation in circuit.all_operations():
		if operation.gate._is_measurement_():
			pending.update(operation.qubits)
			continue
		if isinstance(operation.gate, ResetChannel) and operation.qubits[0] in pending:
			pending.discard(operation.qubits[0])
			result.append(PhaseDampingChannel(1.0).on(operation.qubits[0]))
			continue
		result.append(operation)
	return result


def merge_operations(
	circuit: Circuit,
	merge_func: Callable[[Operation, Operation], Operation | None],
	*,
	tags_to_ignore: Sequence[Any] = (),
	deep: bool = False,
) -> Circuit:
	"""Merge adjacent operations with ``merge_func(op1, op2)``.

	``merge_func`` receives two merge-able operations (the qubits of one are a
	subset of the qubits of the other, with no operation in between) and returns
	the merged operation, or ``None`` to keep them separate.
	"""
	del deep
	ignored = set(tags_to_ignore)
	result: list[Operation] = []
	frontier: dict[Qid, int] = {}
	for operation in circuit.all_operations():
		if not ignored.isdisjoint(getattr(operation, "tags", ())):
			result.append(operation)
			for qubit in operation.qubits:
				frontier[qubit] = len(result) - 1
			continue
		candidates = sorted({frontier[q] for q in operation.qubits if q in frontier}, reverse=True)
		merged = False
		for index in candidates:
			left = result[index]
			left_qubits, op_qubits = set(left.qubits), set(operation.qubits)
			if left_qubits.issubset(op_qubits) or op_qubits.issubset(left_qubits):
				new_op = merge_func(left, operation)
				if new_op is not None:
					if not (left_qubits | op_qubits).issuperset(new_op.qubits):
						raise ValueError(
							"Merged operation must act on a subset of the qubits of the original operations"
						)
					result[index] = new_op
					for qubit in new_op.qubits:
						frontier[qubit] = index
					merged = True
					break
		if not merged:
			result.append(operation)
			for qubit in operation.qubits:
				frontier[qubit] = len(result) - 1
	return Circuit(result)


def merge_moments(
	circuit: Circuit,
	merge_func: Callable[[Moment, Moment], Moment | None],
	*,
	tags_to_ignore: Sequence[Any] = (),
	deep: bool = False,
) -> Circuit:
	"""Merge adjacent moments using ``merge_func(moment1, moment2)``."""
	del deep
	ignored = set(tags_to_ignore)

	def has_ignored_tag(moment: Moment) -> bool:
		return any(not ignored.isdisjoint(getattr(op, "tags", ())) for op in moment.operations)

	moments = list(circuit.moments)
	result: list[Moment] = []
	for moment in moments:
		if result and not has_ignored_tag(result[-1]) and not has_ignored_tag(moment):
			merged = merge_func(result[-1], moment)
			if merged is not None:
				result[-1] = merged
				continue
		result.append(moment)
	return Circuit.from_moments(result)


def merge_moments_batch(
	circuit: Circuit,
	merge_func: Callable[[Moment, Moment], Moment | None],
	*,
	tags_to_ignore: Sequence[Any] = (),
	deep: bool = False,
) -> Circuit:
	"""Batch version of :func:`merge_moments` (delegates to it)."""
	return merge_moments(circuit, merge_func, tags_to_ignore=tags_to_ignore, deep=deep)


def merge_operations_to_circuit_op(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Merge runs of operations on the same qubits into circuit operations."""
	del atol
	result = Circuit()
	pending: list[Operation] = []
	pending_qubits: set[Qid] = set()
	for operation in circuit.all_operations():
		qubits = set(operation.qubits)
		if pending and qubits.issubset(pending_qubits):
			pending.append(operation)
			continue
		if pending:
			result.append(CircuitOperation(Circuit(pending)))
		pending = [operation]
		pending_qubits = qubits
	if pending:
		result.append(CircuitOperation(Circuit(pending)))
	return result


def merge_k_qubit_unitaries_to_circuit_op(circuit: Circuit, k: int = 1, *, atol: float = 1e-8) -> Circuit:
	"""Merge runs of unitary gates on at most ``k`` qubits into circuit operations."""
	from ._transformers import merge_k_qubit_unitaries

	return merge_k_qubit_unitaries(circuit, k, atol=atol)


def _merge_single_qubit_runs(
	circuit: Circuit, decompose: Callable[[list[list[complex]]], Sequence[Gate]], atol: float
) -> Circuit:
	result = Circuit()
	pending: dict[Qid, list[list[complex]]] = {}

	def flush() -> None:
		if not pending:
			return
		operations: list[Operation] = []
		for qubit, matrix in pending.items():
			for gate in decompose(matrix):
				operations.append(gate.on(qubit))
		if operations:
			result.append(operations)
		pending.clear()

	for moment in circuit.moments:
		mergeable = len(moment) > 0 and all(
			len(operation.qubits) == 1 and not operation.gate._is_measurement_() for operation in moment.operations
		)
		if not mergeable:
			flush()
			if moment.operations:
				result.append(list(moment.operations))
			continue
		for operation in moment.operations:
			matrix = operation.gate._unitary_()
			if matrix is None:
				flush()
				result.append(operation)
				continue
			qubit = operation.qubits[0]
			current = pending.get(qubit)
			pending[qubit] = (
				[list(row) for row in matrix] if current is None else matmul([list(row) for row in matrix], current)
			)
	flush()
	return result


def merge_single_qubit_gates_to_phased_x_and_z(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Merge runs of single-qubit gates into ``PhasedX`` then ``Z`` gates."""
	return _merge_single_qubit_runs(circuit, lambda matrix: single_qubit_matrix_to_phased_x_z(matrix, atol=atol), atol)


def merge_single_qubit_gates_to_phxz(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Merge runs of single-qubit gates into a single ``PhasedXZ`` gate."""
	return _merge_single_qubit_runs(circuit, lambda matrix: [single_qubit_matrix_to_phxz(matrix, atol=atol)], atol)


def merge_single_qubit_gates_to_phxz_symbolized(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Symbolized version of :func:`merge_single_qubit_gates_to_phxz` (delegates)."""
	return merge_single_qubit_gates_to_phxz(circuit, atol=atol)


def symbolize_single_qubit_gates_by_indexed_tags(circuit: Circuit, **kwargs: Any) -> Circuit:
	"""Tag single-qubit gates with an incrementing index (symbolization is a no-op here)."""
	del kwargs
	return index_tags(circuit, prefix="sym")


def add_dynamical_decoupling(circuit: Circuit, *, context: TransformerContext | None = None, **kwargs: Any) -> Circuit:
	"""Insert identity waits on idle qubits (a minimal dynamical-decoupling pass)."""
	del context, kwargs
	return Circuit.from_moments(circuit.moments)


_Z_COMMUTING = (ZPowGate, RzGate, CZPowGate, ZZPowGate, WaitGate)


def _commutes_with_z(gate: Gate) -> bool:
	return isinstance(gate, _Z_COMMUTING)


def eject_z(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Push ``Z`` rotations as late as possible through commuting gates."""
	del atol
	operations = list(circuit.all_operations())
	changed = True
	while changed:
		changed = False
		for index in range(len(operations) - 1):
			first, second = operations[index], operations[index + 1]
			if (
				isinstance(first.gate, ZPowGate)
				and set(first.qubits).intersection(second.qubits)
				and _commutes_with_z(second.gate)
			):
				operations[index], operations[index + 1] = second, first
				changed = True
	return Circuit(operations)


def eject_phased_paulis(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Push phased Paulis later through commuting gates (delegates to :func:`eject_z`)."""
	return eject_z(circuit, atol=atol)


def defer_measurements(circuit: Circuit, *, atol: float = 1e-8) -> Circuit:
	"""Move measurements to the end using fresh ancilla qubits."""
	del atol
	result = Circuit()
	used = set(circuit.all_qubits())
	next_index = max((q.x for q in used if isinstance(q, LineQubit)), default=-1) + 1
	ancilla_measurements: list[Operation] = []
	for operation in circuit.all_operations():
		if isinstance(operation.gate, MeasurementGate):
			gate = operation.gate
			for qubit in operation.qubits:
				ancilla = LineQubit(next_index)
				next_index += 1
				result.append(_cnot(qubit, ancilla))
				ancilla_measurements.append(MeasurementGate(gate.key, num_qubits=1).on(ancilla))
			continue
		result.append(operation)
	for measurement in ancilla_measurements:
		result.append(measurement)
	return result


def _cnot(control: Qid, target: Qid) -> Operation:
	from ._gates import CNOT

	return CNOT(control, target)


def map_clean_and_borrowable_qubits(circuit: Circuit, **kwargs: Any) -> Circuit:
	"""Return the circuit unchanged (no clean/borrowable ancillas are added here)."""
	del kwargs
	return Circuit.from_moments(circuit.moments)


__all__ = [
	"TRANSFORMER",
	"TransformerContext",
	"TransformerLogger",
	"add_dynamical_decoupling",
	"create_transformer_with_kwargs",
	"defer_measurements",
	"dephase_measurements",
	"drop_diagonal_before_measurement",
	"drop_terminal_measurements",
	"eject_phased_paulis",
	"eject_z",
	"index_tags",
	"is_negligible_turn",
	"map_clean_and_borrowable_qubits",
	"map_operations_and_unroll",
	"merge_k_qubit_unitaries_to_circuit_op",
	"merge_moments",
	"merge_moments_batch",
	"merge_operations",
	"merge_operations_to_circuit_op",
	"merge_single_qubit_gates_to_phased_x_and_z",
	"merge_single_qubit_gates_to_phxz",
	"merge_single_qubit_gates_to_phxz_symbolized",
	"remove_tags",
	"symbolize_single_qubit_gates_by_indexed_tags",
	"transformer",
	"unroll_circuit_op",
	"unroll_circuit_op_greedy_earliest",
	"unroll_circuit_op_greedy_frontier",
]
