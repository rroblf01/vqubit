"""Compiling circuits to a target gate set.

A target gate set says which gates a device can run. The optimizer here rewrites
a circuit so that every operation is either one of those gates or a single-qubit
:class:`~vqubit.PhasedXZGate`, expanding composite gates first and merging runs
of single-qubit gates afterwards.

General two-qubit synthesis (turning an arbitrary two-qubit unitary into the
target's entangling gate) needs the KAK decomposition and is not implemented
yet; a two-qubit gate that is not already in the set is reported rather than
silently miscompiled.
"""

from __future__ import annotations

from ._circuits import Circuit
from ._gates import CNOTGate, CZGate, Gate, H, is_plain_gate
from ._ops import Operation
from ._optimizers import single_qubit_matrix_to_phxz
from ._transformers import expand_composite, merge_single_qubit_moments_to_phxz

DEFAULT_ATOL = 1e-8


class TargetGateset:
	"""The base gate set: anything but single-qubit compilation."""

	def __init__(self, *, atol: float = DEFAULT_ATOL) -> None:
		self.atol = atol

	def two_qubit_gates(self) -> tuple[type[Gate], ...]:
		"""The entangling gates the target can run."""
		return ()

	def convert_two_qubit(self, operation: Operation) -> list[Operation] | None:
		"""Rewrite an entangling gate the target cannot run, or return ``None``."""
		return None

	def synthesize_two_qubit(self, operation: Operation) -> list[Operation] | None:
		"""Synthesize an arbitrary two-qubit gate, or return ``None`` if unsupported."""
		return None

	def compile_operation(self, operation: Operation) -> list[Operation]:
		"""Return ``operation`` rewritten in the target gate set."""
		gate = operation.gate
		if gate._is_measurement_():
			return [operation]
		if len(operation.qubits) == 1:
			matrix = gate._unitary_()
			if matrix is None:
				raise TypeError(f"{gate!r} is not a unitary single-qubit gate")
			phxz = single_qubit_matrix_to_phxz([list(row) for row in matrix], atol=self.atol)
			return [phxz.on(operation.qubits[0])]
		if isinstance(gate, self.two_qubit_gates()) and is_plain_gate(gate):
			return [operation]
		converted = self.convert_two_qubit(operation)
		if converted is not None:
			return converted
		synthesized = self.synthesize_two_qubit(operation)
		if synthesized is not None:
			return synthesized
		raise TypeError(f"{gate!r} is not in the {type(self).__name__} target gate set")

	def optimize_circuit(self, circuit: Circuit) -> Circuit:
		"""Return ``circuit`` compiled to this gate set."""
		expanded = expand_composite(circuit)
		result = Circuit()
		for operation in expanded.all_operations():
			for compiled in self.compile_operation(operation):
				result.append(compiled)
		return merge_single_qubit_moments_to_phxz(result, atol=self.atol)


class CZTargetGateset(TargetGateset):
	"""A gate set whose only entangling gate is ``CZ``."""

	def two_qubit_gates(self) -> tuple[type[Gate], ...]:
		return (CZGate,)

	def convert_two_qubit(self, operation: Operation) -> list[Operation] | None:
		# CNOT = (I x H) CZ (I x H), the H on the target.
		if isinstance(operation.gate, CNOTGate) and is_plain_gate(operation.gate):
			control, target = operation.qubits
			return [H.on(target), CZGate().on(control, target), H.on(target)]
		return None

	def synthesize_two_qubit(self, operation: Operation) -> list[Operation] | None:
		from ._kak import two_qubit_matrix_to_cz_operations

		matrix = operation.gate._unitary_()
		if matrix is None:
			return None
		return two_qubit_matrix_to_cz_operations(
			operation.qubits[0],
			operation.qubits[1],
			[list(row) for row in matrix],
			atol=self.atol,
		)


class CNOTTargetGateset(TargetGateset):
	"""A gate set whose only entangling gate is ``CNOT``."""

	def two_qubit_gates(self) -> tuple[type[Gate], ...]:
		return (CNOTGate,)

	def convert_two_qubit(self, operation: Operation) -> list[Operation] | None:
		# CZ = (I x H) CNOT (I x H), the H on the target.
		if isinstance(operation.gate, CZGate) and is_plain_gate(operation.gate):
			control, target = operation.qubits
			return [H.on(target), CNOTGate().on(control, target), H.on(target)]
		return None

	def synthesize_two_qubit(self, operation: Operation) -> list[Operation] | None:
		from ._kak import two_qubit_matrix_to_operations

		matrix = operation.gate._unitary_()
		if matrix is None:
			return None
		return two_qubit_matrix_to_operations(
			[list(row) for row in matrix],
			qubits=(operation.qubits[0], operation.qubits[1]),
			atol=self.atol,
		)


def optimize_for_target_gateset(circuit: Circuit, gateset: TargetGateset | None = None) -> Circuit:
	"""Compile ``circuit`` to ``gateset`` (a ``CZ`` target by default)."""
	return (gateset or CZTargetGateset()).optimize_circuit(circuit)


__all__ = [
	"CNOTTargetGateset",
	"CZTargetGateset",
	"TargetGateset",
	"optimize_for_target_gateset",
]
