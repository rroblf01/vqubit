"""Extra ``cirq.circuits`` classes: alignment, point optimizers, QASM output
and the text-diagram drawer.
"""

from __future__ import annotations

import dataclasses
import enum
from collections.abc import Callable, Iterable, Sequence

from ._circuits import AbstractCircuit, Circuit, Moment
from ._ops import Operation
from ._qubits import Qid


class Alignment(enum.Enum):
	"""How a new moment is aligned relative to the existing circuit."""

	EARLIEST = "earliest"
	NEW_THEN_INLINE = "new_then_inline"
	INLINE = "inline"
	APPEND = "append"
	LEFT = "left"
	RIGHT = "right"


class CircuitFunction(Circuit):
	"""A circuit whose moments are produced lazily by a function."""

	def __init__(self, function: Callable[[Sequence[Qid]], Iterable[Moment]], qubits: Sequence[Qid]) -> None:
		self._function = function
		self._qubits = tuple(qubits)
		super().__init__(function(self._qubits))

	@property
	def qubits(self) -> tuple[Qid, ...]:
		return self._qubits


@dataclasses.dataclass
class PointOptimizationSummary:
	"""A description of a local change to a circuit."""

	clear_span: int
	clear_qubits: Sequence[Qid]
	new_operations: Sequence[Operation]


class PointOptimizer:
	"""Optimizes a circuit one moment at a time."""

	def __call__(self, circuit: AbstractCircuit) -> Circuit:
		return self.optimize_circuit(circuit)

	def optimize_circuit(self, circuit: AbstractCircuit) -> Circuit:
		result = Circuit()
		for moment in circuit.moments:
			result.append(list(moment.operations))
		return result

	def optimization_at(self, circuit: AbstractCircuit, index: int, op: Operation) -> PointOptimizationSummary | None:
		return None

	def convert(self, op: Operation) -> Sequence[Operation]:
		return [op]


class QasmOutput:
	"""The OpenQASM 2.0 text of a circuit."""

	def __init__(self, circuit: AbstractCircuit) -> None:
		self.circuit = circuit

	def qasm(self) -> str:
		from ._qasm import qasm

		return qasm(self.circuit)

	def __str__(self) -> str:
		return self.qasm()


class TextDiagramDrawer:
	"""A simple text diagram drawer: a grid of cells with wires drawn across it."""

	def __init__(self) -> None:
		self._cells: dict[tuple[int, int], str] = {}
		self._horizontal: set[tuple[int, int]] = set()
		self._vertical: set[tuple[int, int]] = set()
		self.width = 0
		self.height = 0

	def write(self, x: int, y: int, text: str) -> None:
		self._cells[(x, y)] = text
		self.width = max(self.width, x + 1)
		self.height = max(self.height, y + 1)

	def horizontal_line(self, x1: int, x2: int, y: int) -> None:
		for x in range(min(x1, x2), max(x1, x2)):
			self._horizontal.add((x, y))

	def vertical_line(self, y1: int, y2: int, x: int) -> None:
		for y in range(min(y1, y2), max(y1, y2)):
			self._vertical.add((x, y))

	def render(self, *, horizontal_spacing: int = 1, vertical_spacing: int = 0) -> str:
		lines: list[str] = []
		for y in range(self.height):
			cells: list[str] = []
			for x in range(self.width):
				cells.append(self._cells.get((x, y), "-" if (x, y) in self._horizontal else " "))
			lines.append(" " * horizontal_spacing + (" " * horizontal_spacing).join(cells))
			if vertical_spacing and y < self.height - 1:
				lines.extend([""] * vertical_spacing)
		return "\n".join(lines)


__all__ = [
	"AbstractCircuit",
	"Alignment",
	"CircuitFunction",
	"PointOptimizationSummary",
	"PointOptimizer",
	"QasmOutput",
	"TextDiagramDrawer",
]
