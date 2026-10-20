"""Text diagrams for circuits.

The renderer draws one line per qubit, one column per moment, with ``│``
connecting the qubits of a multi-qubit operation. It is deliberately simple
and stable: the same circuit always produces the same text.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
	from ._circuits import Circuit

_COLUMN_GAP = "───"
_CONNECTOR_GAP = "   "
_CONNECTOR = "│"


def to_text_diagram(circuit: Circuit) -> str:
	"""Return an ASCII diagram of ``circuit``."""
	qubits = circuit.all_qubits()
	if not qubits:
		return ""

	index = {q: i for i, q in enumerate(qubits)}
	widths: list[int] = []
	cells: list[dict[int, str]] = []
	connectors: list[set[int]] = []

	for moment in circuit.moments:
		moment_cells: dict[int, str] = {}
		moment_connectors: set[int] = set()
		width = 1
		for operation in moment.operations:
			symbols = operation._circuit_diagram_info_()
			indices = [index[q] for q in operation.qubits]
			for position, qubit_index in enumerate(indices):
				symbol = symbols[position] if position < len(symbols) else "?"
				moment_cells[qubit_index] = symbol
				width = max(width, len(symbol))
			if indices:
				for gap in range(min(indices), max(indices)):
					moment_connectors.add(gap)
		widths.append(width)
		cells.append(moment_cells)
		connectors.append(moment_connectors)

	label_width = max(len(str(q)) for q in qubits)
	lines: list[str] = []
	for qubit_index, qubit in enumerate(qubits):
		row = f"{qubit!s:>{label_width}}: "
		for width, moment_cells in zip(widths, cells, strict=True):
			row += _COLUMN_GAP
			cell: str | None = moment_cells.get(qubit_index)
			# The cell is left aligned and the rest of the cell is wire, so a bare
			# ``@`` in a five-wide cell is ``@────`` and not ``  @  ``: that is how
			# Cirq draws ``@`` next to ``@^0.5``.
			row += cell.ljust(width, "─") if cell is not None else "─" * width
		row += _COLUMN_GAP
		lines.append(row.rstrip())

		if qubit_index < len(qubits) - 1:
			connector_row = " " * (label_width + 2)
			for width, moment_connectors in zip(widths, connectors, strict=True):
				connector_row += _CONNECTOR_GAP
				# The connector lines up with the start of the cell, as in Cirq.
				connector_row += _CONNECTOR.ljust(width) if qubit_index in moment_connectors else " " * width
			connector_row += _CONNECTOR_GAP
			stripped = connector_row.rstrip()
			if stripped:
				lines.append(stripped)

	return "\n".join(lines)
