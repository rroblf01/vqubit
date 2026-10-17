"""The ``vqubit.circuits`` namespace, mirroring ``cirq.circuits``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.circuits`` so a
deep import (``from cirq.circuits import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	AbstractCircuit,
	Alignment,
	Circuit,
	CircuitOperation,
	FrozenCircuit,
	InsertStrategy,
	Moment,
	PointOptimizationSummary,
	PointOptimizer,
	QasmOutput,
	TextDiagramDrawer,
)

__all__ = [
	"AbstractCircuit",
	"Alignment",
	"Circuit",
	"CircuitOperation",
	"FrozenCircuit",
	"InsertStrategy",
	"Moment",
	"PointOptimizationSummary",
	"PointOptimizer",
	"QasmOutput",
	"TextDiagramDrawer",
]
