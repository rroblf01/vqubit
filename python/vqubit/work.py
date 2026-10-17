"""The ``vqubit.work`` namespace, mirroring ``cirq.work``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.work`` so a
deep import (``from cirq.work import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	CircuitSampleJob,
	Collector,
	PauliSumCollector,
	Sampler,
	ZerosSampler,
)

__all__ = [
	"CircuitSampleJob",
	"Collector",
	"PauliSumCollector",
	"Sampler",
	"ZerosSampler",
]
