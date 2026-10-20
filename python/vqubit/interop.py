"""The ``vqubit.interop`` namespace, mirroring ``cirq.interop``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.interop`` so a
deep import (``from cirq.interop import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	circuit_from_qasm,
	from_cirq,
	qasm,
	quirk_json_to_circuit,
	quirk_url,
	quirk_url_to_circuit,
	to_cirq,
)

__all__ = [
	"circuit_from_qasm",
	"from_cirq",
	"qasm",
	"quirk_json_to_circuit",
	"quirk_url",
	"quirk_url_to_circuit",
	"to_cirq",
]
