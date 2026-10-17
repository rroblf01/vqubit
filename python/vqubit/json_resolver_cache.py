"""The ``vqubit.json_resolver_cache`` namespace, mirroring ``cirq.json_resolver_cache``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.json_resolver_cache`` so a
deep import (``from cirq.json_resolver_cache import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	annotations,
)

__all__ = [
	"annotations",
]
