"""The ``vqubit.linalg`` namespace, mirroring ``cirq.linalg``."""

from __future__ import annotations

from ._linalg_public import *  # noqa: F403
from ._linalg_public import __all__ as _public_all
from ._qis import partial_trace

__all__ = [*_public_all, "partial_trace"]
