"""vqubit: a dependency-free quantum computing framework for Python.

vqubit lets you build quantum circuits and run them on local simulators, with
no third-party dependencies at all. The numerically heavy parts -- linear
algebra, gate matrices and simulation -- are written in V and compiled into the
native extension ``vqubit_native``; this package is the public, dependency-free
Python API on top of it.

The public API mirrors Cirq closely::

    import vqubit as vq

    q0, q1 = vq.LineQubit.range(2)
    circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1))
    print(circuit)
"""

from __future__ import annotations

from vqubit_native import version as _native_version

__version__ = _native_version()

__all__ = ["__version__"]
