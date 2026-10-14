"""vqubit: A Python extension written in V.

The native code lives in the extension module `vqubit_native`; this package
re-exports it so that `import vqubit` is the only name to remember.
"""

from vqubit_native import *  # noqa: F401,F403
