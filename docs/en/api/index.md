# API reference

This section is generated automatically from the *docstrings* in the code.

## The `vqubit` package

The public API is flat: everything is reached as `vqubit.Something`. The
implementation is organised into internal modules, but you only need to remember
the package name.

::: vqubit
    options:
      members: false

## The native module `vqubit_native`

The heavy lifting (linear algebra, gate matrices and simulation) lives in the
native extension `vqubit_native`. You do not need to import it directly: the
`vqubit` package wraps it.

::: vqubit_native
