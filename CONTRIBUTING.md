# Contributing to vqubit

Thanks for your interest in vqubit. This document explains how to set up a
development environment, the conventions used in the project, and how to run
the checks before opening a pull request.

## Setting up

vqubit is a Python extension written in [V](https://vlang.io) and built with
[vcraft](https://github.com/rroblf01/vcraft). You need:

- Python 3.11 or newer;
- the V compiler and a C compiler (gcc or clang);
- `vcraft`;

```console
$ git clone https://github.com/rroblf01/vqubit
$ cd vqubit
$ python -m venv .venv && source .venv/bin/activate
$ pip install vcraft
$ pip install -r requirements-dev.txt
$ vcraft develop
```

`vcraft develop` builds the native extension and installs an editable pointer
into the virtualenv, so edits to the V sources or the Python package are picked
up on the next rebuild.

## Project layout

```
src/            V sources (numeric kernels and the native bindings)
python/vqubit/  dependency-free Python package (the public API)
tests/          test suite
docs/           bilingual documentation (Spanish and English)
examples/       runnable, dependency-free examples
```

The public API is flat: everything is reached as `vqubit.Something`, even
though the implementation is organised into private modules under
`python/vqubit/`.

## Conventions

- **No runtime dependencies.** The published package must import with nothing
  but the standard library. Numeric work belongs in V.
- **Type hints everywhere**, and the code is checked with `ty`.
- **Lint and format with `ruff`** (line length 120).
- **Tests for everything that makes sense**, written with `pytest`.
- **Bilingual docs.** Every user-facing feature needs a Spanish page and an
  English page under `docs/es` and `docs/en`.
- **Commits** are written in English, one logical change each, with no
  co-author trailer.

## Checks

Run these before opening a pull request; CI runs the same ones:

```console
$ vcraft develop
$ pytest
$ ruff check .
$ ruff format --check .
$ ty check
```

## Reporting issues

Please open an issue at https://github.com/rroblf01/vqubit/issues. If the
problem is in the build tool rather than vqubit itself, report it at
https://github.com/rroblf01/vcraft/issues.
