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
- [uv](https://docs.astral.sh/uv/) for the development tooling.

```console
$ git clone https://github.com/rroblf01/vqubit
$ cd vqubit
$ python -m venv .venv && source .venv/bin/activate
$ pip install vcraft
$ uv pip install --group dev
$ vcraft develop
```

The dev and docs tooling are `[dependency-groups]` in `pyproject.toml`. Use
`uv pip install --group docs` for the documentation toolchain. `uv pip install`
adds to the virtualenv without pruning it, unlike `uv sync`, which would remove
the editable extension `vcraft develop` installs.

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
- **Tests for everything that makes sense**, written with `pytest`. Property
  tests use Hypothesis (a dev dependency; skip when it is missing) and the
  parser fuzzers use a seeded `random` so they run everywhere.
- **Every public symbol has a docstring.** `tools/check_api_docs.py` keeps a
  baseline of the ones that do not yet; a new symbol not in the baseline fails
  the check. Prune the baseline with `--write` when you document something.
- **Bilingual docs.** Every user-facing feature needs a Spanish page and an
  English page under `docs/es` and `docs/en`.
- **Commits** are written in English, one logical change each, with no
  co-author trailer.

## Checks

Run these before opening a pull request; CI runs the same ones:

```console
$ vcraft develop
$ pytest                        # the whole suite
$ pytest --cov --cov-fail-under=95   # CI also enforces 95 % coverage
$ ruff check .
$ ruff format --check .
$ ty check
$ python benchmarks/verify.py   # 619 numerical checks
$ python benchmarks/parity.py   # 1867 API checks against cirq-core 1.7.0
$ python tools/check_api_docs.py # new public API must be documented
$ mkdocs build --strict         # both languages
```

The test suite must also pass **without** numpy, sympy or matplotlib: the
matrix installs only the wheel and pytest, so import those lazily or guard
them with `pytest.importorskip`. 

## Reporting issues

Please open an issue at https://github.com/rroblf01/vqubit/issues. If the
problem is in the build tool rather than vqubit itself, report it at
https://github.com/rroblf01/vcraft/issues.
