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
- **Every change gets a changelog fragment**: one small file in `changelog.d/`,
  named `<issue-or-pr>.<type>.md` (`+name.<type>.md` when there is no issue).
  `towncrier build --version X` folds them into `CHANGELOG.md`; see
  `changelog.d/README.md` for the types. A deprecation also needs the
  `vqubit.deprecated` decorator and a test, per
  [the stability policy](docs/en/guide/stability.md).

## Checks

Run these before opening a pull request; CI runs the same ones:

```console
$ vcraft develop
$ pytest                        # the whole suite
$ pytest --cov --cov-fail-under=95   # CI also enforces 95 % coverage
$ ruff check .
$ ruff format --check .
$ ty check
$ python benchmarks/verify.py   # 650 numerical checks
$ python benchmarks/parity.py   # 1873 API checks against cirq-core 1.7.0
$ python benchmarks/perf.py     # four workloads against the ceilings
$ python tools/check_reproducible.py  # two builds, same bytes
$ python tools/check_api_docs.py # new public API must be documented
$ mkdocs build --strict         # both languages
```

The build jobs also run `python tools/check_wheel.py dist/*.whl`, which reads
the freshly built wheel back and verifies every `RECORD` hash and size; run it
after a local build if you touched packaging (`pip` does not check `RECORD`).

The test suite must also pass **without** numpy, sympy or matplotlib: the
matrix installs only the wheel and pytest, so import those lazily or guard
them with `pytest.importorskip`. 

## Releasing a version

CI checks everything before the tag; the release itself is a tag, a GitHub
release and one click on PyPI's side.

1. **Bump the version** in `vcraft.toml` and the string returned by
   `src/vqubit_native.v` (`version()`); `tests/test_deprecation_and_release.py`
   fails if the two disagree with `vqubit.__version__`.
2. **Close the changelog**: rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD`
   and update `CITATION.cff`.
3. **Run the gate** from the checks above, then let CI finish green on `main`.
4. **Tag and push**: `git tag vX.Y.Z && git push origin vX.Y.Z`.
5. **Publish the GitHub release** for that tag. `release.yml` builds the wheels
   and the sdist, verifies each `RECORD`, attests the artifacts and uploads to
   PyPI with Trusted Publishing.
6. **Verify what PyPI serves**: `python -m venv /tmp/check &&
   /tmp/check/bin/python -m pip install vqubit==X.Y.Z`, then import it and run
   `python tools/check_wheel.py dist/*.whl` on what `pip download` fetched.
7. **Deploy the versioned docs** (below) and start the next `## [Unreleased]`
   section.

A release candidate is optional: the `sdist` job already builds a wheel from the
published sdist and installs it in a clean virtualenv, and the matrix installs
the wheel on three systems and five Pythons, so `1.0.0rc1` only adds a manual
look at PyPI's metadata. The roadmap keeps it as a fallback.

## Publishing the documentation

The site publishes from `main` with `actions/deploy-pages`. To keep a copy per
release as well (**versioned docs**, the `mike` toolchain is already in the
`docs` group):

```console
$ uv pip install --group docs
$ mike deploy --push --update-aliases 1.0.0 latest
$ mike set-default --push latest
```

`mike` writes to the `gh-pages` branch, so GitHub Pages has to be switched from
"GitHub Actions" to "Deploy from a branch: `gh-pages`" in the repository
settings; that is the one step a maintainer does in the browser. Until then the
site keeps publishing `main` as it does today.

## Reporting issues

Please open an issue at https://github.com/rroblf01/vqubit/issues. If the
problem is in the build tool rather than vqubit itself, report it at
https://github.com/rroblf01/vcraft/issues.
