## What does this change?

<!-- One paragraph. If it fixes an issue, write "Fixes #123". -->

## Checklist

- [ ] `pytest -q` is green (and `pytest` still passes with numpy, sympy and
      matplotlib absent — the test matrix installs only the wheel and pytest).
- [ ] `ruff check .` and `ruff format --check .` are clean.
- [ ] `ty check` is clean.
- [ ] `python benchmarks/verify.py` (619 numerical checks) is green.
- [ ] `python benchmarks/parity.py` (1867 API checks against cirq-core 1.7.0) is
      green; any new tolerated divergence is in `benchmarks/parity_tolerated.toml`
      with a `reason` and a `plan`.
- [ ] `mkdocs build --strict` is green.
- [ ] New public API has a docstring (`python tools/check_api_docs.py`), and the
      baseline in `tools/undocumented.txt` did not grow.
- [ ] If a behaviour changed, `CHANGELOG.md` has an entry and the docs say so.

## Notes for the reviewer

<!-- Anything surprising: a decision, a divergence, a vcraft limitation. -->
