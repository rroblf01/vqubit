# Stability and versions

vqubit follows [SemVer 2.0](https://semver.org/) and this page says exactly what
that means here, what counts as public API, and how something is removed without
breaking your code.

!!! info "Since 1.0.0"
    **1.0.0** froze the public surface: from here removing a name or changing its
    behaviour incompatibly needs a **major** version, and anything that goes away
    first spends a deprecation cycle with its `DeprecationWarning`.

## What counts as public API

- Every name in `vqubit.__all__` (`vq.X`, `vq.Circuit`, `vq.simulate`,
  `vq.circuit_from_qasm`…) and the **public members** of its classes
  (`Circuit.moments`, `Gate.num_qubits`…).
- The submodules that mirror Cirq's: `vqubit.ops`, `vqubit.sim`, `vqubit.linalg`,
  `vqubit.protocols`, `vqubit.transformers`, `vqubit.work`.
- The formats: Cirq's JSON, OpenQASM 2.0 text and Quirk links.

Anything starting with `_` (`vqubit._gates`, `_linalg`…) is **internal** and may
change at any time without notice.

## What each number means

| Change | Version |
| ------ | ------- |
| Removing a public name, or changing its behaviour incompatibly | **major** |
| Adding a public name, or a new compatible behaviour | **minor** |
| Fixing a bug without changing the surface | **patch** |

The surface is checkable: `tools/api_snapshot.txt` holds the signature of every
public symbol, and CI fails when one changes without being declared
(`python tools/check_api_surface.py`). A fix that changes a wrong numerical
result is a **patch**, even though the result differs: the contract is what the
reference documents, not the wrong value.

## How something is removed

1. **It is announced** in the `CHANGELOG` (`Deprecated` section) and in the
   reference.
2. **It warns at runtime** with `DeprecationWarning`, saying since which version
   it is deprecated, what to use instead and when it will go away:

   ```python
   import vqubit as vq

   @vq.deprecated(
       "vqubit.old_helper",
       since="0.9",
       replacement="vqubit.new_helper",
       removed_in="1.0",
   )
   def old_helper(value):
       return new_helper(value)
   ```

   For something that is not a function or a class (a renamed parameter, say),
   call `vq.warn_deprecated(...)` where you detect it.
3. **It stays for at least one minor cycle** (or a patch cycle when the minor
   number does not move). Nothing is removed in the same cycle it is announced.
4. **It is removed** only in a major version, and the `CHANGELOG` says so.

Warnings are raised with the right `stacklevel`, so they point at your line and
not at the library. The tests pin them with `pytest.warns(DeprecationWarning)`,
so a change of message is noticed.

## Supported platforms

| | |
| --- | --- |
| Python | CPython **3.11 – 3.15** (an `abi3` wheel for 3.11) |
| Systems | Linux, macOS and Windows, x86_64 and arm64 |
| Dependencies | **none** at runtime: only the compiled extension |
| Optional | `numpy`, `sympy`, `matplotlib` and `cirq`, installed separately: vqubit declares no extras |

Every release is tested on the three systems and five Python versions
(`test.yml`), and the numbers are compared with Cirq in a separate job. A
platform that is not in that matrix is not supported, even if it happens to
work.

## How the project checks it

- `benchmarks/verify.py`: 650 numerical checks against closed forms.
- `benchmarks/parity.py`: 1873 contract checks against cirq-core 1.7.0.
- `tools/check_api_surface.py`: the public surface cannot change undeclared.
- `tools/check_api_docs.py`: every public symbol has a docstring.
- `tools/check_wheel.py`: the published wheel holds exactly what its `RECORD` says.
- `tools/check_reproducible.py`: two builds of the same commit give the same bytes.
- `CHANGELOG.md`: what changed and in which version.

**Versioned** documentation per release (`mike`) is prepared but waiting on one
repository setting (publishing from the `gh-pages` branch); the site publishes
from `main` today.
