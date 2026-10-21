# Upstream: issues to report to `vcraft`

Things vqubit needs from [`rroblf01/vcraft`](https://github.com/rroblf01/vcraft).
Each section below is a ready-to-paste issue: title, what happens, why it matters
for vqubit, a reproduction, and the expected behaviour.

Environment for every report: **vcraft 1.4.5**, Linux `x86_64`, Python 3.14/3.15,
project `vqubit` (zero runtime dependencies, single `abi3` wheel per platform).

| # | Title | Blocks |
| - | ----- | ------ |
| 1 | `.pyi` emits enum members as annotations, so the stub is unusable | `ty`, `mypy`, `stubtest` (roadmap 1.4, 7.3) |
| 2 | `.pyi` writes `__init__` for classes built in `tp_new` | `stubtest` (roadmap 7.3) |
| 3 | No way to declare optional dependencies / extras | `vqubit[vis]`, `vqubit[interop]` (roadmap 3.5) |
| 4 | Huge pages for the buffers (`buffer_new_*`) | 20-qubit simulation (roadmap 8) |
| 5 | `WHEEL: Generator:` carries the project version; `audit` rejects `.pyd` | packaging hygiene |
| 6 | Runtime dispatch (function multiversioning) for the hot kernels | SIMD from a portable wheel (roadmap 8) |
| 7 | `NativeTableau` has no `copy()`; `CliffordOp` is missing gates; `NativeDensity` has no `set_density_matrix` | per-step Clifford snapshots, stim's gate set, density-matrix trial states (roadmap 0.3) |

---

## 1. The generated `.pyi` declares enum members as annotations

> **Workaround in the project (0.9):** `pyproject.toml` tells mypy to
> `ignore_errors` and `ignore_missing_imports` for the `vqubit_native` module,
> so the enum diagnostic (and the case where the generated stub is not in the
> checkout) does not mask the rest of the check; `ty` reads the same stub and
> passes. The fix still belongs upstream (`member = value`), because the stub is
> also what users see when they inspect the extension.

**What happens.** For a V enum returned across the FFI boundary, the generated
stub writes one annotated attribute per member:

```python
class CliffordOp(enum.IntEnum):
    h: CliffordOp
    s: CliffordOp
    ...
```

Under the typing spec an enum member in a stub must be an assignment (`h = 0` or
`member = value`), not an annotation. Modern `mypy` rejects the file outright:

```
error: Detected enum "vqubit_native.CliffordOp" in a type stub with zero members.
  There is a chance this is due to a recent change in the semantics of enum
  membership. If so, use `member = value` to mark an enum member, instead of
  `member: type`  [misc]
```

**Why it matters.** Any project that generates a stub for a V enum cannot type
check its own code, and `mypy.stubtest` cannot even start (see issue 2). For
vqubit it blocks the "audit `__all__` ↔ docstrings ↔ annotations" and
"`stubtest` over `vqubit_native.pyi`" workstreams.

**Reproduction.** Build any vcraft project that returns a V enum and run mypy on
the generated stub, or `python -m mypy.stubtest <module>`.

**Expected.** `h: Final = ...`/`h = ...` — the annotation form is not a member.

---

## 2. The generated `.pyi` writes `__init__` for a class built in `tp_new`

**What happens.** vcraft generates

```python
class NativeState:
    num_qubits: int
    def __init__(self, num_qubits: int) -> None: ...
```

but at runtime the constructor lives in `tp_new` (the Python class has no
parameters in `__init__`) — which is the right implementation. The stub is what
is wrong. `stubtest` (after working around issue 1) reports, for each such class:

```
error: vqubit_native.NativeState.__init__ is inconsistent, runtime does not have
  parameter "num_qubits". You may need to write stubs for __new__ instead of __init__.
error: vqubit_native.NativeState is a disjoint base at runtime, but isn't marked
  with @disjoint_base in the stub
```

(3 classes → 6 errors.)

**Why it matters.** The stub is shipped in the wheel and is the only description
of the native API a type checker or IDE sees. As written, `NativeState(4)` is
reported as an error by strict checkers even though it is the documented use.

**Expected.** `def __new__(cls, num_qubits: int) -> NativeState: ...`, and, if the
runtime truly is a disjoint base, the corresponding marker (or a targeted
`# type: ignore`), so `stubtest` passes cleanly.

---

## 3. No way to declare optional dependencies (`extras`)

> **Update (1.0):** vqubit 1.0 ships **no extras by design**, like Cirq, so
> nothing is blocked any more: the optional pieces are documented as plain
> `pip install matplotlib` / `cirq-core` / `sympy` lines and detected at run
> time. The request below stands for any project that does want extras — a
> `METADATA` patched in the pipeline would make the published wheel differ from
> a local build, which is worse than not having them.

**What happens.** `vcraft.toml` has `[package] dependencies`, rendered as
`Requires-Dist`. There is no way to declare an *optional* dependency, and
`extras` in the wheel builder means "extra files to ship", not Python extras. So
no `Provides-Extra:` / `Requires-Dist: numpy; extra == "vis"` can be emitted.

**Why it matters.** vqubit's core promise is *zero runtime dependencies*, with
heavyweights (`numpy`, `sympy`, `matplotlib`, `cirq-core` for the interoperability
bridge) offered as opt-in extras. Without extras we cannot publish
`pip install vqubit[vis]` / `vqubit[interop]`, which is a documented 1.0.0
deliverable.

**Suggested shape.** Accept PEP 621's `[project.optional-dependencies]` (or a
`[package.optional-dependencies]` table) and render one `Provides-Extra:` per
key plus the conditional `Requires-Dist` entries; mirror the same in `sdist`.

---

## 4. Huge pages for the buffers (`buffer_new_*` / `buffer_new_f64`)

**What happens.** Buffers are allocated with a plain `malloc`-style allocation.
For a 20-qubit state vector (2²⁰ `complex128` = **16 MiB**) every gate application
streams 16 MiB in and 16 MiB out. With 4 KiB pages that is ~8192 distinct page
translations per gate, which exceeds the L2 TLB (~1.5–2k entries): every kernel
pass thrashes the TLB, and the kernels end up at ~25 GB/s — the machine's memory
bandwidth — so `-march` (more SIMD) does not help.

**Why it matters.** This is the remaining 20-qubit gap versus Cirq. It is a
memory-access problem, not arithmetic, so it cannot be fixed from vqubit.

**Request.** Allocate the big buffers with an explicit huge-page hint behind an
opt-in/auto threshold:

- Linux: `mmap(MAP_PRIVATE | MAP_ANONYMOUS)` + `madvise(MADV_HUGEPAGE)`, with an
  optional `MAP_HUGETLB` path when `/proc/meminfo` advertises free huge pages;
- Windows: `VirtualAlloc(..., MEM_LARGE_PAGES)` when the privilege is available,
  plain `VirtualAlloc` otherwise;
- macOS: documented as best-effort (`vm_allocate` + superpage advice).

A `VCRAFT_HUGEPAGES=1|0` toggle (auto by default above some size) would let us
measure and roll back, and would keep the behaviour explicit rather than magic.
Ideally the same threshold applies to the state, density and tableau buffers.

---

## 5. `WHEEL` metadata `Generator:` and `audit` on Windows

**5a. `Generator:` reports the project's version, not vcraft's.**
`vlib/vcraft_wheel/build.v:144` passes `'vcraft ${input.version}'`, where
`input.version` is the *project* version. The built wheel for vqubit 0.2.0
therefore contains:

```
Wheel-Version: 1.0
Generator: vcraft 0.2.0        # should be vcraft 1.4.5 (PEP 427: the tool's version)
```

It should be the version of vcraft that generated the wheel.

**5b. `vcraft audit` calls a Windows extension "not a shared object".**
`vlib/vcraft_project/env.v:318` rejects the name unless it contains `.so`:

```v
extension := result.extension
if !extension.contains('.so') {
    problems << '${extension} is not a shared object'
}
```

On Windows the extension is `.pyd`, so a perfectly valid `win_amd64` wheel is
reported as broken. Accepting `.so`/`.pyd`/`.dylib` (or `.pyd` on Windows names)
would fix it.

**5c. Reproducibility is an accident of the ZIP writer, and undocumented.** A
wheel built from the same commit is byte-identical to another, which is what lets
a third party verify the hash published on PyPI: `vlib/vcraft_wheel/zip.v:100`
writes `modification time 0`, `modification date 0x21` (1980-01-01) for every
entry and `build.v` adds the files in a fixed order, and neither `METADATA` nor
`WHEEL` carries a date. Nothing in the documentation promises it, and there is no
`SOURCE_DATE_EPOCH` handling, so when the date is needed for some other reason
there is no knob. Expected: document the fixed timestamp as a guarantee and, if
it ever becomes configurable, derive it from `SOURCE_DATE_EPOCH`
(`max(epoch, 1980-01-01)`). vqubit checks the guarantee on every push with
`tools/check_reproducible.py`.

---

## 6. Runtime dispatch (function multiversioning) for the hot kernels

**Request, lower priority.** Ship SIMD-optimised variants of the hot kernels and
select one at load time (CPUID on x86-64, `HWCAP` on ARM), so a *portable* wheel
gets SIMD without baking `-march` into it. Today `-march` is a documented opt-in
for a local build, because baking it in would break portability and diverge from
the tested artifact.

This is a "nice to have" for vqubit 1.0.0: issue 4 is much more valuable per unit
of effort, because the 20-qubit workload is bandwidth-bound rather than
compute-bound.

---

## 7. `NativeTableau` has no `copy()`, `CliffordOp` is missing gates, `NativeDensity` has no `set_density_matrix`

> **Update (0.9):** the missing `copy()` is no longer a blocker.
> `CliffordSimulator.simulate_moment_steps` rebuilds the tableau per moment by
> replaying the operations seen so far (O(moments²) rather than O(moments)); a
> `copy()` would make that linear and is still the right fix upstream.

**Request, two related gaps in the stabilizer kernels.**

**7a. No `copy()` on `NativeTableau`.** `vqubit.CliffordSimulator` wants a
`simulate_moment_steps` that yields one snapshot per moment, the way
`cirq.CliffordSimulator` does. Every other simulator can snapshot its state (a
list of amplitudes, a density matrix), but a tableau cannot be copied: the type
exposes only `apply`, `measure`, `reset`, `state_vector` and
`stabilizer_tableau` (`vlib/vcraft_project` `NativeTableau`). Without a copy, the
snapshots of a step-by-step simulation all alias the final state, which is worse
than not offering the method.

```python
import vqubit_native as native
tableau = native.NativeTableau(2)
print([name for name in dir(tableau) if not name.startswith("_")])
# ['apply', 'measure', 'num_qubits', 'reset', 'stabilizer_tableau', 'state_vector']
```

Expected: `NativeTableau.copy()` (or a `__copy__`), so
`CliffordSimulator.simulate_moment_steps` can hand out independent snapshots.

**7c. No `set_density_matrix` on `NativeDensity`.** A density-matrix trial result
should be able to hand back a simulation state holding it (as
`cirq.DensityMatrixTrialResult.get_state_containing_qubit` does), but the native
type only accepts amplitudes (`set_from_amplitudes`) or ``|0...0>``
(`set_zero_state`), so a *mixed* state cannot be loaded into it. Expected:
`NativeDensity.set_density_matrix(flat)` (the interleaved form
`buffer.density_matrix()` already returns) so states round-trip in and out.

**7b. `CliffordOp` is missing instructions that stim has.**

> **Update (1.0):** closed on vqubit's side without new instructions. The six
> gates are now rewritten as short words over the nine native ones
> (`python/vqubit/_clifford.py`, `_decomposed_word`), so
> `has_stabilizer_effect` and `CliffordSimulator` answer like Cirq's,
> `benchmarks/verify.py` compares them gate by gate (31 checks in the
> `stabilizer` section) and `tests/test_clifford_cliffords.py` compares the
> states with both simulators. Native opcodes would still be faster — a word is
> up to nine instructions where stim applies one — so the request below stands
> as a performance item, not a correctness one.

vqubit's Clifford gate
set is `h, s, sdg, x, y, z, cnot, cz, swap, measure`. `cirq` (through stim) also
treats `X ** 0.5`, `Y ** 0.5`, `iSWAP`, `XX`, `YY` and `ZZ` as stabilizer gates,
and vqubit reports `has_stabilizer_effect(...) is False` for them because the
kernels cannot apply them:

```python
import cirq, vqubit as vq
q0, q1 = cirq.LineQubit.range(2)
print(cirq.has_stabilizer_effect(cirq.XX(q0, q1)), vq.has_stabilizer_effect(vq.XX(q0, q1)))
# True False
```

(The ones vqubit *can* express are already accepted: `X ** 3`, `CNOT ** -1`,
`Z ** 1.5`; only the instructions are missing.)

Expected: `sx`, `sy`, `iswap`, `xx`, `yy`, `zz` (and their inverses, where stim
has them) in `CliffordOp`, applied by `NativeTableau.apply`. With those, vqubit's
`has_stabilizer_effect` and `CliffordSimulator` would answer like Cirq's for the
whole stim gate set, and `benchmarks/verify.py` could compare them gate by gate.

**7d. `NativeTableau` cannot apply an arbitrary Clifford.** It applies one of the
ten instruction codes, so a `vqubit.CliffordGate` — a gate kept as a unitary,
the way `cirq.CliffordGate` is — cannot be run on the stabilizer simulator, even
though `has_stabilizer_effect` reports `True` for it. Cirq's simulator accepts it
because its tableau is updated from the gate's own tableau. Workaround in place:
those gates are simulated with `vqubit.Simulator`, and the limitation is
documented in `guide/cirq-differences` (it is a divergence, not silence).

```python
import vqubit as vq
q0, q1 = vq.LineQubit.range(2)
print(vq.has_stabilizer_effect(vq.CXSWAP))  # True
vq.CliffordSimulator().simulate(vq.Circuit(vq.CXSWAP(q0, q1)))
# TypeError: CliffordSimulator doesn't support ...
```

Expected: a way to apply a Clifford given its tableau or its action on the Pauli
generators, e.g. `NativeTableau.apply_tableau(xs, zs, rs)` for the gate's
`2n`-by-`2n` symplectic block, or `CliffordOp.from_unitary`. That would let the
stabilizer simulator take any `CliffordGate` (and let `CXSWAP`/`CZSWAP` — vqubit's
own extras — run there).
