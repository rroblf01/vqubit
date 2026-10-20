# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

An audit of the numerical core (driven by the 0.5 coverage work) found several
helpers that returned the wrong answer rather than an error. They are fixed
here, each with a regression test, and the test suite grew from 690 to 1391
tests with coverage at 95 %.

### Fixed

- **`unitary_eig` and `map_eigenvalues` on non-Hermitian normal matrices.**
  The shifted-QR solver used real shifts, so a rotation (`[[0,-1],[1,0]]`)
  returned eigenvalues `[0, 0]` and a permutation matrix stalled forever. The
  eigenvalues now come from the characteristic polynomial (Durand-Kerner) and
  the eigenvectors from the null spaces of `A - λI`, which also keeps the basis
  unitary when an eigenvalue repeats.
- **`sub_state_vector`** returned the factor of the wrong subsystem and failed
  on valid product states such as `|10>`. It now builds the reduced density
  matrix, keeps the eigenvector with eigenvalue one and randomises the global
  phase, exactly as `cirq.sub_state_vector` does, and validates the input.
- **`density_matrix_kronecker_product`** was computing `vec(t1) vec(t2)^T`
  instead of the Kronecker product (wrong shape, zero trace).
- **`targeted_conjugate_about`** raised `IndexError` whenever NumPy was
  installed, because the shape and conjugation helpers only understood lists.
- **`kron`** rejected a scalar mixed with matrices (`kron(2, X)`).
- **`is_diagonal`** accepted rectangular matrices; `is_unitary`, `is_hermitian`
  and `match_global_phase` indexed empty input; `so4_to_magic_su2s` did not
  check that its input is 4x4.
- **`parameter_names` / `parameter_symbols`** missed parameters inside
  expressions (`X ** (2*t + 1)`) and inside circuits; they now walk the
  expression's names and the value's state.
- **`ControlledGate`** reported itself as not parameterized and returned itself
  unchanged when asked to resolve its sub-gate's parameters.
- **`PauliString`** is iterable over its qubits (iterating used to loop
  forever), has `len`, `values` and `get`, and `ps ** -1` returns the string
  instead of the identity.
- **`PauliSum`** has value equality and a hash, as `cirq.PauliSum` does.
- **`Collector.collect`** now stops dispensing jobs once `max_total_samples` is
  spent, instead of draining the whole batch.
- **Mixed kinds of qid sort deterministically**: `sorted([NamedQubit('a'),
  LineQubit(0), GridQubit(0, 0)])` used to raise `TypeError`.
- **`qid_shape`** honours the declared dimension of a gate, so a qutrit channel
  reports `(3,)` rather than `(2,)`.
- **`unitary` of a decomposition that leaves qubits idle** (a uniform
  superposition over a subset, say) now spans the gate's full width.
- **`CliffordTableau.from_flat`** validates the length of the flat layout.

- **`ControlledGate` accepts a control-value formula.** A `SumOfProducts` or
  `ProductOfSums` condition used to collapse to its first satisfying
  assignment, so a gate guarded by `q0 XOR q1` fired on one branch only; every
  satisfying assignment now fires it, and `ControlledOperation` passes the
  formula through instead of reducing it.
- **`CliffordTableauSimulationState` accepts identity gates.** `I` and `X ** 2`
  were rejected as "not a Clifford gate" because they share the `None`
  instruction code with unsupported gates.
- **Impossible measurement outcomes are rejected** with `ValueError` instead of
  silently returning an all-zero state: the native collapse kernel returns a
  zero vector rather than an empty one, so the probability is now checked
  first.
- **Qudit shapes are rejected explicitly** by `measure_state_vector`,
  `sample_state_vector` and the density-matrix variants, which used to treat a
  qutrit as a qubit and truncate it.
- **`SingleQubitCliffordGate.from_pauli` works.** It built a three-entry map and
  handed it to `from_double_map`, which requires two, so every call raised
  `ValueError`. `**` on a Cliffords gate no longer recurses forever on negative
  powers, and the power is taken modulo the element's order (twelve) instead of
  four, which was wrong for the eight elements of order three.
- **QASM import covers more of qelib1**: `cswap`, `crz`, `cu1` and `cu3` work
  (the last three asked for one qubit instead of two), and an indexed classical
  assignment such as `c[0] = measure q[0];` is dispatched to the measurement
  handler instead of being parsed as a gate.
- **Quirk import reads a controlled `Y`, `Z` or `H`**, which `quirk_url` has
  always exported: the control-target lookup only knew `X`.
- **Parameter expressions print like Cirq's**: `X ** t` is drawn `X^t` instead of
  `X^(1.0 * t)`, because the identity elements of `*`, `+` and `/` are no longer
  left in the expression tree.

- **The Quirk and JSON importers only raise `ValueError`.** A malformed
  document used to leak a `KeyError` (`'cols'`, `'qubits'`, a missing field of a
  named type), an `AttributeError` (an entry that is not a string) or an
  `IndexError` (a URL without `circuit=`, an empty `cols`). They now validate
  the document's shape and wrap the rest in `ValueError`, which is what the
  fuzzers in `tests/test_fuzz_parsers.py` pin down.
- **`PauliString`'s constructor accepts a mapping keyed by any qid type.**
  `dict[LineQubit, Pauli]` did not match `Mapping[Qid, …]` under an invariant
  key type, so a literal raised a type error in strict checkers.

- **`ConstantQubitNoiseModel` no longer drops its noise.** When a moment already
  used a qubit it caught the overlapping-operation error and returned the
  moment untouched, silently un-noising the circuit; it now returns the moment
  followed by a noise moment, and `NoiseModel.noisy_moments` flattens that.
  `NoiseModelFromNoiseProperties` was a no-op and now inserts a depolarizing
  channel per operation with the device's gate error.
- **`Operation.num_qubits()`** exists, as in Cirq (the docs promised it, the
  method did not), so `vq.H(q0).num_qubits() == 1` works while `vq.H.num_qubits`
  stays the attribute.

### Tests (0.5)

- Property tests with Hypothesis for the algebraic laws (powers, rotation sums,
  unitarity, the circuit unitary against the simulator, QASM/JSON round trips,
  Pauli multiplication and expectation values, sweep lengths) and
  deterministic fuzzers for the three parsers.

### Changed

- `sub_state_vector` raises `ValueError` for a state that is not normalized,
  as Cirq does, instead of returning the default.

### Tests

- New suites for the sampling collectors, the experiments, the plotting
  helpers, the stabilizer tableau, parameterized gates, Pauli operations, the
  compilation and transformer extras, the protocol helpers and the numerical
  regressions above. Coverage of the published wheel: 79 % → 95 %, with every
  module above 80 % and the CI threshold raised to match.

## [0.3.0] - 2026-10-19

Behaviour parity with Cirq 1.7.0: the plain gates *are* the power gates, value
equality is periodic, the simulators take Cirq's arguments and answer its calls
(sweeps, moment steps, expectation values, batches), the Pauli observables are
backed by real expectation-value maths, and the four `NotImplementedError` gaps
are closed. Two signatures change: `simulate`/`run` take `param_resolver` second
and `Linspace`/`Points` take Cirq's calling form and endpoint convention.

### Added

- **`**` on the plain gates.** `vq.CZ ** 0.5`, `vq.CNOT ** 0.5`,
  `vq.SWAP ** 0.5`, `vq.ISWAP ** 0.5`, `vq.SQRT_ISWAP ** 0.5`, `vq.CCX ** 0.5`,
  `vq.TOFFOLI ** 0.5` and `vq.I ** 0.5` now return the matching power gate, as
  the Cirq idioms do, and `vq.rz(0.3) ** 0.5` scales the angle to `vq.rz(0.15)`.
  `vq.FREDKIN` and `vq.CSWAP` stay non-powerable, like `cirq.CSWAP`. The plain
  gates themselves are the power gates now: see *Changed*.
- **OpenQASM 2.0 import.** `circuit_from_qasm` is a hand-written parser (no
  `ply`, no dependencies at all) that reads what `qasm` writes -- so the two
  round-trip -- plus the `qelib1` gate set, `pi` expressions, custom `gate`
  definitions (expanded inline), measurements, `reset`, `barrier` and `if(c==N)`
  conditions; a malformed program raises `QasmParseError` with the line number.
  The measurement key is now the classical register verbatim (`creg m[2]`)
  instead of carrying an `m_` prefix, which is what lets the round trip keep the
  key.
- **An optional Cirq bridge.** `to_cirq` and `from_cirq` convert circuits both
  ways through the JSON format the two libraries share; Cirq has to be installed
  (`ImportError` otherwise), since it is not a dependency.
- **Quirk covers more.** Measurements are written as `Measure`, rotations as
  powers (`X^t`, with `t = rads/π`), and arbitrary `X`/`Y`/`Z` powers and their
  controlled versions round-trip; `iSWAP` and friends still raise `TypeError`.
  The JSON also writes rotations as Cirq's `Rx`/`Ry`/`Rz` with `rads` (they used
  to be `RxGate`/`RyGate`/`RzGate`), so Cirq can read them and the bridge keeps
  the class, and it writes `Exponent`… whatever the library writes.
- **The async entry points**, `run_async`, `run_sweep_async` and `run_batch_async`,
  as thin coroutines over the synchronous ones (vqubit's simulation is local, so
  there is nothing to await, but Cirq's call sites work).
- **`get_state_containing_qubit`** on a state-vector trial result and
  `StepResult.sample_measurement_ops`, the two result attributes that were still
  missing.
- **`Linspace` and `Points` take Cirq's calling form** (`Linspace("t", 0, 1, 4)`,
  `Points("t", [0.5, 2.0])`) with Cirq's convention that a `Linspace` **includes**
  both ends (`0, 1/3, 2/3, 1`), where vqubit's own form excluded `stop`; both keep
  vqubit's mapping forms and carry `metadata`, as Cirq's sweeps do.
- **A Clifford tableau decomposes into a circuit.**
  `decompose_clifford_tableau_to_operations` implements the Aaronson–Gottesman
  synthesis (the same one Cirq uses), so a tableau can be turned back into
  `H`/`X`/`Z`/`S`/`CNOT`/`SWAP` operations. `CliffordTableau` also gained
  `kron`, `factor` and `reindex` (and `supports_factor` is `True` now), which
  Cirq's tableau leaves unimplemented; the CH-form simulation state applies
  operations and measures instead of raising.
- **`ResetChannel` takes the qudit's *dimension*,** as Cirq's does
  (`ResetChannel(2)` is a qubit, `ResetChannel(3)` a qutrit), where the argument
  used to mean a number of qubits. Its Kraus operators work for any dimension.
- **Expectation values.** `PauliString.expectation_from_state_vector` and
  `expectation_from_density_matrix`, and the same on `PauliSum`, with Cirq's
  `qubit_map`, `atol` and `check_preconditions`;
  `simulate_expectation_values`, `simulate_expectation_values_sweep` and
  `simulate_expectation_values_sweep_iter` on the state-vector and density-matrix
  simulators; and `sample_expectation_values` on all three, which estimates them
  by measuring each Pauli string in its own basis.
- **A `PauliString` takes a Pauli gate or its letter** (`vq.PauliString({q0: vq.Z})`,
  as Cirq does), and the Pauli arithmetic grew the operators Cirq has: `ps + ps`,
  `s - ps`, `-s`, `s * ps`, `s / 2` and the rest build the `PauliString` or
  `PauliSum` Cirq builds.
- **`bloch_vector_of` and `density_matrix_of`** on the state-vector trial and
  step results, `qubit_map` and `measurements` on every trial and step result,
  and `state` on a Clifford step result, under the names Cirq gives them.
- **The simulator entry points Cirq code reaches for.** `run_sweep`,
  `run_sweep_iter`, `simulate_sweep`, `simulate_sweep_iter` and `run_batch` on
  all three simulators; `simulate_moment_steps` (the state after each moment) on
  the state-vector and density-matrix ones; `Simulator.compute_amplitudes`; and
  `CliffordSimulator.is_supported_operation`. `simulate` and `run` take a
  `param_resolver` now, and the sweeps come from the existing
  `Points`/`Linspace`/`Product`/`Zip` classes. The trial and sample results carry
  `params`, and the trial results also expose Cirq's `qubit_map`, `measurements`
  and `final_state_vector`.
- **A Cirq contract harness**, `benchmarks/parity.py`, run in the `cirq-parity`
  CI job: it checks names, `**`, `isinstance`, gate classes, text diagrams and
  value equality against Cirq, and keeps the accepted divergences in
  `benchmarks/parity_tolerated.toml`, where an entry that stops silencing a real
  divergence fails the build.

### Changed

- **Value equality is periodic, the way Cirq's is.** A power gate now compares
  its eigenphases modulo the gate's period, so `vq.X ** 3 == vq.X`,
  `vq.CNOT ** -1 == vq.CNOT` and `vq.Z ** 0.5 == vq.Z ** 2.5`, exactly as in
  Cirq, while `vq.ISWAP ** -1 != vq.ISWAP` (its second eigenspace sits at half a
  turn) and `vq.rx(0.3) ** 3 != vq.rx(0.3)`, as before.
  `PhasedXPowGate` compares the canonical exponent modulo 2, and
  `vqubit.PeriodicValue` normalizes its value on the way in, so
  `PeriodicValue(3, 2) == PeriodicValue(1, 2)` and both read `1`.
- **The gate families accept a gate *type*,** as Cirq's do:
  `ParallelGateFamily(vqubit.ZPowGate)` now means *any* `ZPowGate` (an instance
  still means that one gate), `AnyIntegerPowerGateFamily(vqubit.CZPowGate)` takes
  the class as well as an instance and no longer raises on a symbolic exponent,
  and `ParallelGateFamily` takes Cirq's `max_parallel_allowed` (its `name` and
  `description` are keyword-only now, as in Cirq).
- **`AnyUnitaryGateFamily` takes `num_qubits` first, as Cirq does.**
  `AnyUnitaryGateFamily(2)` is now the two-qubit unitary gates, where it used to
  set the family's name to the integer `2` and accept everything.
- **`simulate` and `run` take `param_resolver` as their second argument,** where
  `run` used to take `repetitions`. The Cirq idiom
  `simulator.simulate(circuit, resolver)` works now; a positional
  `simulator.run(circuit, 100)` does not, so pass `repetitions=` (as the
  documentation always did). A circuit whose symbols are still unresolved raises
  `ValueError` instead of failing somewhere deeper.
- **`simulate` and `simulate_moment_steps` run measurements**, sampling them,
  collapsing the state and reporting the outcomes in `result.measurements`, as
  Cirq's simulator does; a `ValueError` used to be raised instead.
  `permit_terminal_measurements` in the expectation-value methods therefore does
  what it says.
- **The Clifford simulator takes integer powers of its gates**: `X ** 0` is
  nothing, `X ** 3` is `X` and `CNOT ** -1` is `CNOT`, as Cirq's does, and when it
  refuses an operation it names it (`TypeError`), the way Cirq's message does.
- **`to_resolvers(None)` is a single empty resolver and `to_resolvers([])` is
  none at all**, both as in Cirq, which used to raise `TypeError`.

- **The plain gates *are* the power gates, as in Cirq.** `vq.CZ` is now
  `CZPowGate(exponent=1)`, `vq.CNOT`/`vq.CX` are `CXPowGate`, `vq.SWAP` is
  `SwapPowGate`, `vq.ISWAP` is `ISwapPowGate`, `vq.SQRT_ISWAP` is
  `ISwapPowGate(exponent=0.5)`, and `vq.CCX`/`vq.TOFFOLI` are `CCXPowGate`, so
  `isinstance(vq.CZ, vq.CZPowGate)` and `vq.CZ == vq.CZPowGate()` are `True`
  where they used to be `False`. Two visible consequences: `type(vq.CZ).__name__`
  is `"CZPowGate"`, and `vq.SQRT_ISWAP` is no longer an instance of
  `SqrtISwapGate` (the class survives as an `ISwapPowGate` subclass that
  defaults to exponent 0.5). The old class names — `CZGate`, `CNOTGate`,
  `SWAPGate`, `CCXGate`, `ISWAPGate` — remain as **aliases** of the power
  classes, so code that imports or compares them keeps working; the instruction
  tables (QASM, Clifford, target gatesets, neutral atoms) were audited so that
  only the exponent-1 gate maps to the plain instruction.
- **Text diagrams left-align their symbols and annotate the exponent**, as Cirq
  does: `vq.CZ ** 0.5` draws `@` over `@^0.5` instead of `@` over `Z`.

### Fixed

- **Resolving a parameter *to zero* raised `ZeroDivisionError`.** The helper
  evaluating a symbolic expression built a dictionary of all five operations at
  once, so `t * 12` with `t = 0` also evaluated `t / 12`; every sweep whose first
  point was `0.0` crashed. It dispatches on the operation now.
- **`has_stabilizer_effect` raised `TypeError`** for anything that is not a
  Clifford gate (`has_stabilizer_effect(vq.T(q0))`) instead of returning
  `False`, and said `False` for the identity and for measurements where Cirq says
  `True`.
- **`global_shift` did not do what Cirq documents.** It was applied as a phase
  that depended on the eigenvalue (`exp(iπ·s)` on the `+1` block, `exp(-iπ·s)` on
  the `-1` one) instead of the **global phase** `exp(iπ·s·e)`, so
  `vq.ZPowGate(exponent=0.5, global_shift=0.25)` — and every controlled gate built
  on it — differed from Cirq, and `XX`/`YY`/`ZZ` had no `global_shift` argument at
  all. All 17 power-gate families now agree with Cirq for every exponent and
  shift; the new `global shift` section of `benchmarks/verify.py` compares them
  exactly (510 checks).
- **`vq.CSWAP ** 1` and `** -1` raised**, even though every gate has those two
  powers: `gate ** 1` is the gate itself and `gate ** -1` its inverse, as in
  Cirq's `Gate.__pow__`. Fractional powers still raise, like `cirq.CSWAP`.
- **`is_native_neutral_atom_gate` accepted any `CZPowGate`** — the square root
  included — and rejected `X`. It now answers against Cirq's gateset: integer
  powers of `CNOT`/`CCX`/`CZ`/`CCZ`, the parallel families of
  `X`/`Y`/`Z`/`PhasedX`, and measurement and identity.
- **JSON wrote the plain gates under their old type names.** `to_json(vq.CNOT)`
  produced `{"cirq_type": "CNOTGate"}`, where Cirq writes
  `{"cirq_type": "CXPowGate", "exponent": 1.0, "global_shift": 0.0}`, so the two
  files were not interchangeable. The sixteen power gates now round-trip under
  their Cirq names and parameters, and the legacy names still load.
- **`vq.inverse(vq.ISWAP)` returned iSWAP, not its inverse.** The shared
  `_FixedGate._inverse_` returned `type(self)()`, which is only correct for the
  self-inverse gates, so `U(iSWAP) @ U(iSWAP⁻¹)` was not the identity. It now
  returns `ISWAP_INV`; `vq.inverse(vq.SQRT_ISWAP)` returns `SQRT_ISWAP_INV`.

## [0.2.1] - 2026-10-19

Correctness fixes in the state-vector simulator and the noise channels, plus
documentation. No name changed.

### Fixed

- **Measurements now collapse the state in `Simulator`.** When every measurement
  sat at the end of the circuit, the simulator sampled each `measure` operation
  independently from its own marginal, so separate measurements came out
  decorrelated: the two qubits of a Bell state measured with two `measure`
  calls agreed only half the time, and measuring the same qubit twice could
  report two different bits. The terminal measurements are now sampled jointly
  from the union of the measured qubits, so they agree as they must.
  `DensityMatrixSimulator` and `CliffordSimulator` were already correct.
- **Channels reject probabilities outside `[0, 1]`.** `bit_flip`, `phase_flip`,
  `bit_phase_flip`, `depolarize`, `amplitude_damp`, `phase_damp` and their
  channel classes used to accept, say, `bit_flip(2.0)` and silently return a
  density matrix that did not add up to one. They now raise `ValueError`, as
  Cirq does.

### Documentation

- New **Differences from Cirq** page (English and Spanish), listing the call
  conventions that still differ: renamed arguments, different argument order,
  `**` on operations, the missing `noise=`/`dtype=` in the simulators, and the
  names that mean something else.
- Fixed the platform list on the documentation home page, which omitted Windows
  and claimed an Intel macOS wheel that is not published.
- Fixed the channel table in the noise chapter, which showed argument names the
  constructors do not accept.
- Fixed the Spanish documentation link in the README (`/es/` does not exist; the
  Spanish documentation is served at the site root).

## [0.2.0] - 2026-10-18

Platforms, packaging and Cirq's submodule layout: Windows wheels, a source
distribution, smaller (stripped) wheels and deep imports, built with vcraft
1.4.5. No existing name changed.

### Added

- **Windows support.** `win_amd64` wheels (CPython 3.11 and newer through the
  stable ABI) are now built and published, so `pip install vqubit` works on
  Windows without a compiler. The `Operating System :: Microsoft :: Windows`
  classifier is now accurate.
- **Source distribution.** A `vqubit-0.2.0.tar.gz` sdist is published alongside
  the wheels.
- **Cirq's submodules.** `vqubit.ops`, `vqubit.sim`, `vqubit.work`,
  `vqubit.vis`, `vqubit.protocols`, `vqubit.transformers` and the rest now
  mirror `cirq`'s submodules, so deep imports such as `from cirq.sim import
  Simulator` keep working when `vqubit` is used as a drop-in. `vqubit.linalg`
  and `vqubit.testing` gained the few names they were missing.

### Changed

- **Smaller wheels.** The native extension is stripped, cutting about a fifth
  off the wheel size.
- **Native records and enums.** The internal `jacobi_eig`, `simultaneous_diag`
  and `kron_factor` helpers now return named records, and `NativeTableau.apply`
  takes a `CliffordOp` enum, instead of concatenated buffers and bare integers.
  Both are implementation details of `vqubit_native`.
- Built with **vcraft 1.4.5**.

## [0.1.0] - 2026-10-15

The first release: a dependency-free, Cirq-compatible quantum computing
framework for Python, with the numerical core written in V.

### Added

- **Full `cirq.CliffordTableau` API**: the tableau stores the stabilizer
  generators and implements the Aaronson-Gottesman updates (`apply_x/y/z/h/cx/cz`,
  `stabilizers`, `destabilizers`, `measure`, `then`, `inverse`, `matrix`, ...).
- **`vqubit.work`**: `ZerosSampler`, `CircuitSampleJob`, `Collector` and
  `PauliSumCollector`; **`vqubit.vis`**: `get_state_histogram`,
  `plot_state_histogram`, `integrated_histogram`, `Heatmap`,
  `TwoQubitInteractionHeatmap`, `plot_density_matrix` and
  `scatter_plot_normalized_kak_interaction_coefficients` (Matplotlib/NumPy are
  imported lazily, so the core stays dependency-free); and `with_debug`.
- **Metadata**: authors and maintainers reach the wheel's `METADATA`, and the
  `LICENSE` file travels with the wheel.
- **Qubits**: `Qid`, `Qubit`, `LineQubit`, `GridQubit`, `NamedQubit`.
- **Gates**: `I`, `X`, `Y`, `Z`, `H`, `S`, `T`; the rotation gates `Rx`, `Ry`,
  `Rz`; the power gates `XPowGate`, `YPowGate`, `ZPowGate`, `HPowGate`; the
  phased gates `PhasedXPowGate` and `PhasedXZGate`; the two-qubit gates `CNOT`,
  `CZ`, `SWAP`, `ISWAP` and the Ising gates `XX`, `YY`, `ZZ`; the three-qubit
  `CCX` (Toffoli) and `CSWAP` (Fredkin); `MeasurementGate`; arbitrary unitaries
  through `MatrixGate`; and controlled gates through `Gate.controlled_by`.
- **Circuits**: `Circuit` and `FrozenCircuit` with Cirq's earliest-moment
  packing, `Moment`, `GateOperation`, and ASCII text diagrams.
- **Pauli algebra**: `Pauli`, `PauliString` (multiplication, powers,
  commutation, matrix) and `PauliSum`.
- **Symbolic parameters**: `Parameter`, `ParamResolver`, `is_parameterized` and
  `resolve_parameters`, with symbolic angles on the rotation and power gates.
- **Protocols**: `unitary`, `has_unitary`, `num_qubits`, `qid_shape`, `inverse`,
  `pow`, `decompose`, `commutes`, `definitely_commutes`, `is_measurement`,
  `measurement_key_names`, `pauli_expansion`, `phase_by`, `apply_unitary`,
  `kraus`, `has_kraus`, and `Circuit.unitary` / `Circuit.to_text_diagram`.
- **Simulators**: a state-vector `Simulator`, a density-matrix
  `DensityMatrixSimulator`, and a stabilizer-tableau `CliffordSimulator`
  (Gottesman-Knill), all with sampling and `Result.histogram`.
- **Noise channels**: bit-, phase- and bit-phase-flip, depolarizing, amplitude-
  and phase-damping, and reset, with their Kraus operators.
- **Optimizers and transformers**: `drop_empty_moments`,
  `drop_negligible_operations`, `expand_composite`, `stratified_circuit`,
  `align_left`/`align_right`, `merge_single_qubit_moments_to_phxz`,
  `single_qubit_matrix_to_gates`, `single_qubit_matrix_to_phased_x_z`,
  `single_qubit_matrix_to_phxz`, and the KAK-based
  `two_qubit_matrix_to_operations` / `two_qubit_matrix_to_cz_operations`.
- **Target gate sets**: `TargetGateset`, `CZTargetGateset`,
  `CNOTTargetGateset` and `optimize_for_target_gateset`, which compile
  arbitrary single- and two-qubit gates to a device's gate set.
- **Serialization**: `to_json`/`read_json` (and gzip variants) following Cirq's
  `cirq_type` convention, and OpenQASM 2.0 export with `qasm`.
- **Study**: `sample`, `run`, `run_sweep`, the `Points` and `Linspace` sweeps,
  and a `Sampler`.
- **Devices and local execution**: `Device`, `GridDevice`, `LineDevice`, local
  emulations of Google's `Sycamore` (device and `SycamoreGate`), an ion trap, a
  neutral-atom register and AQT, and an `Engine` that runs circuits locally.
- **Quantum information**: `state_vector`, `density_matrix`, the
  `to_valid_state_vector` / `to_valid_density_matrix` helpers,
  `bloch_vector_from_state_vector`, `partial_trace`, `fidelity`,
  `von_neumann_entropy`, `random_state_vector`, `random_density_matrix` and
  `kraus_channel`.
- **More protocols**: `approx_eq`, `equal_up_to_global_phase`, `MeasurementKey`,
  `measurement_key_obj` / `measurement_key_objs`, `control_keys`,
  `with_measurement_key_mapping` and `trace_distance_bound`.
- **More operations and values**: `TaggedOperation`, `CircuitOperation`,
  `GlobalPhaseGate` / `global_phase`, `WaitGate`, `QubitOrder` and
  `DensePauliString`.
- **Study**: the `Sweep` hierarchy (`SingleSweep`, `Zip`, `Product`, `Concat`),
  `to_sweep`, `to_resolvers`, `TrialResult` and `run_batch`.
- **More transformers**: `map_operations`, `map_moments`, `toggle_tags`,
  `merge_k_qubit_unitaries` and `synchronize_terminal_measurements`.
- **More gates and devices**: the `FSimGate`, `SqrtISwapGate` (`SQRT_ISWAP`) and
  `MSGate` two-qubit gates, the `XmonDevice`, and
  `DeviceMetadata.duration_of`.
- **Testing helpers**: a `vqubit.testing` module with `assert_equivalent_repr`,
  `assert_same_circuits`,
  `assert_circuits_with_terminal_measurements_are_equivalent`,
  `assert_allclose_up_to_global_phase` and `random_circuit`, and eval-able gate
  `repr`s.
- **Interop**: `QasmArgs` and `quirk_url`, which renders a circuit as a Quirk
  link.
- **Linear algebra**: the `cirq.linalg` helpers, also available as
  `vqubit.linalg`: predicates (`is_unitary`, `is_hermitian`, `is_normal`,
  `is_cptp`, `matrix_commutes`, ...), combinators (`dot`, `kron`, `block_diag`,
  `kron_with_controls`), operator spaces (`PAULI_BASIS`, `kron_bases`,
  `pow_pauli_combination`, `expand_matrix_in_orthogonal_basis`),
  transformations (`match_global_phase`, `targeted_left_multiply`,
  `sub_state_vector`, `to_special`), diagonalization, and the KAK family
  (`KakDecomposition`, `kak_vector`, `kak_canonicalize_vector`,
  `num_cnots_required`, `unitary_eig`, `map_eigenvalues`). `kak_decomposition`
  now returns a `KakDecomposition`.
- **Optional numpy/sympy interop**: when numpy is installed the public API
  returns `numpy.ndarray` results (from `unitary`, `kraus`, `PauliString.matrix`,
  `Circuit.unitary`, `state_vector`, `density_matrix` and the simulation
  results) and accepts numpy inputs; the core stays dependency-free otherwise.
  A dedicated CI job runs the whole suite with numpy installed.
- **More gates**: the controlled power gates (`CXPowGate`, `CYPowGate`,
  `CZPowGate`, `CCXPowGate`, `CCYPowGate`, `CCZPowGate`, plus the
  `CNotPowGate`/`CCNotPowGate` aliases), `SwapPowGate`, `ISwapPowGate`,
  `DiagonalGate`/`TwoQubitDiagonalGate`/`ThreeQubitDiagonalGate`,
  `QubitPermutationGate`, `QuantumFourierTransformGate`, `PhasedISwapPowGate`,
  and the instances `CY`, `CCZ`, `CCY`, `CCNOT`, `ISWAP_INV`, `SQRT_ISWAP_INV`.
- **More factories**: `q`, `rx`/`ry`/`rz`, `ms`, `cphase`, `givens`, `wait`,
  `qft`, `reset_each`, `identity_each` and `global_phase_operation`.
- **Qids**: `LineQid`, `GridQid` and `NamedQid`.
- **Op trees and gatesets**: `OP_TREE`, `flatten_to_ops`, `flatten_op_tree`,
  `flatten_to_ops_or_moments`, `freeze_op_tree`, `transform_op_tree`,
  `Gateset` and `GateFamily`.
- **Complete `cirq.ops` surface**: the controlled power gates and swap family,
  `PhasedISwapPowGate`, `PhasedFSimGate`, `PhaseGradientGate`,
  `UniformSuperpositionGate`, `DiagonalGate` and its two/three-qubit variants,
  `QubitPermutationGate`, `QuantumFourierTransformGate`; the Pauli-string
  operations (`PauliStringGateOperation`, `SingleQubitPauliStringGateOperation`,
  `MutablePauliString`, `MutableDensePauliString`, `BaseDensePauliString`,
  `PauliStringPhasor(Gate)`, `PauliMeasurementGate`, `PauliSumExponential`,
  `ProjectorString`, `ProjectorSum`, `LinearCombinationOfGates/Operations`,
  `PauliInteractionGate`, `measure_single_paulistring`,
  `measure_paulistring_terms`); the gate families (`AnyUnitaryGateFamily`,
  `AnyIntegerPowerGateFamily`, `ParallelGateFamily`, `ParallelGate`,
  `parallel_gate_op`, `InterchangeableQubitsGate`); the Clifford gates
  (`CliffordGate`, `SingleQubitCliffordGate`, `CXSWAP`, `CZSWAP`); the extra
  channels (`AsymmetricDepolarizingChannel`, `GeneralizedAmplitudeDampingChannel`,
  `KrausChannel`, `MixedUnitaryChannel`, `RandomGateChannel`,
  `StatePreparationChannel`); control (`ControlledOperation`,
  `AbstractControlValues`, `SumOfProducts`, `ProductOfSums`); classical control
  (`ClassicallyControlledOperation`, `If`, `SetVariable`, `Condition`,
  `KeyCondition`); the qubit managers (`QubitManager`, `SimpleQubitManager`,
  `GreedyQubitManager`, `BorrowableQubit`, `CleanQubit`); the variable qids
  (`VariableQid`, `VariableLineQid`, `VariableGridQid`); `ArithmeticGate` and
  `BooleanHamiltonianGate`; the tags `VirtualTag`/`RoutingSwapTag`; and the
  aliases `CSwapGate`, `M`, `R`, `QubitOrderOrList`, `PAULI_GATE_LIKE`,
  `PAULI_STRING_LIKE` and `PauliSumLike`. The factories `q`, `rx`/`ry`/`rz`,
  `ms`, `cphase`, `givens`, `riswap`, `wait`, `qft`, `reset_each`,
  `identity_each` and `global_phase_operation` round it out.
- **Quantum information, completed**: the state validators
  (`validate_qid_shape`, `validate_normalized_state_vector`,
  `validate_density_matrix`, `validate_indices`), the Choi/superoperator
  conversions (`kraus_to_choi`, `choi_to_kraus`, `kraus_to_superoperator`,
  `superoperator_to_kraus`, `choi_to_superoperator`, `superoperator_to_choi`,
  `operation_to_choi`, `operation_to_superoperator`),
  `density_matrix_from_state_vector`, `dirac_notation`, `one_hot`,
  `eye_tensor`, `entanglement_fidelity`, and the `QuantumState`/
  `QuantumStateRepresentation`/`quantum_state`/`StabilizerState` types.
- **Step-by-step simulation**: the simulator interfaces (`SimulatorBase`,
  `SimulatesSamples`, `SimulatesAmplitudes`, `SimulatesExpectationValues`,
  `SimulatesFinalState`, `SimulatesIntermediateState`,
  `SimulatesIntermediateStateVector`, `SimulationStateBase`, `StepResultBase`,
  `SimulationTrialResultBase`, `ClassicalStateSimulator`, `CIRCUIT_LIKE`), the
  concrete simulation states (`StateVectorSimulationState`,
  `DensityMatrixSimulationState`, `CliffordTableauSimulationState`,
  `StabilizerSimulationState`, `StabilizerChFormSimulationState`,
  `StabilizerStateChForm`, `SimulationState`, `SimulationProductState`), the
  step and trial results (`StepResult`, `StateVectorStepResult`,
  `DensityMatrixStepResult`, `CliffordSimulatorStepResult`,
  `SparseSimulatorStep`, `SimulationTrialResult`, `StateVectorTrialResult`,
  `CliffordTrialResult`, `StateVectorMixin`), `StabilizerSampler`, and the
  helpers `measure_state_vector`, `measure_density_matrix`,
  `final_state_vector`, `final_density_matrix`, `sample_state_vector`,
  `sample_density_matrix` and `sample_sweep`.
- **Transformers, completed**: the transformer decorator and
  `TransformerContext`/`TransformerLogger`, `merge_operations`/`merge_moments`
  (and their circuit-op variants), `unroll_circuit_op` (and greedy variants),
  `map_operations_and_unroll`, `remove_tags`/`index_tags`,
  `drop_terminal_measurements`, `drop_diagonal_before_measurement`,
  `dephase_measurements`, `eject_z`/`eject_phased_paulis`, `is_negligible_turn`,
  `defer_measurements`, the single-qubit merge-to-phased-x-and-z variants,
  `add_dynamical_decoupling`, `map_clean_and_borrowable_qubits`; the analytical
  decompositions (`single_qubit_matrix_to_pauli_rotations`,
  `single_qubit_op_to_framed_phase_form`, `decompose_multi_controlled_x`/
  `_rotation`, `two_qubit_matrix_to_sqrt_iswap_operations`/`_ion_operations`/
  `_cz_isometry`/`_diagonal_and_cz_operations`, `parameterized_2q_op_to_sqrt_iswap_operations`,
  `quantum_shannon_decomposition`, `three_qubit_matrix_to_operations`,
  `prepare_two_qubit_state_using_*`, `decompose_cphase_into_two_fsim`,
  `decompose_two_qubit_interaction_into_four_fsim_gates`,
  `TwoQubitGateTabulation`/`TwoQubitGateTabulationResult`); the compilation
  gatesets (`CompilationTargetGateset`, `TwoQubitCompilationTargetGateset`,
  `SqrtIswapTargetGateset`); and routing (`AbstractInitialMapper`,
  `HardCodedInitialMapper`, `LineInitialMapper`, `MappingManager`, `RouteCQC`,
  `routed_circuit_with_mapping`).
- **Devices and values, completed**: noise models (`NoiseModel`,
  `ConstantQubitNoiseModel`, `NoiseModelFromNoiseProperties`, `NoiseProperties`,
  `SuperconductingQubitsNoiseProperties`, `NO_NOISE`), `OpIdentifier`,
  `GridDeviceMetadata`, `UNCONSTRAINED_DEVICE`, the topologies (`NamedTopology`,
  `LineTopology`, `TiltedSquareLattice`) with `get_placements`/
  `is_valid_placement`/`draw_gridlike`/`draw_placements`; the `cirq.value`
  helpers (`big_endian_*`, `canonicalize_half_turns`, `chosen_angle_to_*`,
  `state_vector_to_probabilities`, `validate_probability`, `alternative`,
  `ABCMetaImplementAnyOneOf`, `LinearDict`), the `Duration`/`Timestamp`/
  `PeriodicValue`/`MeasurementType` values, the conditions
  (`BitMaskKeyCondition`, `SympyCondition`), the classical data stores
  (`ClassicalDataStore`, `ClassicalDataDictionaryStore`,
  `ClassicalDataStoreReader`), the `KET_*`/`PAULI_STATES` constants and
  `ProductState`.
- **Protocols and interop, completed**: the apply arguments
  (`ApplyUnitaryArgs`, `ApplyChannelArgs`, `ApplyMixtureArgs`), `act_on`,
  `apply_channel`, `apply_mixture`, `apply_unitaries`, the diagram helpers
  (`circuit_diagram_info`, `CircuitDiagramInfo`, `CircuitDiagramInfoArgs`),
  `DecompositionContext`/`decompose_once_with_qubits`, `DEFAULT_RESOLVERS`,
  `has_mixture`/`mixture`/`validate_mixture`, `has_stabilizer_effect`, the JSON
  helpers (`HasJSONNamespace`, `JsonResolver`, `json_namespace`,
  `dataclass_json_dict`, `obj_to_dict_helper`, `SerializableByKey`),
  `LabelEntity`, `measurement_key_name`, `mul`, `parameter_names`/
  `parameter_symbols`/`resolve_parameters_once`, `trace_distance_from_angle_list`,
  the key-path helpers (`with_key_path`, `with_key_path_prefix`,
  `with_rescoped_keys`) and the `Supports*` protocol classes; the Quirk import
  (`quirk_json_to_circuit`, `quirk_url_to_circuit`); `is_native_neutral_atom_gate`/
  `is_native_neutral_atom_op`; and more `vqubit.testing` helpers
  (`assert_implements_consistent_protocols`,
  `assert_has_consistent_apply_unitary`,
  `assert_decompose_is_consistent_with_unitary`,
  `assert_pauli_expansion_is_consistent_with_unitary`,
  `assert_qasm_is_consistent_with_unitary`, `assert_has_consistent_qid_shape`,
  `equals_tester`).
- **Full cirq-core surface**: the last circuits classes (`AbstractCircuit`,
  `Alignment`, `CircuitFunction`, `PointOptimizationSummary`, `PointOptimizer`,
  `QasmOutput`, `TextDiagramDrawer`), the study helpers (`Sweepable`, `symbol`,
  `flatten*`, `to_sweeps`, `ListSweep`, `ZipLongest`, `UnitSweep`/`UNIT_SWEEP`,
  `dict_to_*_sweep`, `list_of_dicts_to_zip`, `ResultDict`, `ExpressionMap`) and
  the experiments helpers (`linear_`/`log_xeb_fidelity*`, `xeb_fidelity`,
  `hog_score_xeb_fidelity_from_probabilities`,
  `estimate_*_readout_errors`, `measure_confusion_matrix`,
  `TensoredConfusionMatrices`). With these, **every one of cirq-core's 583
  top-level names is exported** (checked by a parity test).
- **Benchmarks**: a `benchmarks/compare.py` harness that runs the same tasks
  with vqubit and Cirq in fresh subprocesses and reports wall time and peak
  memory, plus a bilingual [benchmarks guide](docs/en/guide/benchmarks.md). The
  comparison dependency lives in the `bench` dependency group.
- **Packaging and tooling**: an MIT license, abi3 packaging covering CPython
  3.11-3.15 in a single wheel per platform, development tooling (`ruff`, `ty`,
  `pytest`), CI for build/test/docs/release, and bilingual documentation
  (Spanish and English) with runnable examples.

### Changed

- **Toolchain: built with vcraft 1.3.0.** The native extension now keeps its
  instance state inside the object instead of behind a pointer, so construction
  and native method calls are faster (a no-op `apply_gate` on a two-qubit state
  went from ~0.34 µs to ~0.31 µs) and a failed native constructor no longer leaks
  the object. The `NativeState`/`NativeDensity` buffers are exposed as typed
  `memoryview`s (`format "d"`, `itemsize 8`), and the generated CI build
  containers were regenerated for 1.3.0.
- **Performance: gate-heavy and density-matrix workloads are much faster.** The
  single-qubit power gates, `PhasedXPowGate` and `PhasedXZGate` now build their
  matrix in closed form (a `PhasedXZ` matrix went from ~34 µs to ~3 µs); gate
  unitaries and their native flat form are memoised on the gate, and
  `value_equality` caches the hash. On the density-matrix side the channel
  kernels were reworked: the 2x2 Kraus coefficients are hoisted, diagonal Kraus
  operators (CZ, CCZ, ...) get a direct scaling path, and the general block
  product is factored into two matrix products. A 250-layer, 12-qubit circuit
  now simulates in ~30 ms (Cirq: ~184 ms) and `density_large` is ~2x faster than
  Cirq.
- **Performance: vqubit now matches or beats Cirq on every benchmark** (see the
  [benchmarks guide](docs/en/guide/benchmarks.md)). The hot V kernels were
  rewritten (`apply_gate` as a block gather/apply/scatter with in-place and
  contiguous 1-/2-qubit paths plus a permutation fast path; reordered
  `matmul`/`sandwich`; an in-place subspace density channel), the simulator
  results are handed back as zero-copy NumPy views of the native buffers, the
  flatten/rebuild helpers, sampling, the Bloch vector and the circuit unitary
  use NumPy when it is installed, and `|0...0>` needs no Python list. vqubit
  imports ~2.7x less memory than Cirq (77 MB vs 208 MB).
- **Circuit construction is now O(1) amortised per operation.** `Circuit`
  caches each moment's qubit set and `append` scans backwards, stopping at the
  first conflicting moment, instead of re-deriving `Moment.qubits` for every
  earlier moment. Building a 300-moment circuit went from ~4.6 s to ~16 ms.
- **Sampling terminal measurements no longer re-simulates per repetition.** When
  every measurement is terminal, `Simulator.run` computes the final state once
  and draws each repetition's outcomes from it; sampling 1000 shots of a
  12-qubit circuit went from ~2.7 s to ~0.1 s.
- **Simulator state now lives in V.** The state-vector and density-matrix
  simulators and the Clifford simulator keep their state in `@[vc_buffer]`-backed
  native classes (`NativeState`, `NativeDensity`, `NativeTableau`) for the
  simulation's lifetime, so gates are applied to native memory without a Python
  round trip between calls. `memoryview` of a state is a zero-copy view.
- **Development tooling moved into `pyproject.toml`.** The dev and docs
  dependencies are `[dependency-groups]` and the ruff, pytest and ty settings are
  `[tool.*]` tables, installed with `uv pip install --group` (vcraft merges the
  file instead of overwriting it). The standalone requirements and config files
  are gone.

### Fixed

- **Numerical results are now cross-validated against Cirq.** A new harness
  (`benchmarks/verify.py`, run by the `cirq parity` CI job) compares gate
  unitaries, channels, circuit simulation, sampling, decompositions, Pauli and
  Clifford algebra, and the public API surface. It found and this release fixes:
  - `DensityMatrixSimulator` dropped the sign of diagonal two-qubit gates such
    as `CZ`/`CCZ` (the permutation fast path treated any real entry of modulus
    one as an unsigned permutation), so superposed states came out wrong.
  - `PhasedISwapPowGate` conjugated the wrong off-diagonal phase factor.
  - `MSGate`/`ms` used `exp(-i theta/2 XX)` instead of Cirq's `exp(-i rads XX)`;
    `MSGate` is now `MSGate(*, rads)`.
  - `QuantumFourierTransformGate`/`qft` used the wrong qubit reversal; it now
    returns the standard DFT and `qft` returns an `Operation`.
  - `DiagonalGate`/`TwoQubitDiagonalGate`/`ThreeQubitDiagonalGate` took the
    diagonal entries instead of the diagonal angles.
  - `QubitPermutationGate` permuted basis states instead of qubits.
  - `kak_vector` now accepts a batch of matrices and returns real values.
- **Signatures now match Cirq** for `two_qubit_matrix_to_cz_operations` and its
  isometry/diagonal variants, and for `merge_operations`/`merge_moments` (whose
  `merge_func` now takes two operations/moments).

[0.2.0]: https://github.com/rroblf01/vqubit/releases/tag/v0.2.0
[0.1.0]: https://github.com/rroblf01/vqubit/releases/tag/v0.1.0
