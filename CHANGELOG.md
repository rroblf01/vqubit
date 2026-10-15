# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

The first release: a dependency-free, Cirq-compatible quantum computing
framework for Python, with the numerical core written in V.

### Changed

- **Simulator state now lives in V.** The state-vector and density-matrix
  simulators and the Clifford simulator keep their state in `@[vc_buffer]`-backed
  native classes (`NativeState`, `NativeDensity`, `NativeTableau`) for the
  simulation's lifetime, so gates are applied to native memory without a Python
  round trip between calls. `memoryview` of a state is a zero-copy view.
- **Development tooling moved into `pyproject.toml`.** The dev and docs
  dependencies are `[dependency-groups]` and the ruff, pytest and ty settings are
  `[tool.*]` tables, installed with `uv pip install --group` (vcraft 1.2.0 merges
  the file instead of overwriting it). The standalone requirements and config
  files are gone.

### Added

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
  returns `numpy.ndarray` results and accepts numpy inputs; the core stays
  dependency-free otherwise.
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
- **Packaging and tooling**: an MIT license, abi3 packaging covering CPython
  3.11-3.15 in a single wheel per platform, development tooling (`ruff`, `ty`,
  `pytest`), CI for build/test/docs/release, and bilingual documentation
  (Spanish and English) with runnable examples.

[Unreleased]: https://github.com/rroblf01/vqubit/commits/main
