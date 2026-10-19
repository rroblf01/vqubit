# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`**` on the plain gates.** `vq.CZ ** 0.5`, `vq.CNOT ** 0.5`,
  `vq.SWAP ** 0.5`, `vq.ISWAP ** 0.5`, `vq.SQRT_ISWAP ** 0.5`, `vq.CCX ** 0.5`,
  `vq.TOFFOLI ** 0.5` and `vq.I ** 0.5` now return the matching power gate, as
  the Cirq idioms do, and `vq.rz(0.3) ** 0.5` scales the angle to `vq.rz(0.15)`.
  `vq.FREDKIN` and `vq.CSWAP` stay non-powerable, like `cirq.CSWAP`. The plain
  gates themselves are the power gates now: see *Changed*.
- **A Cirq contract harness**, `benchmarks/parity.py`, run in the `cirq-parity`
  CI job: it checks names, `**`, `isinstance`, gate classes and text diagrams
  against Cirq, and keeps the accepted divergences in
  `benchmarks/parity_tolerated.toml`, where an entry that stops silencing a real
  divergence fails the build.

### Changed

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
