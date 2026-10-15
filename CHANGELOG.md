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
- **Packaging and tooling**: an MIT license, abi3 packaging covering CPython
  3.11-3.15 in a single wheel per platform, development tooling (`ruff`, `ty`,
  `pytest`), CI for build/test/docs/release, and bilingual documentation
  (Spanish and English) with runnable examples.

[Unreleased]: https://github.com/rroblf01/vqubit/commits/main
