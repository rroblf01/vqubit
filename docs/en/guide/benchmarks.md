# Benchmarks

This page compares vqubit and [Cirq](https://quantumai.google/cirq) on the same
tasks, so you can see the trade-offs in speed and memory between a
dependency-free, pure-Python library and a NumPy-based one.

The harness lives in `benchmarks/compare.py` and only uses the standard library.
It runs each task in a **fresh subprocess** (so peak memory is not polluted by
earlier tasks) and reports:

- **time**: the best wall time over several runs;
- **extra RAM**: the peak resident memory (`VmHWM`) during the task minus the
  baseline after importing the package.

## Running it

```bash
# Install the comparison dependency (Cirq) into your environment.
uv pip install --group bench

# Print the table (Cirq column is blank if Cirq is not installed).
python benchmarks/compare.py

# Or get the raw numbers as JSON.
python benchmarks/compare.py --json
```

## Methodology

Each task is implemented twice — once with vqubit and once with Cirq — doing the
same amount of work. Every task runs in its own process:

1. import the package and run the task once to warm up;
2. record the baseline RSS;
3. run the task `N` times, keeping the fastest time;
4. read `VmHWM` (peak RSS) for the memory figure.

Times are wall-clock seconds on a single core. The machine used for the numbers
below is an ordinary Linux x86-64 development machine; **run the harness on your
own hardware** for figures that matter to you.

## Results

Representative run (Python 3.15, Cirq 1.7.0, NumPy 2.x):

| Task | vqubit | Cirq | Speedup | vqubit extra RAM | Cirq extra RAM |
|---|---:|---:|---:|---:|---:|
| Build a 5-qubit, 300-moment circuit | 18.4 ms | 18.7 ms | 1.02× | 0 MB | 0 MB |
| Unitary of a 6-qubit circuit | 45.7 ms | 0.7 ms | 0.02× | 1 MB | 0 MB |
| State-vector simulate 12 qubits | 4.0 ms | 1.7 ms | 0.41× | 0 MB | 0 MB |
| Sample 1000 shots of 12 qubits | 96.7 ms | 6.5 ms | 0.07× | 0 MB | 0 MB |
| Density-matrix simulate 7 qubits | 22.7 ms | 8.3 ms | 0.37× | 7 MB | 0 MB |
| State-vector simulate 20 qubits | 1453 ms | 42 ms | 0.03× | 199 MB | 39 MB |
| Density-matrix simulate 10 qubits | 2366 ms | 387 ms | 0.16× | 593 MB | 49 MB |
| KAK decomposition of a 4×4 unitary | 0.3 ms | 0.9 ms | 3.2× | 0 MB | 0 MB |
| Single-qubit matrix to PhasedXZ | <0.1 ms | <0.1 ms | 5.6× | 0 MB | 0 MB |
| Multiply Pauli strings 1000× | 6.2 ms | 9.7 ms | 1.6× | 0 MB | 0 MB |
| OpenQASM export of a 50-gate circuit | 0.5 ms | 1.3 ms | 2.5× | 0 MB | 0 MB |
| Bloch vector of 12 qubits | 8.2 ms | 1.7 ms | 0.21× | 0 MB | 0 MB |

A **speedup above 1×** means vqubit is faster; below 1× means Cirq is faster.

## How to read this

- **Circuit building, decompositions, Pauli algebra and QASM** are where vqubit
  is competitive or faster. These are pure-Python, allocation-light operations,
  and vqubit's compact representations pay off.
- **Large simulations and the unitary of a whole circuit** favour Cirq. vqubit's
  simulators do the linear algebra in the native V kernels, but the result is
  materialised as Python `complex` objects, which costs time and memory at 20
  qubits. Cirq keeps everything in NumPy arrays, so its conversions are cheap.
  When NumPy is installed, vqubit's own `unitary`/`kraus`/state accessors return
  `numpy.ndarray` (see [Installation](../installation.md)), which narrows the
  gap for consumers but does not change the simulator's internal cost.
- **Memory** is dominated by the import footprint for small tasks (Cirq pulls in
  NumPy, SciPy, SymPy, pandas and matplotlib), and by the simulation buffer for
  large tasks. vqubit's wheel has no runtime dependencies, so importing it is
  far cheaper; the harness reports the *extra* memory during a task so the
  import cost does not hide the simulation cost.

## Reproducing the numbers

```bash
uv pip install --group bench
python benchmarks/compare.py --json > benchmarks.json
```

The `bench` dependency group is only for the comparison; the published vqubit
wheel never depends on Cirq or NumPy.
