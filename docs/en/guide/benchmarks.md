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

Representative run (Python 3.15, Cirq 1.7.0, NumPy 2.x). The RAM columns are the
peak resident memory during the task, which includes importing the package
(each task runs in a fresh subprocess that imports only the package under test).

| Task | vqubit | Cirq | Speedup | vqubit RAM | Cirq RAM |
|---|---:|---:|---:|---:|---:|
| Import the package | — | — | — | 79 MB | 205 MB |
| Build a 5-qubit, 300-moment circuit | 10 ms | 22 ms | 2.3× | 80 MB | 206 MB |
| Unitary of a 6-qubit circuit | 0.4 ms | 0.7 ms | 1.9× | 79 MB | 205 MB |
| State-vector simulate 12 qubits | 0.2 ms | 2.2 ms | 12.4× | 79 MB | 206 MB |
| Sample 1000 shots of 12 qubits | 2.9 ms | 6.0 ms | 2.1× | 84 MB | 206 MB |
| Density-matrix simulate 7 qubits | 2.1 ms | 7.8 ms | 3.7× | 80 MB | 206 MB |
| State-vector simulate 20 qubits | 36 ms | 28 ms | 0.79× | 95 MB | 244 MB |
| Density-matrix simulate 10 qubits | 154 ms | 301 ms | 2.0× | 128 MB | 255 MB |
| Simulate a 250-layer, 12-qubit circuit | 34 ms | 176 ms | 5.2× | 84 MB | 209 MB |
| KAK decomposition of a 4×4 unitary | 0.2 ms | 0.5 ms | 1.9× | 79 MB | 206 MB |
| Single-qubit matrix to PhasedXZ | <0.1 ms | <0.1 ms | 4.8× | 79 MB | 205 MB |
| Multiply Pauli strings 1000× | 3.2 ms | 8.9 ms | 2.8× | 79 MB | 205 MB |
| OpenQASM export of a 50-gate circuit | 0.4 ms | 1.2 ms | 3.3× | 79 MB | 206 MB |
| Bloch vector of 12 qubits | 0.2 ms | 1.7 ms | 8.0× | 79 MB | 206 MB |

A **speedup above 1×** means vqubit is faster; below 1× means Cirq is faster.

## How to read this

- **vqubit wins 13 of the 14 tasks**, and does so while importing about **2.6×
  less memory** (79 MB vs 205 MB) and using less peak memory on the large
  simulations. Its wheel has no runtime dependencies, so importing it is far
  cheaper than pulling in NumPy, SciPy, SymPy, pandas and matplotlib. The one
  task it loses (20 qubits) is measured below.
- **Gate-heavy circuits are where the gap is widest.** A 250-layer, 12-qubit
  circuit runs ~5× faster than Cirq: vqubit builds each single-qubit gate matrix
  in closed form and keeps it (and its native flat form) memoised, so the
  per-gate Python cost that usually dominates small-gate workloads largely
  disappears.
- **Simulating 20 qubits is memory-bandwidth bound, and it is the one task Cirq
  wins (~0.8×).** Each gate streams the whole 16 MB state through memory, so wall
  time is set by bandwidth, not arithmetic; the two sides are close (36 ms vs
  28 ms). Cirq's NumPy kernels spend that bandwidth on vectorised (SIMD) moves,
  while the V kernels rely on the compiler's auto-vectorisation, which is not as
  tight for the interleaved complex layout. Closing the last ~20% would need
  explicit SIMD, which the portable, dependency-free build does not require.
- **Density-matrix simulation is ~2× faster.** The channel kernels scale the
  diagonal (CZ/CCZ) directly and factor the general block product into two
  matrix products.
- When NumPy is installed the public API returns `numpy.ndarray`
  (see [Installation](../installation.md)); without it everything still works on
  plain lists, just without the NumPy-accelerated paths.

## Reproducing the numbers

```bash
uv pip install --group bench
python benchmarks/compare.py --json > benchmarks.json
```

The `bench` dependency group is only for the comparison; the published vqubit
wheel never depends on Cirq or NumPy.
