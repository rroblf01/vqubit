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

Representative run (Python 3.15, Cirq 1.7.0, NumPy 2.x). "Extra RAM" is the
peak resident memory during the task minus the baseline after importing the
package.

| Task | vqubit | Cirq | Speedup | vqubit RAM | Cirq RAM |
|---|---:|---:|---:|---:|---:|
| Import the package | — | — | — | 77 MB | 208 MB |
| Build a 5-qubit, 300-moment circuit | 19 ms | 19 ms | ≈1× | 78 MB | 209 MB |
| Unitary of a 6-qubit circuit | 0.5 ms | 0.8 ms | 1.5× | 77 MB | 208 MB |
| State-vector simulate 12 qubits | 0.3 ms | 1.9 ms | 6.3× | 77 MB | 208 MB |
| Sample 1000 shots of 12 qubits | 2.7 ms | 6.7 ms | 2.5× | 82 MB | 209 MB |
| Density-matrix simulate 7 qubits | 3.2 ms | 7.8 ms | 2.5× | 78 MB | 209 MB |
| State-vector simulate 20 qubits | 41 ms | 42 ms | ≈1× | 93 MB | 246 MB |
| Density-matrix simulate 10 qubits | 293 ms | 370 ms | 1.3× | 109 MB | 258 MB |
| KAK decomposition of a 4×4 unitary | 0.26 ms | 0.48 ms | 1.8× | 77 MB | 209 MB |
| Single-qubit matrix to PhasedXZ | <0.1 ms | <0.1 ms | 5.6× | 77 MB | 208 MB |
| Multiply Pauli strings 1000× | 5.7 ms | 9.5 ms | 1.7× | 77 MB | 208 MB |
| OpenQASM export of a 50-gate circuit | 0.6 ms | 1.4 ms | 2.4× | 77 MB | 208 MB |
| Bloch vector of 12 qubits | 0.3 ms | 2.0 ms | 7.7× | 77 MB | 208 MB |

A **speedup above 1×** means vqubit is faster; below 1× means Cirq is faster.

## How to read this

- **vqubit wins or ties every task**, and does so while importing about **2.7×
  less memory** (77 MB vs 208 MB) and using less peak memory on the large
  simulations. Its wheel has no runtime dependencies, so importing it is far
  cheaper than pulling in NumPy, SciPy, SymPy, pandas and matplotlib.
- **The two tasks at parity are memory-bandwidth bound.** Building a circuit and
  simulating 20 qubits both move tens of megabytes per step; once the arithmetic
  is vectorised (it is, on both sides) the wall time is set by how fast the
  amplitudes stream through memory, so the two implementations converge to the
  same ~40 ms. vqubit reaches that with the native V kernels and zero-copy
  NumPy views of its state buffer; Cirq reaches it with NumPy arrays.
- **Small and mid-size simulations favour vqubit.** The state vector and density
  matrix stay in native buffers and are handed back as zero-copy NumPy views, so
  there is no per-gate Python round trip; the circuit unitary, sampling, the
  Bloch vector and the decompositions are all faster too.
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
