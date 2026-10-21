# Migrating from Cirq

If you already program with Cirq, almost everything you know still applies:
vqubit copies its names, its signatures and its behaviour. This guide is the
checklist for porting a real program, in the order the problems show up.

!!! tip "The short version"
    Replace `import cirq` with `import vqubit as vq` and run it. Then: use the
    full argument names (`probability=`, `radians=`, `num_qubits=`), treat
    `num_qubits` as an **attribute** on gates, and review the four behavioural
    differences in the last section.

## 1. The import

```python
# before
import cirq

# after
import vqubit as vq
```

vqubit exposes **606 of the 607 top-level names** of cirq-core 1.7.0 (the one
missing is `contrib`). That includes the submodules that mirror Cirq's:

```python
from vqubit import ops, sim, linalg, protocols, transformers, work
```

Where Cirq had `cirq.ops.X`, vqubit has `vq.ops.X`; the flat name `vq.X` also
exists and is the same object.

## 2. Arguments change their names

vqubit spells parameters out. Positional calls are unaffected; keyword calls
are not:

| Cirq                     | vqubit                          |
| ------------------------ | ------------------------------- |
| `p=0.1`                  | `probability=0.1`               |
| `rads=math.pi / 2`       | `radians=math.pi / 2`           |
| `n_qubits=2`             | `num_qubits=2`                  |
| `val=vq.X`               | `value=vq.X`                    |
| `mat=u`                  | `matrix=u`                      |
| `contents=[...]`         | `operations=[...]`              |
| `tolerance=1e-8`         | `atol=1e-8`                     |

A wrong name raises `TypeError` rather than doing something else quietly.

## 3. `num_qubits` is an attribute on gates

```python
# Cirq
cirq.H.num_qubits()      # method
cirq.H(q0).num_qubits()  # method

# vqubit
vq.H.num_qubits          # attribute
vq.H(q0).num_qubits()    # operations keep the method
```

The rule: **gate → attribute**, **operation → method**. It is the most common
change when porting loops over operations.

## 4. Measuring and sampling

`sample` returns a `vqubit.Result`, not a `pandas.DataFrame`. If you only want
the counts, use `histogram`:

```python
result = vq.Simulator().run(circuit, repetitions=1000)
print(result.histogram(key="m"))          # Counter({0: 512, 1: 488})
rows = result.measurements["m"]           # list of rows
```

`run_sweep`, `simulate_sweep`, `run_batch`, `simulate_moment_steps` and
`simulate_expectation_values` take the same arguments as in Cirq, with the
parameter resolver in **second** position.

## 5. Parameters are not SymPy symbols

```python
# Cirq
import sympy
t = sympy.Symbol("t")

# vqubit
t = vq.Parameter("t")
```

Arithmetic (`t / 2`, `2 * t + 1`, `-t`), `resolve_parameters` and symbolic
gates work the same way, but equality is vqubit's (two expressions are equal
when they are structurally equal). If you really need SymPy, `vqubit.interop`
documents the conversion: see [Serialization](serialization.md).

## 6. The plain gates are the power gates

```python
vq.CZ is vq.CZPowGate()        # True, as in Cirq
vq.CZ ** 0.5                   # CZPowGate(exponent=0.5)
vq.CNOT ** -1 == vq.CNOT       # True: equality is periodic
vq.ISWAP ** -1 != vq.ISWAP     # ISWAP does distinguish -1
```

The historical aliases (`vq.CZGate`, `vq.CNOTGate`, `vq.SWAPGate`, `vq.CCXGate`,
`vq.ISWAPGate`) still exist, so an `isinstance` check keeps working.

## 7. What is not there

- **`cirq.contrib`**: it does not exist (its QASM importer depends on `ply`).
  Instead, vqubit ships `vq.circuit_from_qasm`, dependency-free.
- **The `EigenGate` tree**: multi-qubit power gates derive from `_PowGate`, not
  from `EigenGate`.
- **Qudit `SimulationState`s**: vqubit rejects non-qubit shapes explicitly in
  the measurement helpers.

Everything else —qubits, gates, circuits, simulators, transformers, devices,
`cirq.work`— is there, with the same design.

## 8. The four behavioural differences

They are explained, with examples and closing criteria, in
[Differences from Cirq](cirq-differences.md). In short:

1. **Text diagrams**: two-qubit Clifford gates are drawn with their name
   (`CliffordGate(2q)`), and exponents are decimal (`^0.5`).
2. **Class equality**: it requires the same concrete class (Cirq compares
   against the decorated class). Only `SqrtISwapGate` shows it, and that class
   has no Cirq counterpart.
3. **`sample`**: returns `Result`, not a `DataFrame`.
4. **`Parameter`**: is not a `sympy.Symbol`.

## 9. A program ported end to end

The same program in both libraries: prepare a three-qubit GHZ state, measure it
and estimate the expectation value of `XXX`.

```python
# Cirq
import cirq

q0, q1, q2 = cirq.LineQubit.range(3)
circuit = cirq.Circuit(cirq.H(q0), cirq.CNOT(q0, q1), cirq.CNOT(q1, q2))
circuit.append(cirq.measure(q0, q1, q2, key="m"))
result = cirq.Simulator().run(circuit, repetitions=100)
print(result.histogram(key="m"))
observable = cirq.X(q0) * cirq.X(q1) * cirq.X(q2)
print(cirq.Simulator().simulate_expectation_values(circuit, observables=[observable]))
```

```python
# vqubit
import vqubit as vq

q0, q1, q2 = vq.LineQubit.range(3)
circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.CNOT(q1, q2))
circuit.append(vq.measure(q0, q1, q2, key="m"))
result = vq.Simulator().run(circuit, repetitions=100)
print(result.histogram(key="m"))
observable = vq.PauliString({q0: vq.Pauli.X, q1: vq.Pauli.X, q2: vq.Pauli.X})
print(vq.Simulator().simulate_expectation_values(circuit, observables=[observable]))
```

Both print a histogram with only `0` and `7` (the two GHZ states) and an
expectation value of `1.0`: `X ⊗ X ⊗ X` is the stabilizer of the state, while
`Z ⊗ Z ⊗ Z` would give `0` (the two GHZ terms have opposite sign).

## 10. Checking parity

If you doubt a specific call, the repository ships two harnesses that compare
vqubit with cirq-core 1.7.0:

```bash
python benchmarks/parity.py    # 1867 API checks
python benchmarks/verify.py    # 650 numerical checks
```

`benchmarks/parity_tolerated.toml` is the closed list of deliberate
divergences: if something is not there, it is a bug.
