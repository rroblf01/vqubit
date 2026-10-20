# Example: the three-qubit code and stabilizers

The **three-qubit repetition code** is the smallest example of quantum error
correction: one logical qubit is spread over three physical ones, and a
*syndrome* says which one flipped without looking at the state. It is also the
way into the **stabilizer simulator**, which tracks the state as a *tableau*
instead of a vector of amplitudes.

## Encoding

```python
import vqubit as vq

data = vq.LineQubit.range(3)
encode = vq.Circuit(vq.CNOT(data[0], data[1]), vq.CNOT(data[0], data[2]))
encoded = vq.CliffordSimulator().simulate(encode).final_state
print([str(stabilizer) for stabilizer in encoded.stabilizers()])
```

```text
['Z(0)', 'Z(0)*Z(1)', 'Z(0)*Z(2)']
```

A stabilizer state is not stored as amplitudes: it is stored by the operators
that leave it unchanged. For `|000>` those are the `Z` of each qubit (and their
products). The code space also contains `|111>`:

```python
import numpy as np

one = vq.Circuit(vq.X(data[0]), *encode.all_operations())
print(np.round(np.asarray(vq.CliffordSimulator().simulate(one).final_state.state_vector()), 3))
```

```text
[0.j 0.j 0.j 0.j 0.j 0.j 0.j 1.j]
```

## An error and its syndrome

An `X` on qubit 1 does not leave the code space, but it **flips the sign of one
stabilizer**: that is the syndrome.

```python
flipped = vq.CliffordSimulator().simulate(vq.Circuit(*encode.all_operations(), vq.X(data[1]))).final_state
print([str(stabilizer) for stabilizer in flipped.stabilizers()])
```

```text
['Z(0)', '-1*Z(0)*Z(1)', 'Z(0)*Z(2)']
```

The `-1` on `Z(0)*Z(1)` says qubits 0 and 1 **disagree**, which is exactly the
pair the first stabilizer covers: the error is on one of the two (and the intact
second stabilizer tells which). Undoing it is applying the same `X` again:

```python
fixed = vq.CliffordSimulator().simulate(
    vq.Circuit(*encode.all_operations(), vq.X(data[1]), vq.X(data[1]))
).final_state
print([str(stabilizer) for stabilizer in fixed.stabilizers()])
```

```text
['Z(0)', 'Z(0)*Z(1)', 'Z(0)*Z(2)']
```

## Why the stabilizer simulator

A stabilizer state is described by `2n` rows of `2n + 1` bits, not by `2ⁿ`
amplitudes: a Clifford gate costs polynomially in `n`. These circuits —Hadamard,
`S`, `X`, `Y`, `Z`, `CNOT`, `CZ`, `SWAP`— can be simulated that way with far
more qubits than the state-vector simulator.

```python
print(vq.CliffordSimulator().simulate(vq.Circuit(vq.H(data[0]), vq.CNOT(data[0], data[1]))).final_state.dirac_notation())
```

```text
0.707|00⟩ + 0.707|11⟩
```

## What you have learned

- Stabilizers are the operators that leave the state unchanged, and their signs
  are the syndrome.
- `CliffordSimulator` simulates Clifford circuits exactly and polynomially.
- A detectable error does not leave the code space: it flips a sign you can read
  and correct.

For the state-vector simulator see [Simulation](../guide/simulation.md); for the
Clifford gates, [Quantum gates](../guide/quantum-gates.md).
