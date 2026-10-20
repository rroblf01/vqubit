# Example: Grover's search

**Grover's** algorithm finds a marked item in an unordered list with ~√N queries,
where a classical search would need N. With two qubits the list has four items
and **a single iteration** is enough to find the marked one every time.

## The circuit

We mark `|11>` and apply one iteration: the oracle (a `CZ`, which flips the sign
of `|11>`) and the diffuser (a reflection about the mean amplitude).

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
grover = vq.Circuit(
    vq.H(q0), vq.H(q1),                      # 1. uniform superposition
    vq.CZ(q0, q1),                            # 2. oracle: mark |11>
    vq.H(q0), vq.H(q1),                       # 3. diffuser: H, X, CZ, X, H
    vq.X(q0), vq.X(q1),
    vq.CZ(q0, q1),
    vq.X(q0), vq.X(q1),
    vq.H(q0), vq.H(q1),
)
```

The amplitude of each state after the circuit:

```python
import numpy as np

state = vq.Simulator().simulate(grover).state_vector()
print(np.round(np.abs(state), 3))
```

```text
[0. 0. 0. 1.]
```

All the amplitude sits on `|11>` (index 3), so measuring gives `11` every time:

```python
grover.append(vq.measure(q0, q1, key="m"))
print(vq.Simulator(seed=1).run(grover, repetitions=1000).histogram(key="m"))
```

```text
Counter({3: 1000})
```

## Why it works

The oracle flips the sign of the marked state and the diffuser reflects about
the mean amplitude. Together they are a rotation in the collective state space
that, with one iteration and four states, moves all the amplitude onto the
marked one. With more qubits it takes ~π√N/4 iterations, and the success
probability is no longer exact.

## What you have learned

- An oracle is a conditional phase: `CZ` already marks `|11>` with no extra
  gates.
- The diffuser is Hadamards and `X`s around that same `CZ`.
- The histogram of 1000 runs lands entirely on the marked state: the simulation
  is exact and the circuit deterministic.

For the gate catalogue see [Quantum gates](../guide/quantum-gates.md); for the
stabilizer simulator (which runs this circuit too), the
[three-qubit code example](stabilizers.md).
