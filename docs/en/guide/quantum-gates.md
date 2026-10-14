# Quantum gates

A **quantum gate** is an operation that changes the state of one or more qubits.
It is the quantum version of a classical logic operation, but reversible: every
gate (except measurement) can be undone.

In vqubit, gates are called on the qubits they act on:

```python
import vqubit as vq

q0 = vq.LineQubit(0)

vq.H(q0)      # Hadamard
vq.X(q0)      # quantum NOT (flips 0 <-> 1)
vq.measure(q0)
```

## The most common single-qubit gates

| Gate | Name       | What it does                                            |
| ---- | ---------- | ------------------------------------------------------- |
| `I`  | Identity   | Changes nothing.                                        |
| `X`  | Pauli X    | Flips the qubit, like a NOT: `0` → `1` and `1` → `0`.   |
| `Y`  | Pauli Y    | Like `X`, but with an imaginary phase factor.           |
| `Z`  | Pauli Z    | Flips the sign of the `1` part.                         |
| `H`  | Hadamard   | Creates or destroys superposition.                      |
| `S`  | Phase      | Adds a quarter turn of phase.                           |
| `T`  | π/8 phase  | Adds an eighth of a turn of phase.                      |

## Rotations

Instead of fixed gates, you can rotate a qubit by an angle you choose:

```python
import math

vq.Rx(q0, radians=math.pi / 2)   # half turn around the X axis
vq.Ry(q0, radians=math.pi / 2)   # around the Y axis
vq.Rz(q0, radians=math.pi / 4)   # around the Z axis
```

These are the gates you will need for real algorithms, where the angles are
computed at run time.

## Two-qubit gates

Entanglement needs gates that act on two qubits at once:

| Gate   | Name          | What it does                                                  |
| ------ | ------------- | ------------------------------------------------------------- |
| `CNOT` | Controlled-NOT | Flips the second qubit only if the first is `1`.             |
| `CZ`   | Controlled-Z  | Applies `Z` to the second qubit controlled by the first.      |
| `SWAP` | Swap          | Exchanges the states of the two qubits.                       |

```python
vq.CNOT(q0, q1)   # q0 controls, q1 is the target
vq.CZ(q0, q1)
vq.SWAP(q0, q1)
```

## Measurement

Measurement is not reversible: it collapses the state and returns a classical
bit.

```python
vq.measure(q0)                 # measure one qubit
vq.measure(q0, q1)             # measure several
vq.measure(q0, key="result")   # give the outcome a name
```

## How are gates represented?

Internally, each gate is a complex matrix (its **unitary**). For example, the
`X` gate is:

$$
X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
$$

You do not need to work with those matrices directly, but you can inspect them:

```python
matrix = vq.unitary(vq.X)
```

## Next steps

Combine gates into [circuits](circuits.md).
