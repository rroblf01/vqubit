# Example: expectation values and QAOA

An **expectation value** answers "what does this observable average to if I
repeat the experiment many times?". This example uses it to score a **QAOA**
(Quantum Approximate Optimization Algorithm) ansatz for MaxCut on a single edge,
and to look for the parameters that minimise it.

## The observable

For two qubits, the MaxCut cost is `Z ⊗ Z`: it is `+1` when the two qubits agree
and `-1` when they disagree. We want to minimise it.

```python
import math

import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
cost = vq.PauliString({q0: vq.Pauli.Z, q1: vq.Pauli.Z})
```

## The ansatz

One QAOA layer: a Hadamard on each qubit, the phase separator `exp(-i γ Z Z)`,
and an `Rx(2β)` mixer. The separator is written with `CZ`, because
`exp(-i γ Z Z)` and `CZ ** (2γ/π)` differ only by a global phase:

```python
def ansatz(gamma, beta):
    return vq.Circuit(
        vq.H(q0),
        vq.H(q1),
        (vq.CZ ** (2 * gamma / math.pi))(q0, q1),
        vq.rx(2 * beta)(q0),
        vq.rx(2 * beta)(q1),
    )

print(ansatz(math.pi / 2, math.pi / 8))
```

```text
0: ───H───@───Rx───
          │
1: ───H───@───Rx───
```

## The parameter grid

`simulate_expectation_values` takes a list of observables and returns one value
per observable:

```python
for gamma in (0.0, math.pi / 4, math.pi / 2):
    for beta in (0.0, math.pi / 8, math.pi / 4):
        value = vq.Simulator().simulate_expectation_values(
            ansatz(gamma, beta), observables=[cost]
        )[0].real
        print(f"gamma={gamma:.3f} beta={beta:.3f} -> <ZZ>={value:+.4f}")
```

```text
gamma=0.000 beta=0.000 -> <ZZ>=+0.0000
gamma=0.000 beta=0.393 -> <ZZ>=+0.0000
gamma=0.000 beta=0.785 -> <ZZ>=+0.0000
gamma=0.785 beta=0.000 -> <ZZ>=+0.0000
gamma=0.785 beta=0.393 -> <ZZ>=-0.2500
gamma=0.785 beta=0.785 -> <ZZ>=+0.5000
gamma=1.571 beta=0.000 -> <ZZ>=+0.0000
gamma=1.571 beta=0.393 -> <ZZ>=+0.5000
gamma=1.571 beta=0.785 -> <ZZ>=+1.0000
```

The minimum on the grid is `-0.25`, at `γ = π/4, β = π/8`: the combination that
brings the state closest to `|01>`/`|10>` (the two sides of the cut).

## A Hamiltonian with several terms

A real cost has many terms. They add with `+`, and the simulation accepts a
single `PauliString` or a whole `PauliSum`:

```python
hamiltonian = (
    vq.PauliString({q0: vq.Pauli.Z, q1: vq.Pauli.Z})
    + 0.5 * vq.PauliString({q0: vq.Pauli.Z})
    - 0.25
)
values = vq.Simulator().simulate_expectation_values(ansatz(math.pi / 4, math.pi / 8), observables=[hamiltonian])
print(f"<H> = {values[0].real:+.4f}")
```

```text
<H> = -0.3232
```

## What you have learned

- `PauliString` is an observable; `+`, `*` and `-` build a `PauliSum`.
- `simulate_expectation_values(circuit, observables=[...])` returns a list of
  complex numbers; the imaginary part is zero for a Hermitian observable.
- Powers of `CZ` are phase separators: `CZ ** t` applies
  `exp(i π t (1 - Z Z) / 2)`.

For the simulator details see [Simulation](../guide/simulation.md); for the list
of observables, [Quantum gates](../guide/quantum-gates.md).
