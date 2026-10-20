# Ejemplo: valores esperados y QAOA

Un **valor esperado** responde a la pregunta "¿cuánto vale este observable, de
media, si repito el experimento muchas veces?". En este ejemplo lo usamos para
evaluar un ansatz de **QAOA** (Quantum Approximate Optimization Algorithm) para
MaxCut sobre una arista, y buscamos los parámetros que lo minimizan.

## El observable

Para dos qubits, el coste de MaxCut es `Z ⊗ Z`: vale `+1` si los dos qubits
están de acuerdo y `-1` si están en desacuerdo. Queremos minimizarlo.

```python
import math

import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
cost = vq.PauliString({q0: vq.Pauli.Z, q1: vq.Pauli.Z})
```

## El ansatz

Una capa de QAOA: Hadamard en los dos qubits, un separador de fase
`exp(-i γ Z Z)` y un mezclador `Rx(2β)` en cada uno. El separador se escribe con
`CZ`, porque `exp(-i γ Z Z)` y `CZ ** (2γ/π)` solo difieren en una fase global:

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

## La rejilla de parámetros

`simulate_expectation_values` acepta una lista de observables y devuelve una
lista de valores, uno por observable:

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

El mínimo de la rejilla es `-0.25`, en `γ = π/4, β = π/8`: es la combinación que
más acerca el estado a `|01>`/`|10>` (los dos cortes de la arista).

## Un Hamiltoniano de varios términos

Un coste real tiene muchos términos. Se suman con `+` y la simulación acepta
tanto una `PauliString` como una `PauliSum`:

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

## Qué has aprendido

- `PauliString` es un observable; `+`, `*` y `-` construyen una `PauliSum`.
- `simulate_expectation_values(circuito, observables=[...])` devuelve una lista
  de complejos; la parte imaginaria es cero para un observable hermítico.
- Las potencias de `CZ` son separadores de fase: `CZ ** t` aplica
  `exp(i π t (1 - Z Z) / 2)`.

Para los detalles del simulador, mira [Simulación](../guide/simulation.md); para
la lista de observables, [Puertas cuánticas](../guide/quantum-gates.md).
