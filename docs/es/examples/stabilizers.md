# Ejemplo: el código de 3 qubits y los estabilizadores

El **código de repetición de 3 qubits** es el ejemplo más pequeño de corrección
de errores cuánticos: un qubit de información se reparte en tres qubits físicos
y un *síndrome* dice cuál se ha dado la vuelta, sin mirar el estado. Es también
la puerta de entrada al **simulador de estabilizadores**, que sigue el estado
como un *tableau* en vez de como un vector de amplitudes.

## Codificar

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

Un estado de estabilizador no se guarda como amplitudes: se guarda por los
operadores que lo dejan igual. Para `|000>` son los `Z` de cada qubit (y sus
productos). El código también contiene `|111>`:

```python
import numpy as np

one = vq.Circuit(vq.X(data[0]), *encode.all_operations())
print(np.round(np.asarray(vq.CliffordSimulator().simulate(one).final_state.state_vector()), 3))
```

```text
[0.j 0.j 0.j 0.j 0.j 0.j 0.j 1.j]
```

## Un error y su síndrome

Un `X` en el qubit 1 no cambia el estado del código, pero **sí cambia el signo
de un estabilizador**: eso es el síndrome.

```python
flipped = vq.CliffordSimulator().simulate(vq.Circuit(*encode.all_operations(), vq.X(data[1]))).final_state
print([str(stabilizer) for stabilizer in flipped.stabilizers()])
```

```text
['Z(0)', '-1*Z(0)*Z(1)', 'Z(0)*Z(2)']
```

El `-1` en `Z(0)*Z(1)` dice que los qubits 0 y 1 **discrepan**, que es exactamente
el par que comparten el primer estabilizador: el error está en uno de los dos (y
el otro estabilizador, intacto, lo desambigua). Deshacerlo es aplicar el mismo
`X` otra vez:

```python
fixed = vq.CliffordSimulator().simulate(
    vq.Circuit(*encode.all_operations(), vq.X(data[1]), vq.X(data[1]))
).final_state
print([str(stabilizer) for stabilizer in fixed.stabilizers()])
```

```text
['Z(0)', 'Z(0)*Z(1)', 'Z(0)*Z(2)']
```

## Por qué el simulador de estabilizadores

Un estado de estabilizador se describe con `2n` filas de `2n + 1` bits, no con
`2ⁿ` amplitudes: el coste de una puerta Clifford es polinómico en `n`. Estos
circuitos —Hadamard, `S`, `X`, `Y`, `Z`, `CNOT`, `CZ`, `SWAP`— pueden simularse
así con muchos más qubits que con el vector de estado.

```python
print(vq.CliffordSimulator().simulate(vq.Circuit(vq.H(data[0]), vq.CNOT(data[0], data[1]))).final_state.dirac_notation())
```

```text
0.707|00⟩ + 0.707|11⟩
```

## Qué has aprendido

- Los estabilizadores son los operadores que dejan el estado invariante, y sus
  signos son el síndrome.
- `CliffordSimulator` simula circuitos Clifford de forma exacta y polinómica.
- Un error detectable no rompe el estado del código: cambia un signo que se puede
  leer y corregir.

Para el simulador de vectores de estado, mira [Simulación](../guide/simulation.md);
para las puertas Clifford, [Puertas cuánticas](../guide/quantum-gates.md).
