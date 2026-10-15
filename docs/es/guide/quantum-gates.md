# Puertas cuánticas

Una **puerta cuántica** es una operación que cambia el estado de uno o varios
qubits. Es la versión cuántica de las operaciones lógicas de un ordenador
clásico, pero reversible: toda puerta (salvo la medición) puede deshacerse.

En vqubit, las puertas se llaman con los qubits sobre los que actúan:

```python
import vqubit as vq

q0 = vq.LineQubit(0)
q1 = vq.LineQubit(1)

vq.H(q0)      # Hadamard
vq.X(q0)      # NOT cuántico (invierte 0 <-> 1)
vq.measure(q0)
```

## Las puertas de un solo qubit más comunes

| Puerta | Nombre          | Qué hace                                                    |
| ------ | --------------- | ----------------------------------------------------------- |
| `I`    | Identidad       | No cambia nada.                                             |
| `X`    | Pauli X         | Invierte el qubit, como un NOT: `0` → `1` y `1` → `0`.      |
| `Y`    | Pauli Y         | Como `X`, pero con un factor de fase imaginario.            |
| `Z`    | Pauli Z         | Cambia el signo de la parte `1`.                            |
| `H`    | Hadamard        | Crea o destruye superposición.                              |
| `S`    | Fase            | Añade un cuarto de vuelta de fase.                          |
| `T`    | Fase π/8        | Añade un octavo de vuelta de fase.                          |

## Rotaciones

En lugar de puertas fijas, puedes girar un qubit un ángulo que elijas:

```python
import math

vq.Rx(math.pi / 2)(q0)   # media vuelta alrededor del eje X
vq.Ry(math.pi / 2)(q0)   # alrededor del eje Y
vq.Rz(math.pi / 4)(q0)   # alrededor del eje Z
```

Estas son las puertas que necesitarás para algoritmos reales, donde los ángulos
se calculan en tiempo de ejecución.

## Puertas de dos qubits

Para crear entrelazamiento hacen falta puertas que actúen sobre dos qubits a la
vez:

| Puerta  | Nombre             | Qué hace                                                        |
| ------- | ------------------ | --------------------------------------------------------------- |
| `CNOT`  | Controlled-NOT     | Invierte el segundo qubit solo si el primero está en `1`.       |
| `CZ`    | Controlled-Z       | Aplica `Z` al segundo qubit controlado por el primero.          |
| `SWAP`  | Intercambio        | Intercambia los estados de los dos qubits.                      |

```python
vq.CNOT(q0, q1)   # q0 controla, q1 es el objetivo
vq.CZ(q0, q1)
vq.SWAP(q0, q1)
```

## Medición

La medición no es reversible: colapsa el estado y devuelve un bit clásico.

```python
vq.measure(q0)                 # mide un qubit
vq.measure(q0, q1)             # mide varios
vq.measure(q0, key="result")   # le da un nombre al resultado
```

## ¿Cómo se representan las puertas?

Internamente, cada puerta es una matriz compleja (su **unitario**). Por ejemplo,
la puerta `X` es:

$$
X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
$$

No necesitas trabajar con esas matrices directamente, pero puedes consultarlas:

```python
matriz = vq.unitary(vq.X)
```

## Siguientes pasos

Combina puertas para formar [circuitos](circuits.md).
