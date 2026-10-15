# Optimización

Un **transformador** toma un circuito y devuelve otro equivalente, más limpio o
más barato. Los optimizadores viven en el espacio de nombres `vqubit`.

## Limpiar circuitos

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
c = vq.Circuit(vq.H(q0), vq.I(q0), vq.CNOT(q0, q1))

vq.drop_empty_moments(c)              # quita momentos vacíos
vq.drop_negligible_operations(c)      # quita puertas que no hacen nada
vq.expand_composite(c)                # expande puertas compuestas (CCX, ...)
```

## Agrupar operaciones

`merge_single_qubit_moments_to_phxz` combina rachas de puertas de un solo qubit
en una sola `PhasedXZGate`, que es la puerta nativa de muchos procesadores:

```python
c = vq.Circuit(vq.X(q0), vq.H(q0), vq.Rz(0.3)(q0))
c = vq.merge_single_qubit_moments_to_phxz(c)
print(c)   # una sola puerta PhXZ
```

## Descomponer un unitario de un qubit

Cualquier puerta de un qubit se puede escribir como una secuencia de rotaciones:

```python
matriz = vq.unitary(vq.H)
print(vq.single_qubit_matrix_to_gates(matriz))          # [Rz, Ry, Rz]
print(vq.single_qubit_matrix_to_phased_x_z(matriz))     # [PhasedX, Z]
print(vq.single_qubit_matrix_to_phxz(matriz))           # una PhasedXZGate
```

## Compilar a un conjunto de puertas

Un **conjunto objetivo** (target gate set) dice qué puertas puede ejecutar un
dispositivo. `optimize_for_target_gateset` reescribe el circuito para usarlas:

```python
q0, q1 = vq.LineQubit.range(2)
c = vq.Circuit(vq.CNOT(q0, q1), vq.H(q0))

compilado = vq.optimize_for_target_gateset(c, vq.CZTargetGateset())
print(compilado)   # solo puertas CZ y PhasedXZ
```

Una puerta de dos qubits arbitraria se sintetiza con la descomposición **KAK**:

```python
import vqubit as vq
from vqubit import two_qubit_matrix_to_cz_operations

q0, q1 = vq.LineQubit.range(2)
ops = two_qubit_matrix_to_cz_operations(q0, q1, vq.unitary(vq.ISWAP))
```

## Siguientes pasos

- Guarda circuitos con [serialización](serialization.md).
- Ejecútalos en [dispositivos](devices.md) locales.
