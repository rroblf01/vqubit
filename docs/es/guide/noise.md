# Ruido

Un **canal** es una puerta que no es unitaria: describe ruido y lleva un estado
puro a uno mezclado. vqubit define los canales con sus **operadores de Kraus**,
las matrices $K_i$ tales que el estado evoluciona como
$\rho \to \sum_i K_i \rho K_i^\dagger$.

## Canales incluidos

| Canal | Constructor | Qué modela |
| ----- | ----------- | ---------- |
| Bit flip | `vq.bit_flip(p)` | Invierte el qubit con probabilidad `p`. |
| Phase flip | `vq.phase_flip(p)` | Aplica `Z` con probabilidad `p`. |
| Bit-phase flip | `vq.bit_phase_flip(p)` | Aplica `Y` con probabilidad `p`. |
| Depolarizante | `vq.depolarize(p, n=1)` | Aplica un Pauli aleatorio con probabilidad `p`. |
| Amortiguamiento de amplitud | `vq.amplitude_damp(gamma)` | Pérdida de energía (`|1>` → `|0>`). |
| Amortiguamiento de fase | `vq.phase_damp(gamma)` | Pérdida de coherencia. |
| Reset | `vq.reset(n=1)` | Lleva los qubits a `|0>`. |

## Usarlos en un circuito

Los canales se colocan como cualquier operación:

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuito = vq.Circuit(
    vq.H(q0),
    vq.depolarize(0.05)(q0),     # ruido de un qubit
    vq.CNOT(q0, q1),
    vq.bit_flip(0.02)(q1),
)

# Solo el simulador de matriz de densidad ejecuta canales
resultado = vq.DensityMatrixSimulator().simulate(circuito)
print(resultado.density_matrix()[0][0])   # población de |00>
```

## Inspeccionar los operadores

```python
import vqubit as vq

operadores = vq.kraus(vq.depolarize(0.1))
print(len(operadores))            # 4: I, X, Y, Z
print(vq.has_unitary(vq.depolarize(0.1)))   # False: un canal no es unitario
```

Una puerta unitaria normal tiene un único operador de Kraus, su propio unitario:

```python
print(len(vq.kraus(vq.X)))        # 1
```

## Siguientes pasos

- Aprende a [optimizar](optimization.md) circuitos.
- Lee sobre [dispositivos](devices.md) y su ejecución local.
