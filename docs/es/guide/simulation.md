# Simulación

vqubit trae tres simuladores locales. Todos funcionan sin conexión y sin
dependencias.

## Vector de estado

El `Simulator` evoluciona el vector de estado y es el que usarás casi siempre
para circuitos sin ruido.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuito = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1))

# Obtener el estado final
resultado = vq.Simulator().simulate(circuito)
print(resultado.dirac_notation())   # 0.707|00⟩ + 0.707|11⟩

# Muestrear medidas
resultado = vq.Simulator(seed=0).run(
    vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.measure(q0, q1, key="m")),
    repetitions=1000,
)
print(resultado.histogram(key="m"))  # Counter({0: ..., 3: ...})
```

`simulate` devuelve un `SimulationResult` con `final_state` y `state_vector()`;
`run` devuelve un `Result` con `measurements` y `histogram`.

## Matriz de densidad

El `DensityMatrixSimulator` evoluciona una matriz de densidad, así que puede
simular circuitos **con ruido**.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuito = vq.Circuit(vq.H(q0), vq.depolarize(0.1)(q0), vq.CNOT(q0, q1))

resultado = vq.DensityMatrixSimulator().simulate(circuito)
rho = resultado.density_matrix()
print(rho[0][3])  # la coherencia, reducida por el ruido
```

## Clifford (estabilizadores)

Para circuitos de Clifford (H, S, CNOT, CZ, SWAP y Paulis) el
`CliffordSimulator` sigue un *tableau* de estabilizadores en lugar del vector
de estado, así que escala mucho mejor.

```python
import vqubit as vq

q0, q1, q2 = vq.LineQubit.range(3)
circuito = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.CNOT(q1, q2))

estado = vq.CliffordSimulator().simulate(circuito).final_state
print(estado.stabilizers())
# ['X(0)*X(1)*X(2)', 'Z(0)*Z(1)', 'Z(1)*Z(2)']
```

## Siguientes pasos

- Añade [ruido](noise.md) a tus circuitos.
- Optimiza tus circuitos en el capítulo de [optimización](optimization.md).
