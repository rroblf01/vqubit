# Ejemplo: búsqueda de Grover

El algoritmo de **Grover** encuentra un elemento marcado en una lista no
ordenada con ~√N consultas, donde una búsqueda clásica necesitaría N. Con dos
qubits la lista tiene cuatro elementos y **una sola iteración** basta para
encontrar el marcado con certeza.

## El circuito

Marcamos `|11>` y aplicamos una iteración: el oráculo (un `CZ`, que invierte el
signo de `|11>`) y el difusor (que refleja sobre la media).

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
grover = vq.Circuit(
    vq.H(q0), vq.H(q1),                      # 1. superposición uniforme
    vq.CZ(q0, q1),                            # 2. oráculo: marca |11>
    vq.H(q0), vq.H(q1),                       # 3. difusor: H, X, CZ, X, H
    vq.X(q0), vq.X(q1),
    vq.CZ(q0, q1),
    vq.X(q0), vq.X(q1),
    vq.H(q0), vq.H(q1),
)
```

La amplitud de cada estado después del circuito:

```python
import numpy as np

state = vq.Simulator().simulate(grover).state_vector()
print(np.round(np.abs(state), 3))
```

```text
[0. 0. 0. 1.]
```

Toda la amplitud está en `|11>` (el índice 3), así que medir da `11` siempre:

```python
grover.append(vq.measure(q0, q1, key="m"))
print(vq.Simulator(seed=1).run(grover, repetitions=1000).histogram(key="m"))
```

```text
Counter({3: 1000})
```

## Por qué funciona

El oráculo invierte el signo del estado marcado y el difusor lo refleja sobre la
media de las amplitudes. La combinación es una rotación de la esfera de Bloch
colectiva que, con una sola iteración y cuatro estados, lleva toda la amplitud
al marcado. Con más qubits hacen falta ~π√N/4 iteraciones, y la probabilidad de
éxito deja de ser exacta.

## Qué has aprendido

- Un oráculo es una fase condicional: `CZ` ya marca `|11>` sin ninguna puerta
  extra.
- El difusor se escribe con Hadamards y `X`s alrededor del mismo `CZ`.
- El histograma de 1000 ejecuciones sale entero en el estado marcado: la
  simulación es exacta y el circuito determinista.

Para el catálogo de puertas, mira [Puertas cuánticas](../guide/quantum-gates.md);
para el simulador de estabilizadores (que también corre este circuito), el
[ejemplo del código de 3 qubits](stabilizers.md).
