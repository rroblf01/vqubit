# Inicio rápido

Esta página te lleva de cero a un circuito cuántico ejecutándose en un minuto.

## 1. Crea dos qubits

Un **qubit** es la unidad básica de información cuántica, el equivalente
cuántico de un bit. En vqubit describes dónde vive cada qubit usando un tipo de
"identificador" (`Qid`). `LineQubit.range(2)` crea dos qubits en una línea:

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
```

## 2. Construye un circuito

Un **circuito** es una secuencia de operaciones ordenadas en el tiempo. Aquí
preparamos un **estado de Bell**, un estado entrelazado de dos qubits:

```python
circuito = vq.Circuit(
    vq.H(q0),        # puerta Hadamard sobre q0: crea superposición
    vq.CNOT(q0, q1), # puerta CNOT: entrelaza q0 y q1
)
print(circuito)
```

## 3. Añade una medida

Medir convierte el estado cuántico en bits clásicos que puedes leer:

```python
circuito.append(vq.measure(q0, q1))
```

## 4. Simúlalo

vqubit trae simuladores locales. El más sencillo es el `Simulator`, que
evoluciona el vector de estado:

```python
resultado = vq.Simulator().run(circuito, repetitions=1000)
print(resultado.histogram(key="result"))
```

Obtendrás algo parecido a:

```
Counter({0: 503, 3: 497})
```

Las claves `0` y `3` son los resultados posibles de medir dos qubits: `00` y
`11`. Solo aparecen esos dos, nunca `01` ni `10`, porque los qubits están
entrelazados. ¡Eso es el entrelazamiento cuántico en acción!

## Siguientes pasos

- Aprende [qué es un qubit](guide/what-is-a-qubit.md).
- Conoce las [puertas cuánticas](guide/quantum-gates.md).
- Construye [circuitos](guide/circuits.md) más grandes.
