# vqubit

**Un framework de computación cuántica para Python, sin dependencias, escrito en [V](https://vlang.io).**

vqubit te permite construir circuitos cuánticos y ejecutarlos en simuladores
locales **sin ninguna dependencia externa**. Todo lo que hace un uso intensivo
del cálculo —álgebra lineal, matrices de puertas y simulación— está escrito en V
y compilado a una extensión nativa; la API pública es Python puro sobre ella.

La API sigue de cerca a [Cirq](https://quantumai.google/cirq), así que si ya
conoces Cirq, ya conoces vqubit.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuito = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.measure(q0, q1, key="result"))
print(circuito)

resultado = vq.Simulator().run(circuito, repetitions=100)
print(resultado.histogram(key="result"))
```

## Por qué vqubit

- **Sin dependencias.** `pip install vqubit` no instala nada más. Todo funciona
  en local, en tu máquina.
- **Velocidad nativa.** El álgebra lineal, las matrices de puertas y la
  simulación están compiladas desde V, no interpretadas.
- **Una rueda para todas las versiones de Python.** Una única rueda *abi3* cubre
  CPython 3.11–3.15 por plataforma.
- **Multiplataforma.** Linux (x86-64 y ARM), macOS (ARM) y Windows (x86-64).
- **API compatible con Cirq.** Construye, transforma y simula circuitos tal como
  documenta Cirq; mira [qué cambia](guide/cirq-differences.md).

## ¿Nuevo en computación cuántica?

No hace falta que sepas física para empezar. Si nunca has oído hablar de un
qubit, empieza por [¿Qué es un qubit?](guide/what-is-a-qubit.md) y sigue la guía
en orden. Cada concepto se explica desde cero, con ejemplos que puedes ejecutar.

## Siguientes pasos

1. [Instala vqubit](installation.md).
2. Sigue el [inicio rápido](quickstart.md) para ejecutar tu primer circuito.
3. Aprende los [conceptos básicos](guide/what-is-a-qubit.md).
4. Consulta la [referencia de la API](api/index.md).
