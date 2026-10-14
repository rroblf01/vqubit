# ¿Qué es un qubit?

Si nunca has estudiado computación cuántica, esta página es tu punto de partida.
No necesitas saber física: solo iremos despacio.

## Bits, pero cuánticos

Un ordenador clásico guarda información en **bits**. Un bit solo puede estar en
uno de dos estados: `0` o `1`, como un interruptor apagado o encendido.

Un **qubit** (bit cuántico) es la versión cuántica de esa idea, pero con dos
diferencias importantes:

1. **Superposición.** Mientras nadie lo mide, un qubit puede estar *a la vez* en
   `0` y en `1`, con unos pesos (amplitudes) para cada opción.
2. **Entrelazamiento.** Varios qubits pueden quedar tan ligados que el estado de
   uno depende del otro, aunque estén separados.

Aunque parezca magia, es simplemente cómo se comporta la materia a escalas muy
pequeñas. Las matemáticas que lo describen son el álgebra lineal: un qubit es un
vector de dos números complejos.

## Cómo se ve en vqubit

No necesitas manejar esos vectores a mano. En vqubit, lo primero es decidir
**dónde** vive cada qubit. Ese "dónde" se llama identificador cuántico o `Qid`:

```python
import vqubit as vq

q0 = vq.LineQubit(0)      # el qubit 0 de una línea
q1 = vq.LineQubit(1)      # el qubit 1 de la misma línea
```

Una línea (`LineQubit`) es solo un modelo sencillo: los qubits se numeran `0, 1,
2, ...`. También existen rejillas (`GridQubit`) para imitar la disposición
física de muchos procesadores reales, y qubits con nombre (`NamedQubit`).

## Superposición con una sola puerta

La puerta **Hadamard** (`H`) pone un qubit en superposición: a partes iguales
entre `0` y `1`.

```python
circuito = vq.Circuit(vq.H(q0))
estado = vq.Simulator().simulate(circuito, qubit_order=[q0])
print(estado)
```

Si mides ese qubit muchas veces, verás `0` aproximadamente la mitad de las veces
y `1` la otra mitad. La superposición no significa que "esté en los dos a la
vez" en el sentido cotidiano: significa que, hasta que lo mides, sus
probabilidades están repartidas.

## Medición

Medir un qubit lo obliga a "decidirse" por `0` o por `1`. El resultado es
aleatorio, pero con probabilidades que dependen del estado. En vqubit:

```python
resultado = vq.Simulator().run(vq.Circuit(H(q0), vq.measure(q0)), repetitions=1000)
print(resultado.histogram(key="result"))
```

## Siguientes pasos

Ahora que sabes qué es un qubit, aprende qué le puedes hacer:
[las puertas cuánticas](quantum-gates.md).
