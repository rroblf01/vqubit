# Circuitos

Un **circuito** es la estructura central de vqubit: una secuencia ordenada de
operaciones cuánticas. Si las puertas son las palabras, el circuito es la frase.

## Un circuito es una lista de momentos

vqubit organiza las operaciones en **momentos** (`Moment`). Un momento agrupa
las operaciones que ocurren *a la vez* (sobre qubits distintos). Esto refleja
cómo funciona el hardware: las operaciones simultáneas deben usar qubits
diferentes.

Cuando creas un circuito a partir de puertas, vqubit decide automáticamente en
qué momento va cada una, respetando el orden y evitando conflictos:

```python
import vqubit as vq

q0, q1, q2 = vq.LineQubit.range(3)

circuito = vq.Circuit(
    vq.H(q0),
    vq.CNOT(q0, q1),
    vq.X(q2),        # puede ir en el mismo momento que H(q0)
)
print(circuito)
```

Se dibuja como texto, con una línea por qubit:

```
          ┌───┐
0: ───H───@───────
          │
1: ───────X───────

2: ───X───────────
```

Lee el diagrama de izquierda a derecha: cada columna es un momento.

## Construir paso a paso

Puedes crear un circuito vacío y añadir operaciones después:

```python
circuito = vq.Circuit()
circuito.append(vq.H(q0))
circuito.append(vq.measure(q0, key="m"))
```

O insertar en un momento concreto:

```python
circuito.insert(0, vq.X(q1))   # inserta una operación temprano
```

## Iterar por el circuito

Un circuito se puede recorrer por momentos o por operaciones:

```python
for momento in circuito:
    print(momento)

for operacion in circuito.all_operations():
    print(operacion)
```

## Ejecutar el circuito

Un circuito es una descripción; para obtener resultados hay que **simularlo**:

```python
resultado = vq.Simulator().run(circuito, repetitions=1000)
print(resultado.histogram(key="m"))
```

## Siguientes pasos

Ve un ejemplo completo y comentado: el [estado de Bell](../examples/bell-state.md).
