# Ejemplo: el estado de Bell

El **estado de Bell** es el ejemplo más pequeño de entrelazamiento cuántico, y
el "hola mundo" de la computación cuántica. Este ejemplo está completo y puedes
copiarlo y ejecutarlo.

## El circuito

Dos qubits, dos puertas:

```python
import vqubit as vq

# 1. Dos qubits en una línea.
q0, q1 = vq.LineQubit.range(2)

# 2. Una puerta Hadamard sobre q0 y una CNOT controlada por q0.
circuito = vq.Circuit(
    vq.H(q0),
    vq.CNOT(q0, q1),
)

print(circuito)
```

Salida:

```
0: ───H───@───
          │
1: ───────X───
```

La puerta `H` pone a `q0` en superposición. Después, `CNOT` invierte `q1` *solo
si* `q0` vale `1`. El resultado es que los dos qubits quedan entrelazados: sus
destinos ya no son independientes.

## Añadir una medida y simular

```python
# 3. Medimos ambos qubits en un resultado llamado "result".
circuito.append(vq.measure(q0, q1, key="result"))

# 4. Ejecutamos el circuito 1000 veces.
resultado = vq.Simulator().run(circuito, repetitions=1000)

# 5. Contamos cuántas veces salió cada resultado.
print(resultado.histogram(key="result"))
```

Salida típica:

```
Counter({0: 498, 3: 502})
```

## Qué significan esos números

Cada qubit medido produce un `0` o un `1`. Con dos qubits hay cuatro resultados
posibles: `00`, `01`, `10` y `11`. vqubit los codifica como un entero, leyendo
`q0` como el bit más significativo:

| Resultado | Entero |
| --------- | ------ |
| `00`      | `0`    |
| `01`      | `1`    |
| `10`      | `2`    |
| `11`      | `3`    |

En la salida solo aparecen `0` (`00`) y `3` (`11`), nunca `1` ni `2`. Es decir:
los dos qubits **siempre coinciden**. Eso es el entrelazamiento. Si mides `q0` y
obtienes `0`, sabes con certeza que `q1` también dará `0`, sin mirarlo.

## Variaciones para experimentar

- Cambia `CNOT` por `CZ` y observa que el resultado es el mismo conjunto.
- Añade una puerta `X(q1)` antes de la medida y mira cómo cambian los resultados.
- Sustituye `H` por `vq.Ry(math.pi / 8)(q0)` y observa cómo se
  desequilibran las probabilidades.
