# Diferencias con Cirq

vqubit sigue a Cirq de cerca. Medido contra **cirq-core 1.7.0**, exporta **606 de
sus 607 nombres de nivel superior** —el que falta es `contrib`, un paquete de
integraciones de terceros—, y sus circuitos, puertas, simuladores,
transformadores y dispositivos siguen el mismo diseño, así que el código que ya
conoces suele funcionar tras `import vqubit as vq`.

Aun así, **no** es un reemplazo directo. Esta página recoge los sitios donde
cambia la convención de llamada, para que el código de Cirq que portes falle con
un error en vez de hacer otra cosa en silencio.

## Argumentos renombrados

La diferencia más habitual es la ortografía: vqubit escribe los nombres de los
parámetros completos. Las llamadas que pasan el valor **por posición** no se
ven afectadas; las que usan el nombre deben actualizarse.

| Cirq        | vqubit        | Ejemplo                                     |
| ----------- | ------------- | ------------------------------------------- |
| `p`         | `probability` | `vq.bit_flip(probability=0.1)`              |
| `rads`      | `radians`     | `vq.rx(radians=math.pi / 2)`                |
| `n_qubits`  | `num_qubits`  | `vq.depolarize(0.1, num_qubits=2)`          |
| `val`       | `value`       | `vq.unitary(value=vq.X)`                    |
| `mat`       | `matrix`      | `vq.single_qubit_matrix_to_gates(matrix=u)` |
| `contents`  | `operations`  | `vq.Moment(operations=[vq.H(q0)])`          |
| `tolerance` | `atol`        | `vq.kak_decomposition(matrix, atol=1e-8)`   |

El tamaño de una puerta es un atributo, no un método: `vq.H.num_qubits` vale `1`
donde Cirq tiene `cirq.H.num_qubits()`. Las operaciones sí conservan el método,
así que `vq.H(q0).num_qubits() == 1`.

## Orden de argumentos distinto

Unas pocas llamadas toman los argumentos en otro orden, o reciben primero la
matriz. Lanzan `TypeError` si portas la llamada tal cual.

| Para                | Cirq                                        | vqubit                                             |
| ------------------- | ------------------------------------------- | -------------------------------------------------- |
| Puerta de medida    | `MeasurementGate(num_qubits, key)`          | `MeasurementGate(key, num_qubits)`                 |
| Momento             | `Moment(op1, op2)`                          | `Moment([op1, op2])`                               |
| Síntesis de 2 qubits | `two_qubit_matrix_to_operations(q0, q1, u)` | `two_qubit_matrix_to_operations(u, qubits=(q0, q1))` |
| Síntesis de 3 qubits | `three_qubit_matrix_to_operations(q0, q1, q2, u)` | `three_qubit_matrix_to_operations([q0, q1, q2], u)` |
| Control clásico     | `ClassicallyControlledOperation(op, conditions)` | `ClassicallyControlledOperation(conditions, op)` |

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)

gate = vq.MeasurementGate("m", 2)   # la clave va primero
assert gate.key == "m"

moment = vq.Moment([vq.H(q0), vq.X(q1)])   # un iterable, no *operaciones
assert len(moment) == 2
```

## Las operaciones no admiten `**`, y las puertas planas no son extrapolables

En Cirq se puede elevar cualquier cosa a una potencia: `cirq.X(q) ** 0.5` deja la
puerta a la mitad. En vqubit `**` vive en las clases de puerta *exponenciadas*,
así que `vq.X ** 0.5` funciona, pero `vq.X(q0) ** 0.5` lanza `TypeError`. Eleva
la puerta y luego llámala.

```python
import vqubit as vq

q0 = vq.LineQubit(0)

assert vq.X ** 0.5 == vq.XPowGate(exponent=0.5)
assert (vq.X ** 0.5)(q0) == vq.XPowGate(exponent=0.5)(q0)
```

Lo mismo ocurre con las puertas cuyo equivalente en Cirq está exponenciado:
`vq.CZ`, `vq.CNOT`, `vq.SWAP`, `vq.ISWAP` y `vq.rz(...)` son las puertas planas
(`CZGate`, `CNOTGate`, `SWAPGate`, `ISWAPGate`, `RzGate`) y no definen `**`. Usa
las clases `Pow` para interpolar entre la identidad y la puerta:

```python
import vqubit as vq

assert vq.CZPowGate(exponent=0.5) == vq.CZPowGate() ** 0.5
assert vq.SwapPowGate(exponent=0.5) == vq.SwapPowGate() ** 0.5
```

## Los simuladores solo aceptan `seed`

`vq.Simulator` y `vq.DensityMatrixSimulator` aceptan `seed` y nada más: no hay
`noise=` ni `dtype=`. El ruido viaja **dentro del circuito**, como un canal, y
necesita el simulador de matrices de densidad.

```python
import vqubit as vq

q0 = vq.LineQubit(0)
ruidoso = vq.Circuit(vq.H(q0), vq.depolarize(0.05)(q0), vq.measure(q0, key="m"))

resultado = vq.DensityMatrixSimulator(seed=0).run(ruidoso, repetitions=100)
assert sum(resultado.histogram(key="m").values()) == 100
```

Los vectores de estado y las matrices de densidad son siempre `complex128`. (El
simulador de matrices de densidad de Cirq usa `complex64` por defecto, así que
vqubit es más preciso y ocupa el doble de memoria para el mismo circuito.)

## Mismo nombre, distinto papel

Un puñado de nombres existen en las dos librerías pero significan otra cosa.
Fallan con un error en vez de en silencio, pero hay que adaptar la llamada.

| Nombre | Cirq | vqubit |
| ------ | ---- | ------ |
| `final_state_vector` | recibe un **circuito** | recibe un **resultado de simulación** |
| `final_density_matrix` | recibe un **circuito** | recibe un **resultado de simulación** |
| `read_json` | `read_json(file_or_fn, json_text)` | `read_json(text)` |
| `qasm` | `qasm(val, args, qubits)` | `qasm(value)` |
| `measure` | `measure(q, invert_mask=…)` | `measure(*qubits, key=None)`; la máscara vive en `MeasurementGate` |
| `entanglement_fidelity` | fidelidad de un **canal** | fidelidad entre dos **estados** |
| `xeb_fidelity` | `xeb_fidelity(circuit, bitstrings, …)` | `xeb_fidelity(bitstrings, hilbert_space_dimension)` |
| `estimate_single_qubit_readout_errors` | recibe un **sampler** | recibe **resultados de medida** |
| `prepare_two_qubit_state_using_cz` | `(q0, q1, state)` | `(gate, qubits)` |

Para el estado final, la llamada idiomática en vqubit pasa por el simulador:

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuito = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1))

estado = vq.Simulator().simulate(circuito).state_vector()
assert abs(abs(estado[0]) - 2**-0.5) < 1e-9
assert abs(abs(estado[3]) - 2**-0.5) < 1e-9
```

## Lo que no cambia

Todo lo demás: los qubits (`LineQubit`, `GridQubit`, `NamedQubit`), las clases de
puertas y su aritmética, el comportamiento y los diagramas de `Circuit`/`Moment`,
los protocolos de simulación, `PauliString` y su álgebra, los transformadores y
optimizadores, `to_json`/`read_json` con su formato `cirq_type`, el escritor de
OpenQASM 2.0, las clases de dispositivos y el `Engine`.

## Siguientes pasos

- Conoce las [puertas](quantum-gates.md) y construye [circuitos](circuits.md).
- Lee la [referencia de la API](../api/index.md).
