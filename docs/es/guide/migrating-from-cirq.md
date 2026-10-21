# Migrar desde Cirq

Si ya programas con Cirq, casi todo lo que sabes sigue valiendo: vqubit copia
sus nombres, sus firmas y su forma de comportarse. Esta guía es la lista de
comprobación para portar un programa real, en el orden en que aparecen los
problemas.

!!! tip "La versión corta"
    Cambia `import cirq` por `import vqubit as vq` y ejecuta. Después:
    pon los nombres de argumento completos (`probability=`, `radians=`,
    `num_qubits=`), usa `num_qubits` como **atributo** en las puertas, y revisa
    las cuatro diferencias de comportamiento de la última sección.

## 1. El import

```python
# antes
import cirq

# después
import vqubit as vq
```

vqubit expone **606 de los 607 nombres de nivel superior** de cirq-core 1.7.0
(el que falta es `contrib`). Eso incluye los submódulos que imitan a Cirq:

```python
from vqubit import ops, sim, linalg, protocols, transformers, work
```

Si en Cirq escribías `cirq.ops.X`, en vqubit puedes escribir `vq.ops.X`; el
nombre plano `vq.X` también existe y es el mismo objeto.

## 2. Los argumentos cambian de nombre

vqubit escribe los parámetros **completos**. Las llamadas posicionales no se
ven afectadas; las que usan el nombre, sí:

| Cirq                     | vqubit                          |
| ------------------------ | ------------------------------- |
| `p=0.1`                  | `probability=0.1`               |
| `rads=math.pi / 2`       | `radians=math.pi / 2`           |
| `n_qubits=2`             | `num_qubits=2`                  |
| `val=vq.X`               | `value=vq.X`                    |
| `mat=u`                  | `matrix=u`                      |
| `contents=[...]`         | `operations=[...]`              |
| `tolerance=1e-8`         | `atol=1e-8`                     |

Un error de nombre lanza `TypeError`, no hace otra cosa en silencio.

## 3. `num_qubits` es un atributo en las puertas

```python
# Cirq
cirq.H.num_qubits()      # método
cirq.H(q0).num_qubits()  # método

# vqubit
vq.H.num_qubits          # atributo
vq.H(q0).num_qubits()    # las operaciones sí conservan el método
```

La regla: **puerta → atributo**, **operación → método**. Es el cambio más
habitual al portar bucles que recorren operaciones.

## 4. Medir y muestrear

`sample` devuelve un `vqubit.Result`, no un `pandas.DataFrame`. Si solo quieres
las cuentas, usa `histogram`:

```python
result = vq.Simulator().run(circuito, repetitions=1000)
print(result.histogram(key="m"))          # Counter({0: 512, 1: 488})
medidas = result.measurements["m"]        # lista de filas
```

`run_sweep`, `simulate_sweep`, `run_batch`, `simulate_moment_steps` y
`simulate_expectation_values` aceptan los mismos argumentos que en Cirq, con el
resolver de parámetros en **segunda** posición.

## 5. Los parámetros no son símbolos de SymPy

```python
# Cirq
import sympy
t = sympy.Symbol("t")

# vqubit
t = vq.Parameter("t")
```

La aritmética (`t / 2`, `2 * t + 1`, `-t`), `resolve_parameters` y las puertas
simbólicas funcionan igual, pero la igualdad es la de vqubit (dos expresiones
son iguales si lo son estructuralmente). Si necesitas SymPy de verdad,
`vqubit.interop` documenta cómo convertir: ver
[Serialización](serialization.md).

## 6. Las puertas planas son las puertas de potencia

```python
vq.CZ is vq.CZPowGate()        # True, como en Cirq
vq.CZ ** 0.5                   # CZPowGate(exponent=0.5)
vq.CNOT ** -1 == vq.CNOT       # True: la igualdad es periódica
vq.ISWAP ** -1 != vq.ISWAP     # ISWAP sí distingue -1
```

Los alias históricos (`vq.CZGate`, `vq.CNOTGate`, `vq.SWAPGate`, `vq.CCXGate`,
`vq.ISWAPGate`) siguen existiendo, así que un `isinstance` no se rompe.

## 7. Lo que no está

- **`cirq.contrib`**: no existe (su importador de QASM depende de `ply`). En su
  lugar, vqubit trae `vq.circuit_from_qasm`, sin dependencias.
- **El árbol `EigenGate`**: las potencias de puertas multi-qubit heredan de
  `_PowGate`, no de `EigenGate`.
- **`cirq.SimulationState` con qudits**: vqubit rechaza explícitamente las
  formas que no son de qubits en los helpers de medida.

Todo lo demás —qubits, puertas, circuitos, simuladores, transformadores,
dispositivos, `cirq.work`— está, con el mismo diseño.

## 8. Las cuatro diferencias de comportamiento

Están explicadas, con ejemplos y con su criterio de cierre, en
[Diferencias con Cirq](cirq-differences.md). Resumen:

1. **Diagramas de texto**: las puertas Clifford de dos qubits se dibujan con su
   nombre (`CliffordGate(2q)`), y los exponentes salen en decimal (`^0.5`).
2. **Igualdad de clases**: exige la misma clase concreta (Cirq compara contra la
   clase decorada). Solo se nota con `SqrtISwapGate`, que no existe en Cirq.
3. **`sample`**: devuelve `Result`, no un `DataFrame`.
4. **`Parameter`**: no es un `sympy.Symbol`.

## 9. Un programa portado de principio a fin

Este es el mismo programa en las dos bibliotecas: prepara un estado GHZ de tres
qubits, lo mide y estima el valor esperado de `ZZZ`.

```python
# Cirq
import cirq

q0, q1, q2 = cirq.LineQubit.range(3)
circuit = cirq.Circuit(cirq.H(q0), cirq.CNOT(q0, q1), cirq.CNOT(q1, q2))
circuit.append(cirq.measure(q0, q1, q2, key="m"))
result = cirq.Simulator().run(circuit, repetitions=100)
print(result.histogram(key="m"))
observable = cirq.X(q0) * cirq.X(q1) * cirq.X(q2)
print(cirq.Simulator().simulate_expectation_values(circuit, observables=[observable]))
```

```python
# vqubit
import vqubit as vq

q0, q1, q2 = vq.LineQubit.range(3)
circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.CNOT(q1, q2))
circuit.append(vq.measure(q0, q1, q2, key="m"))
result = vq.Simulator().run(circuit, repetitions=100)
print(result.histogram(key="m"))
observable = vq.PauliString({q0: vq.Pauli.X, q1: vq.Pauli.X, q2: vq.Pauli.X})
print(vq.Simulator().simulate_expectation_values(circuit, observables=[observable]))
```

Los dos imprimen un histograma con solo `0` y `7` (los dos estados del GHZ) y
un valor esperado de `1.0`: `X ⊗ X ⊗ X` es el estabilizador del estado, mientras
que `Z ⊗ Z ⊗ Z` daría `0` (los dos términos del GHZ tienen signo opuesto).

## 10. Comprobar la paridad

Si dudas de una llamada concreta, el repositorio incluye dos arneses que
comparan vqubit con cirq-core 1.7.0:

```bash
python benchmarks/parity.py    # 1867 comprobaciones de API
python benchmarks/verify.py    # 650 comprobaciones numéricas
```

`benchmarks/parity_tolerated.toml` es la lista cerrada de divergencias
deliberadas: si algo no está ahí, es un fallo.
