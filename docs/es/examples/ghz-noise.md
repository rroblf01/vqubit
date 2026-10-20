# Ejemplo: GHZ con ruido

Un estado **GHZ** de tres qubits es la versión de tres qubits del estado de
Bell: una superposición de `|000>` y `|111>`. En un ordenador real, cada puerta
introduce un poco de ruido, y este ejemplo muestra cómo modelarlo y qué le hace
al estado.

## El estado perfecto

```python
import vqubit as vq

q0, q1, q2 = vq.LineQubit.range(3)
ghz = vq.Circuit(
    vq.H(q0),
    vq.CNOT(q0, q1),
    vq.CNOT(q1, q2),
    vq.measure(q0, q1, q2, key="m"),
)
print(vq.Simulator(seed=7).run(ghz, repetitions=2000).histogram(key="m"))
```

```text
Counter({0: 1013, 7: 987})
```

Solo salen `000` (0) y `111` (7): los tres qubits están de acuerdo, que es lo
que significa "entrelazados".

## Añadir ruido

Un **modelo de ruido** convierte un circuito limpio en uno con canales. El más
simple aplica una puerta de ruido a cada qubit en cada momento:

```python
model = vq.ConstantQubitNoiseModel(vq.depolarize(0.05))
qubits = sorted(ghz.all_qubits())
noisy = vq.Circuit(model.noisy_moments(ghz.moments, qubits))
print(noisy)
```

```text
0: ───H───D───@───D───────────M───D───
              │   │           │   │
1: ───────────X───D───@───D───M───D───
                      │   │   │   │
2: ───────────────────X───D───M───D───
```

La `D` es el canal despolarizante. Como un canal no es unitario, hay que
simularlo con la matriz de densidad:

```python
result = vq.DensityMatrixSimulator(seed=7).run(noisy, repetitions=2000)
print(result.histogram(key="m"))
```

```text
Counter({0: 803, 7: 770, 4: 107, 3: 105, 1: 88, 6: 63, 2: 34, 5: 30})
```

Ahora aparecen resultados con un número impar de unos (`1`, `2`, `4`…): el
ruido ha roto la correlación entre los tres qubits.

## Ruido según las propiedades del dispositivo

Un dispositivo real declara su error por tipo de puerta, no un número único:

```python
properties = vq.SuperconductingQubitsNoiseProperties(
    gate_error_1q=0.01,
    gate_error_2q=0.03,
    readout_error_rate=0.02,
)
device_model = vq.NoiseModelFromNoiseProperties(properties)
realistic = vq.Circuit(device_model.noisy_moments(ghz.moments, qubits))
print(vq.DensityMatrixSimulator(seed=7).run(realistic, repetitions=2000).histogram(key="m"))
```

```text
Counter({0: 995, 7: 924, 4: 30, 3: 26, 5: 8, 6: 8, 2: 6, 1: 3})
```

El reparto es mucho más limpio que con el 5 % constante, porque las puertas de
un qubit solo fallan el 1 % y las de dos el 3 %.

## Qué has aprendido

- `NoiseModel.noisy_moments` devuelve los momentos ya con el ruido; el circuito
  se construye con `vq.Circuit(...)`.
- Los canales necesitan `DensityMatrixSimulator`: el `Simulator` de vector de
  estado solo aplica puertas unitarias y lanza `TypeError` si le das un canal.
- `histogram(key=...)` cuenta los resultados por clave de medida, con el
  registro leído como un entero big-endian.

Para la teoría detrás de los canales, mira [Ruido](../guide/noise.md); para las
propiedades de un dispositivo, [Dispositivos](../guide/devices.md).
