# Dispositivos y ejecución local

Un **dispositivo** describe dónde puede correr un circuito: qué qubits existen y
qué conexiones están permitidas entre ellos. vqubit incluye emulaciones locales
de algunos dispositivos de fabricantes, de modo que el código escrito para esas
APIs funciona sin conexión y sin credenciales.

## Dispositivos genéricos

```python
import vqubit as vq

linea = vq.LineDevice(3)                 # qubits 0, 1, 2 en línea
rejilla = vq.GridDevice.rect(2, 3)       # rejilla de 2x3

print(sorted(device.qubits) for device in (linea, rejilla))
print(rejilla.is_adjacent(vq.GridQubit(0, 0), vq.GridQubit(0, 1)))  # True
```

`validate_operation` comprueba que una operación cabe en el dispositivo:

```python
rejilla.validate_operation(vq.CNOT(vq.GridQubit(0, 0), vq.GridQubit(0, 1)))
# vq.CNOT(vq.GridQubit(0, 0), vq.GridQubit(5, 5)) lanzaría ValueError
```

## Emulaciones de fabricantes

```python
import vqubit as vq

sycamore = vq.Sycamore()          # rejilla de 54 qubits
print(len(sycamore.qubits))
print(vq.SycamoreGate())          # la puerta fSim de Google

ion = vq.IonDevice(5)             # trampa de iones: todos contra todos
```

## El motor local

`Engine` tiene la misma forma que el motor de un fabricante, pero ejecuta todo
en local sobre los simuladores.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
motor = vq.Engine(device=vq.LineDevice(2), seed=0)

circuito = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.measure(q0, q1, key="m"))
resultado = motor.run(circuito, repetitions=1000)
print(resultado.histogram(key="m"))
```

`Engine` también valida el circuito contra su dispositivo antes de ejecutarlo,
ofrece `run_sweep` para barridos de parámetros y `get_sampler` para obtener un
`Sampler`.

## Siguientes pasos

- Repasa la [simulación](simulation.md).
- Consulta la [referencia de la API](../api/index.md).
