# Devices and running locally

A **device** describes where a circuit may run: which qubits exist and which
connections between them are allowed. vqubit includes local emulations of some
vendor devices, so code written for those APIs runs offline with no credentials.

## Generic devices

```python
import vqubit as vq

line = vq.LineDevice(3)                  # qubits 0, 1, 2 on a line
grid = vq.GridDevice.rect(2, 3)          # a 2x3 grid

print([sorted(device.qubits) for device in (line, grid)])
print(grid.is_adjacent(vq.GridQubit(0, 0), vq.GridQubit(0, 1)))  # True
```

`validate_operation` checks that an operation fits the device:

```python
grid.validate_operation(vq.CNOT(vq.GridQubit(0, 0), vq.GridQubit(0, 1)))
# vq.CNOT(vq.GridQubit(0, 0), vq.GridQubit(5, 5)) would raise ValueError
```

## Vendor emulations

```python
import vqubit as vq

sycamore = vq.Sycamore()          # a 54-qubit grid
print(len(sycamore.qubits))
print(vq.SycamoreGate())          # Google's fSim gate

ion = vq.IonDevice(5)             # an ion trap: all-to-all
```

## The local engine

`Engine` has the same shape as a vendor engine, but runs everything locally on
the simulators.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
engine = vq.Engine(device=vq.LineDevice(2), seed=0)

circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.measure(q0, q1, key="m"))
result = engine.run(circuit, repetitions=1000)
print(result.histogram(key="m"))
```

`Engine` also validates the circuit against its device before running, offers
`run_sweep` for parameter sweeps, and `get_sampler` to obtain a `Sampler`.

## Next steps

- Review [simulation](simulation.md).
- Browse the [API reference](../api/index.md).
