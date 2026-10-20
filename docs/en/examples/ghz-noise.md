# Example: GHZ with noise

A three-qubit **GHZ** state is the three-qubit version of the Bell state: a
superposition of `|000>` and `|111>`. On a real machine every gate adds a little
noise, and this example shows how to model it and what it does to the state.

## The perfect state

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

Only `000` (0) and `111` (7) come out: the three qubits always agree, which is
what "entangled" means here.

## Adding noise

A **noise model** turns a clean circuit into one with channels. The simplest one
applies a noise gate to every qubit in every moment:

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

The `D` is the depolarizing channel. A channel is not unitary, so it has to be
simulated with the density matrix:

```python
result = vq.DensityMatrixSimulator(seed=7).run(noisy, repetitions=2000)
print(result.histogram(key="m"))
```

```text
Counter({0: 803, 7: 770, 4: 107, 3: 105, 1: 88, 6: 63, 2: 34, 5: 30})
```

Results with an odd number of ones (`1`, `2`, `4`…) now appear: the noise has
broken the correlation between the three qubits.

## Noise from the device's properties

A real device declares an error per gate kind rather than one number:

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

The spread is much smaller than with the flat 5 %, because one-qubit gates only
fail 1 % of the time and two-qubit gates 3 %.

## What you have learned

- `NoiseModel.noisy_moments` returns the moments with the noise already in
  them; the circuit is built with `vq.Circuit(...)`.
- Channels need `DensityMatrixSimulator`: the state-vector `Simulator` only
  applies unitary gates and raises `TypeError` when handed a channel.
- `histogram(key=...)` counts the outcomes per measurement key, reading the
  register as a big-endian integer.

For the theory behind channels see [Noise](../guide/noise.md); for a device's
properties see [Devices](../guide/devices.md).
