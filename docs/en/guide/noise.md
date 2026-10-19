# Noise

A **channel** is a gate that is not unitary: it describes noise and takes a pure
state to a mixed one. vqubit defines channels with their **Kraus operators**,
the matrices $K_i$ such that the state evolves as
$\rho \to \sum_i K_i \rho K_i^\dagger$.

## The channels that ship

| Channel | Constructor | What it models |
| ------- | ----------- | -------------- |
| Bit flip | `vq.bit_flip(probability)` | Flips the qubit with probability `probability`. |
| Phase flip | `vq.phase_flip(probability)` | Applies `Z` with probability `probability`. |
| Bit-phase flip | `vq.bit_phase_flip(probability)` | Applies `Y` with probability `probability`. |
| Depolarizing | `vq.depolarize(probability, num_qubits=1)` | Applies a random Pauli with probability `probability`. |
| Amplitude damping | `vq.amplitude_damp(gamma)` | Energy loss (`|1>` → `|0>`). |
| Phase damping | `vq.phase_damp(gamma)` | Loss of coherence. |
| Reset | `vq.reset(qubit)` | Drives a qubit to `|0>`. |

The argument names are the ones the functions accept, so you can pass them by
keyword; a probability outside `[0, 1]` raises `ValueError`.

## Using them in a circuit

Channels are placed like any other operation:

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuit = vq.Circuit(
    vq.H(q0),
    vq.depolarize(0.05)(q0),     # single-qubit noise
    vq.CNOT(q0, q1),
    vq.bit_flip(0.02)(q1),
)

# Only the density-matrix simulator runs channels
result = vq.DensityMatrixSimulator().simulate(circuit)
print(result.density_matrix()[0][0])   # population of |00>
```

## Inspecting the operators

```python
import vqubit as vq

operators = vq.kraus(vq.depolarize(0.1))
print(len(operators))            # 4: I, X, Y, Z
print(vq.has_unitary(vq.depolarize(0.1)))   # False: a channel is not unitary
```

An ordinary unitary gate has exactly one Kraus operator, its own unitary:

```python
print(len(vq.kraus(vq.X)))        # 1
```

## Next steps

- Learn to [optimize](optimization.md) circuits.
- Read about [devices](devices.md) and running locally.
