# Simulation

vqubit ships three local simulators. They all run offline, with no
dependencies.

## State vector

`Simulator` evolves the state vector and is the one you will reach for on any
noise-free circuit.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1))

# Get the final state
result = vq.Simulator().simulate(circuit)
print(result.dirac_notation())   # 0.707|00⟩ + 0.707|11⟩

# Sample measurements
result = vq.Simulator(seed=0).run(
    vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.measure(q0, q1, key="m")),
    repetitions=1000,
)
print(result.histogram(key="m"))  # Counter({0: ..., 3: ...})
```

`simulate` returns a `SimulationResult` with `final_state` and `state_vector()`;
`run` returns a `Result` with `measurements` and `histogram`.

!!! note "If you have `numpy` installed"
    `run` draws all the `repetitions` at once with `numpy`; without it draws
    with the running totals (`bisect`) over the same distribution. Both paths
    simulate the circuit **once** when the measurements are terminal, so 2000
    shots of 16 qubits take about 0.02 s either way. `numpy` stays optional.

## Density matrix

`DensityMatrixSimulator` evolves a density matrix, so it can simulate **noisy**
circuits.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuit = vq.Circuit(vq.H(q0), vq.depolarize(0.1)(q0), vq.CNOT(q0, q1))

result = vq.DensityMatrixSimulator().simulate(circuit)
rho = result.density_matrix()
print(rho[0][3])  # the coherence, reduced by the noise
```

## Expectation values

`simulate_expectation_values` gives the exact expectation value of one or more
Pauli observables on the final state, and `sample_expectation_values` estimates
them by measuring (which the Clifford simulator can do too).

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1))
observables = [vq.PauliString({q0: vq.Z, q1: vq.Z}), vq.PauliString({q0: vq.Z})]

print(vq.Simulator().simulate_expectation_values(circuit, observables))  # [1, 0]
```

## Clifford (stabilizers)

For Clifford circuits (H, S, CNOT, CZ, SWAP and Paulis), `CliffordSimulator`
tracks a stabilizer *tableau* instead of the state vector, so it scales far
better.

```python
import vqubit as vq

q0, q1, q2 = vq.LineQubit.range(3)
circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.CNOT(q1, q2))

state = vq.CliffordSimulator().simulate(circuit).final_state
print([str(s) for s in state.stabilizers()])
# ['X(0)*X(1)*X(2)', 'Z(0)*Z(1)', 'Z(1)*Z(2)']
```

## Next steps

- Add [noise](noise.md) to your circuits.
- Optimize your circuits in the [optimization](optimization.md) chapter.
