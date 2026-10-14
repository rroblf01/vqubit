# Quickstart

This page takes you from nothing to a running quantum circuit in a minute.

## 1. Create two qubits

A **qubit** is the basic unit of quantum information, the quantum analogue of a
bit. In vqubit you describe where each qubit lives with an *identifier* type
(`Qid`). `LineQubit.range(2)` creates two qubits on a line:

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
```

## 2. Build a circuit

A **circuit** is a time-ordered sequence of operations. Here we prepare a
**Bell state**, an entangled state of two qubits:

```python
circuit = vq.Circuit(
    vq.H(q0),        # Hadamard on q0: creates superposition
    vq.CNOT(q0, q1), # CNOT: entangles q0 and q1
)
print(circuit)
```

## 3. Add a measurement

Measurement turns the quantum state into classical bits you can read:

```python
circuit.append(vq.measure(q0, q1))
```

## 4. Simulate it

vqubit ships local simulators. The simplest is `Simulator`, which evolves the
state vector:

```python
result = vq.Simulator().run(circuit, repetitions=1000)
print(result.histogram(key="result"))
```

You will get something like:

```
Counter({0: 503, 3: 497})
```

The keys `0` and `3` are the possible outcomes of measuring two qubits: `00` and
`11`. Only those two appear, never `01` or `10`, because the qubits are
entangled. That is quantum entanglement in action.

## Next steps

- Learn [what a qubit is](guide/what-is-a-qubit.md).
- Meet the [quantum gates](guide/quantum-gates.md).
- Build larger [circuits](guide/circuits.md).
