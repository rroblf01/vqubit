# Circuits

A **circuit** is the central structure in vqubit: an ordered sequence of quantum
operations. If gates are the words, the circuit is the sentence.

## A circuit is a list of moments

vqubit organises operations into **moments** (`Moment`). A moment groups
operations that happen *at the same time* (on different qubits). This reflects
how hardware works: simultaneous operations must use different qubits.

When you build a circuit from gates, vqubit decides automatically which moment
each one goes into, respecting order and avoiding conflicts:

```python
import vqubit as vq

q0, q1, q2 = vq.LineQubit.range(3)

circuit = vq.Circuit(
    vq.H(q0),
    vq.CNOT(q0, q1),
    vq.X(q2),        # can share a moment with H(q0)
)
print(circuit)
```

It is drawn as text, one line per qubit:

```
0: ───H───@───
          │
1: ───────X───
2: ───X───────
```

Read the diagram left to right: each column is a moment.

## Building step by step

You can create an empty circuit and add operations later:

```python
circuit = vq.Circuit()
circuit.append(vq.H(q0))
circuit.append(vq.measure(q0, key="m"))
```

Or insert at a specific moment:

```python
circuit.insert(0, vq.X(q1))   # insert an operation early
```

## Iterating over a circuit

A circuit can be traversed by moments or by operations:

```python
for moment in circuit:
    print(moment)

for operation in circuit.all_operations():
    print(operation)
```

## Running the circuit

A circuit is a description; to get results you have to **simulate** it:

```python
result = vq.Simulator().run(circuit, repetitions=1000)
print(result.histogram(key="m"))
```

## Next steps

See a complete, commented example: the [Bell state](../examples/bell-state.md).
