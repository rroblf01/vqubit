# What is a qubit?

If you have never studied quantum computing, this page is your starting point.
You do not need any physics: we will go slowly.

## Bits, but quantum

A classical computer stores information in **bits**. A bit can be in one of two
states: `0` or `1`, like a switch that is off or on.

A **qubit** (quantum bit) is the quantum version of that idea, with two
important differences:

1. **Superposition.** Until it is measured, a qubit can be *both* `0` and `1` at
   once, with weights (amplitudes) for each option.
2. **Entanglement.** Several qubits can become so tightly linked that the state
   of one depends on another, even when they are far apart.

It sounds like magic, but it is simply how matter behaves at very small scales.
The mathematics behind it is linear algebra: a qubit is a vector of two complex
numbers.

## What it looks like in vqubit

You do not have to handle those vectors by hand. In vqubit the first thing is to
decide **where** each qubit lives. That "where" is called a quantum identifier,
or `Qid`:

```python
import vqubit as vq

q0 = vq.LineQubit(0)      # qubit 0 of a line
q1 = vq.LineQubit(1)      # qubit 1 of the same line
```

A line (`LineQubit`) is just a simple model: qubits are numbered `0, 1, 2, ...`.
There are also grids (`GridQubit`), which imitate the physical layout of many
real processors, and named qubits (`NamedQubit`).

## Superposition with a single gate

The **Hadamard** gate (`H`) puts a qubit into superposition: an equal split
between `0` and `1`.

```python
circuit = vq.Circuit(vq.H(q0))
state = vq.Simulator().simulate(circuit, qubit_order=[q0])
print(state)
```

If you measure that qubit many times, you will see `0` roughly half the time and
`1` the other half. Superposition does not mean it is "both at once" in the
everyday sense: it means that, until you measure it, its probabilities are
spread out.

## Measurement

Measuring a qubit forces it to "decide" between `0` and `1`. The outcome is
random, with probabilities that depend on the state. In vqubit:

```python
result = vq.Simulator().run(vq.Circuit(vq.H(q0), vq.measure(q0, key="result")), repetitions=1000)
print(result.histogram(key="result"))
```

## Next steps

Now that you know what a qubit is, learn what you can do to it:
[the quantum gates](quantum-gates.md).
