# Example: the Bell state

The **Bell state** is the smallest example of quantum entanglement, and the
"hello world" of quantum computing. This example is complete and you can copy and
run it.

## The circuit

Two qubits, two gates:

```python
import vqubit as vq

# 1. Two qubits on a line.
q0, q1 = vq.LineQubit.range(2)

# 2. A Hadamard gate on q0 and a CNOT controlled by q0.
circuit = vq.Circuit(
    vq.H(q0),
    vq.CNOT(q0, q1),
)

print(circuit)
```

Output:

```
          ┌───┐
0: ───H───@───
          │
1: ───────X───
```

The `H` gate puts `q0` into superposition. Then `CNOT` flips `q1` *only if* `q0`
is `1`. The result is that the two qubits become entangled: their fates are no
longer independent.

## Add a measurement and simulate

```python
# 3. Measure both qubits into an outcome called "result".
circuit.append(vq.measure(q0, q1, key="result"))

# 4. Run the circuit 1000 times.
result = vq.Simulator().run(circuit, repetitions=1000)

# 5. Count how often each outcome came up.
print(result.histogram(key="result"))
```

Typical output:

```
Counter({0: 498, 3: 502})
```

## What those numbers mean

Each measured qubit produces a `0` or a `1`. With two qubits there are four
possible outcomes: `00`, `01`, `10` and `11`. vqubit encodes them as a single
integer, reading `q0` as the most significant bit:

| Outcome | Integer |
| ------- | ------- |
| `00`    | `0`     |
| `01`    | `1`     |
| `10`    | `2`     |
| `11`    | `3`     |

In the output only `0` (`00`) and `3` (`11`) appear, never `1` or `2`. That is,
the two qubits **always agree**. This is entanglement. If you measure `q0` and
get `0`, you know for certain that `q1` will also give `0`, without looking.

## Variations to experiment with

- Replace `CNOT` with `CZ` and note that the outcome set is the same.
- Add an `X(q1)` gate before the measurement and watch the results change.
- Replace `H` with `vq.Ry(q0, radians=math.pi / 8)` and watch the probabilities
  become uneven.
