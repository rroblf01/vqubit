# vqubit

**A dependency-free quantum computing framework for Python, written in [V](https://vlang.io).**

vqubit lets you build quantum circuits and run them on local simulators with
**no external dependencies at all**. Everything that is numerically heavy —
linear algebra, gate matrices and simulation — is written in V and compiled to a
native extension; the public API is pure Python on top of it.

The API follows [Cirq](https://quantumai.google/cirq) closely, so if you know
Cirq you already know vqubit.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1), vq.measure(q0, q1))
print(circuit)

result = vq.Simulator().run(circuit, repetitions=100)
print(result.histogram(key="result"))
```

## Why vqubit

- **No dependencies.** `pip install vqubit` installs nothing else. Everything
  runs locally, on your machine.
- **Native speed.** Linear algebra, gate matrices and simulation are compiled
  from V, not interpreted.
- **One wheel for every Python.** A single *abi3* wheel covers CPython
  3.11–3.15 per platform.
- **Cross-platform.** Linux and macOS (x86-64 and ARM).
- **Cirq-compatible API.** Build, transform and simulate circuits the way Cirq
  documents them.

## New to quantum computing?

You do not need to know any physics to start. If you have never heard of a
qubit, begin with [What is a qubit?](guide/what-is-a-qubit.md) and follow the
guide in order. Every concept is explained from scratch, with runnable examples.

## Next steps

1. [Install vqubit](installation.md).
2. Follow the [quickstart](quickstart.md) to run your first circuit.
3. Learn the [basics](guide/what-is-a-qubit.md).
4. Browse the [API reference](api/index.md).
