# Differences from Cirq

vqubit follows Cirq closely. Measured against **cirq-core 1.7.0**, it exports
**606 of its 607 top-level names** — the missing one is `contrib`, which is a
package of third-party integrations — and its circuits, gates, simulators,
transformers and devices follow the same design, so code you know usually works
after `import vqubit as vq`.

It is **not** a drop-in replacement, though. This page lists the places where the
call convention differs, so that ported Cirq code raises an error instead of
quietly doing something else.

## Renamed arguments

The most common difference is spelling: vqubit writes parameter names out in
full. Calls that pass values **positionally** are unaffected; calls that use the
keyword must be updated.

| Cirq         | vqubit        | Example                                     |
| ------------ | ------------- | ------------------------------------------- |
| `p`          | `probability` | `vq.bit_flip(probability=0.1)`              |
| `rads`       | `radians`     | `vq.rx(radians=math.pi / 2)`                |
| `n_qubits`   | `num_qubits`  | `vq.depolarize(0.1, num_qubits=2)`          |
| `val`        | `value`       | `vq.unitary(value=vq.X)`                    |
| `mat`        | `matrix`      | `vq.single_qubit_matrix_to_gates(matrix=u)` |
| `contents`   | `operations`  | `vq.Moment(operations=[vq.H(q0)])`          |
| `tolerance`  | `atol`        | `vq.kak_decomposition(matrix, atol=1e-8)`   |

Gate sizes are attributes, not methods: `vq.H.num_qubits` is `1` where Cirq has
`cirq.H.num_qubits()`. Operations still have the method, so
`vq.H(q0).num_qubits() == 1`.

## Different argument order

A few calls take their arguments in another order, or take the matrix first.
These raise `TypeError` if you port the call verbatim.

| Purpose             | Cirq                                        | vqubit                                             |
| ------------------- | ------------------------------------------- | -------------------------------------------------- |
| Measurement gate    | `MeasurementGate(num_qubits, key)`          | `MeasurementGate(key, num_qubits)`                 |
| Moment              | `Moment(op1, op2)`                          | `Moment([op1, op2])`                               |
| Two-qubit synthesis | `two_qubit_matrix_to_operations(q0, q1, u)` | `two_qubit_matrix_to_operations(u, qubits=(q0, q1))` |
| Three-qubit synthesis | `three_qubit_matrix_to_operations(q0, q1, q2, u)` | `three_qubit_matrix_to_operations([q0, q1, q2], u)` |
| Classical control   | `ClassicallyControlledOperation(op, conditions)` | `ClassicallyControlledOperation(conditions, op)` |

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)

gate = vq.MeasurementGate("m", 2)   # the key comes first
assert gate.key == "m"

moment = vq.Moment([vq.H(q0), vq.X(q1)])   # one iterable, not *operations
assert len(moment) == 2
```

## No `**` on operations, and the plain gates are not extrapolatable

In Cirq anything can be raised to a power: `cirq.X(q) ** 0.5` halves the gate.
In vqubit `**` lives on the *exponentiated* gate classes, so `vq.X ** 0.5` works
but `vq.X(q0) ** 0.5` raises `TypeError`. Exponentiate the gate, then call it.

```python
import vqubit as vq

q0 = vq.LineQubit(0)

assert vq.X ** 0.5 == vq.XPowGate(exponent=0.5)
assert (vq.X ** 0.5)(q0) == vq.XPowGate(exponent=0.5)(q0)
```

The same holds for the gates whose Cirq counterparts are exponentiated:
`vq.CZ`, `vq.CNOT`, `vq.SWAP`, `vq.ISWAP` and `vq.rz(...)` are the plain gates
(`CZGate`, `CNOTGate`, `SWAPGate`, `ISWAPGate`, `RzGate`) and define no `**`.
Use the `Pow` classes to interpolate between the identity and the gate:

```python
import vqubit as vq

assert vq.CZPowGate(exponent=0.5) == vq.CZPowGate() ** 0.5
assert vq.SwapPowGate(exponent=0.5) == vq.SwapPowGate() ** 0.5
```

## The simulators take only `seed`

`vq.Simulator` and `vq.DensityMatrixSimulator` accept `seed` and nothing else:
there is no `noise=` and no `dtype=`. Noise travels *inside the circuit*, as a
channel, and needs the density-matrix simulator.

```python
import vqubit as vq

q0 = vq.LineQubit(0)
noisy = vq.Circuit(vq.H(q0), vq.depolarize(0.05)(q0), vq.measure(q0, key="m"))

result = vq.DensityMatrixSimulator(seed=0).run(noisy, repetitions=100)
assert sum(result.histogram(key="m").values()) == 100
```

State vectors and density matrices are always `complex128`. (Cirq's
density-matrix simulator defaults to `complex64`, so vqubit is more precise and
uses twice the memory for the same circuit.)

## Same name, different role

A handful of names exist in both libraries but mean something else. They fail
loudly rather than silently, but you do have to adapt the call.

| Name | Cirq | vqubit |
| ---- | ---- | ------ |
| `final_state_vector` | takes a **circuit** | takes a **simulation result** |
| `final_density_matrix` | takes a **circuit** | takes a **simulation result** |
| `read_json` | `read_json(file_or_fn, json_text)` | `read_json(text)` |
| `qasm` | `qasm(val, args, qubits)` | `qasm(value)` |
| `measure` | `measure(q, invert_mask=…)` | `measure(*qubits, key=None)`, the mask lives on `MeasurementGate` |
| `entanglement_fidelity` | fidelity of a **channel** | fidelity between two **states** |
| `xeb_fidelity` | `xeb_fidelity(circuit, bitstrings, …)` | `xeb_fidelity(bitstrings, hilbert_space_dimension)` |
| `estimate_single_qubit_readout_errors` | takes a **sampler** | takes **measurement results** |
| `prepare_two_qubit_state_using_cz` | `(q0, q1, state)` | `(gate, qubits)` |

For a final state the idiomatic vqubit call goes through the simulator:

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
circuit = vq.Circuit(vq.H(q0), vq.CNOT(q0, q1))

state = vq.Simulator().simulate(circuit).state_vector()
assert abs(abs(state[0]) - 2**-0.5) < 1e-9
assert abs(abs(state[3]) - 2**-0.5) < 1e-9
```

## What is unchanged

Everything else: qubits (`LineQubit`, `GridQubit`, `NamedQubit`), the gate
classes and their arithmetic, `Circuit`/`Moment` behaviour and diagrams, the
simulation protocols, `PauliString` and its algebra, the transformers and
optimizers, `to_json`/`read_json` and its `cirq_type` format, the OpenQASM 2.0
writer, the device classes and the `Engine`.

## Next steps

- Meet the [gates](quantum-gates.md) and build [circuits](circuits.md).
- Read the [API reference](../api/index.md).
