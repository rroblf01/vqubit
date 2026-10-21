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

## No `**` on operations

In Cirq anything can be raised to a power, operations included: `cirq.X(q) ** 0.5`
halves the gate. In vqubit `**` lives on the *gate*, so `vq.X ** 0.5` works but
`vq.X(q0) ** 0.5` raises `TypeError`. Exponentiate the gate, then call it.

```python
import vqubit as vq

q0 = vq.LineQubit(0)

assert vq.X ** 0.5 == vq.XPowGate(exponent=0.5)
assert (vq.X ** 0.5)(q0) == vq.XPowGate(exponent=0.5)(q0)
```

The plain gates **are** the power gates, as in Cirq: `vq.CZ` is a
`CZPowGate(exponent=1)`, so `isinstance(vq.CZ, vq.CZPowGate)` and
`vq.CZ == vq.CZPowGate()` are `True`. Likewise for `CNOT`/`CX` (a `CXPowGate`),
`SWAP`, `ISWAP`, `SQRT_ISWAP`, `CCX`/`TOFFOLI`/`CCNOT`, `CY`, `CCZ` and `CCY`.

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)

assert isinstance(vq.CZ, vq.CZPowGate)
assert vq.CZ == vq.CZPowGate(exponent=1)
assert vq.CZ ** 0.5 == vq.CZPowGate(exponent=0.5)
assert vq.SQRT_ISWAP ** 0.5 == vq.ISwapPowGate(exponent=0.25)
assert vq.rz(0.3) ** 0.5 == vq.RzGate(0.15)
assert vq.I ** 0.5 == vq.I
```

The historical names `vq.CZGate`, `vq.CNOTGate`, `vq.SWAPGate`, `vq.CCXGate` and
`vq.ISWAPGate` live on as **aliases** of the power classes, so code that mentions
them (in an `isinstance`, say) keeps working.

**Equality** is Cirq's too: it compares the eigenphases modulo the gate's period,
so `vq.X ** 3 == vq.X` and `vq.CNOT ** -1 == vq.CNOT`, as in Cirq. `vq.ISWAP` is
the case that shows it is not a blunt "modulo 2" on the exponent: one of its
eigenspaces sits at half a turn, so `vq.ISWAP ** -1 != vq.ISWAP`.

```python
import vqubit as vq

assert vq.X ** 3 == vq.X
assert vq.CNOT ** -1 == vq.CNOT
assert vq.ISWAP ** -1 != vq.ISWAP
```

Equality also requires the **same concrete class**: `vq.SqrtISwapGate()` and
`vq.SQRT_ISWAP` (an `ISwapPowGate`) have the same matrix but are not equal,
because `SqrtISwapGate` is a subclass vqubit adds for its `repr`. Cirq compares
against the decorated class and would let a subclass equal its base; vqubit is
stricter, and this is the only observable difference of that kind
(`SqrtISwapGate` does not exist in Cirq).

`vq.FREDKIN` and `vq.CSWAP` do not take fractional powers, exactly like
`cirq.CSWAP`; they do take the two powers every gate has, `** 1` (the gate
itself) and `** -1` (its inverse, which for CSWAP is CSWAP again).

One exception remains: the **Clifford** gates (`CZSWAP`, `CXSWAP`) compare their
matrix rather than their tableau, so `vq.CZSWAP ** 3` is not equal to
`vq.CZSWAP` even though Cirq says it is (the matrix's floating-point noise keeps
them apart). It is recorded, with its closing criterion, in `ROADMAP.md`.

## Parameters are `Parameter`, not `sympy.Symbol`

Cirq accepts `sympy` in its gates (`cirq.rz(sympy.Symbol("t"))`); vqubit has its
own type, `vq.Parameter` (or `vq.symbol("t")`), so that it does not depend on
`sympy`. A `sympy` symbol in a gate raises `TypeError` rather than converting
silently; it *is* usable as a **key** of a `vq.ParamResolver`, because it reduces
to its name.

```python
import sympy
import vqubit as vq

q0 = vq.LineQubit(0)
t = vq.symbol("t")

circuit = vq.Circuit(vq.rz(t)(q0))
resolved = vq.resolve_parameters(circuit, vq.ParamResolver({sympy.Symbol("t"): 0.5}))
assert vq.allclose_up_to_global_phase(vq.unitary(resolved), vq.unitary(vq.Circuit(vq.rz(0.5)(q0))))
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

## Parameter sweeps and `sample`

`vq.Simulator().simulate(circuit, resolver)` and
`vq.Simulator().run(circuit, resolver, repetitions=100)` take the
`ParamResolver` where Cirq takes it, and `run_sweep`, `simulate_sweep`,
`simulate_moment_steps` and `run_batch` have Cirq's signatures. Two things
change:

- **The sweep classes have a different signature.** vqubit writes
  `vq.Linspace({"t": 0.0}, {"t": 1.0}, 4)` (two assignments to interpolate
  between, and it handles several parameters at once) and
  `vq.Points([{"t": 0.0}, {"t": 1.0}])`, where Cirq writes
  `cirq.Linspace("t", 0.0, 1.0, 4)` and `cirq.Points("t", [0.0, 1.0])`.
- **`sample` returns a `Result`, not a `DataFrame`.** Cirq returns a
  `pandas.DataFrame` and takes `params=`; vqubit does not depend on pandas, so it
  returns its own `Result` (just like `run`) and sweeps are written with
  `run_sweep`.

```python
import vqubit as vq

q0 = vq.LineQubit(0)
t = vq.Parameter("t")
circuit = vq.Circuit(vq.rx(t)(q0), vq.measure(q0, key="m"))

results = vq.Simulator(seed=0).run_sweep(circuit, vq.Points([{"t": 0.0}, {"t": 0.5}]), repetitions=10)
assert [dict(result.params)["t"] for result in results] == [0.0, 0.5]
assert all(len(result) == 10 for result in results)
```

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

## The Clifford simulator and the gates it accepts

`vq.CliffordSimulator` accepts the same gates Cirq treats as stabilizer
operations: the elementary ones (`H`, `S`, `X`, `Y`, `Z`, `CNOT`, `CZ`, `SWAP`)
with an integer exponent and, since 1.0, also `X ** 0.5`, `Y ** 0.5`, `XX`,
`YY`, `ZZ` (integer or half-integer exponent) and `iSWAP` (integer exponent).
The ones the native tableau cannot spell are rewritten as short words over the
nine instructions it does have; the state matches both Cirq's stabilizer
simulator and vqubit's own `Simulator`. `vq.has_stabilizer_effect` answers the
same as `cirq.has_stabilizer_effect` for all of them.

One gate stays out: a **bare `CliffordGate`** (say `vq.CXSWAP`, or one built by
hand with `CliffordGate.from_unitary`). `has_stabilizer_effect` says `True` — it
is one — but the stabilizer simulator refuses it, because the native tableau
only takes instructions and a `CliffordGate` is an arbitrary matrix; simulate it
with `vq.Simulator`. The tableau is still yours: `vq.CliffordTableau` and
`vq.CliffordGate` work with it. It is recorded as a vcraft improvement in
`UPSTREAM.md` §7d.

## The text diagram

`vq.Circuit` draws one line per qubit and one column per moment, and it matches
Cirq's for the shared gates, exponents included: `vq.CZ ** 0.5` draws `@` over
`@^0.5`, just like `cirq.CZ ** 0.5`. Two details still differ, by design of the
renderer:

- **No tableau blocks.** Cirq draws a `CliffordGate` as a multi-line
  "Clifford Gate with Tableau" block; vqubit draws the one-line symbol
  `CliffordGate(2q)`.
- **Decimal exponents.** vqubit writes every exponent as a decimal (`Z^0.125`)
  where Cirq writes the fractions it recognizes as such (`Z^(1/8)`).

Both sit in the harness's list of tolerated divergences, with
`plan = "permanent"`.

## What is unchanged

Everything else: qubits (`LineQubit`, `GridQubit`, `NamedQubit`), the gate
classes and their arithmetic, `Circuit`/`Moment` behaviour and diagrams, the
simulation protocols, `PauliString` and its algebra, the transformers and
optimizers, `to_json`/`read_json` and its `cirq_type` format, the OpenQASM 2.0
writer, the device classes and the `Engine`.

## Next steps

- Meet the [gates](quantum-gates.md) and build [circuits](circuits.md).
- Read the [API reference](../api/index.md).
