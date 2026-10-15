# Optimization

A **transformer** takes a circuit and returns an equivalent, cleaner or cheaper
one. The optimizers live in the `vqubit` namespace.

## Cleaning up circuits

```python
import vqubit as vq

q0, q1 = vq.LineQubit.range(2)
c = vq.Circuit(vq.H(q0), vq.I(q0), vq.CNOT(q0, q1))

vq.drop_empty_moments(c)              # remove empty moments
vq.drop_negligible_operations(c)      # remove gates that do nothing
vq.expand_composite(c)                # expand composite gates (CCX, ...)
```

## Merging operations

`merge_single_qubit_moments_to_phxz` combines runs of single-qubit gates into a
single `PhasedXZGate`, the native gate of many processors:

```python
c = vq.Circuit(vq.X(q0), vq.H(q0), vq.Rz(0.3)(q0))
c = vq.merge_single_qubit_moments_to_phxz(c)
print(c)   # a single PhXZ gate
```

## Decomposing a single-qubit unitary

Any one-qubit gate can be written as a sequence of rotations:

```python
matrix = vq.unitary(vq.H)
print(vq.single_qubit_matrix_to_gates(matrix))          # [Rz, Ry, Rz]
print(vq.single_qubit_matrix_to_phased_x_z(matrix))     # [PhasedX, Z]
print(vq.single_qubit_matrix_to_phxz(matrix))           # a PhasedXZGate
```

## Compiling to a gate set

A **target gate set** says which gates a device can run.
`optimize_for_target_gateset` rewrites the circuit to use them:

```python
q0, q1 = vq.LineQubit.range(2)
c = vq.Circuit(vq.CNOT(q0, q1), vq.H(q0))

compiled = vq.optimize_for_target_gateset(c, vq.CZTargetGateset())
print(compiled)   # only CZ and PhasedXZ gates
```

An arbitrary two-qubit gate is synthesized with the **KAK** decomposition:

```python
import vqubit as vq
from vqubit import two_qubit_matrix_to_cz_operations

q0, q1 = vq.LineQubit.range(2)
ops = two_qubit_matrix_to_cz_operations(q0, q1, vq.unitary(vq.ISWAP))
```

## Next steps

- Save circuits with [serialization](serialization.md).
- Run them on local [devices](devices.md).
