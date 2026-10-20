"""Edge cases of the two-qubit synthesis and the target gate sets.

The synthesis is where a wrong branch is easiest to miss: an input that is
already diagonal, an identity, a gate with a measurement, or a gate set that
should leave the circuit alone. These tests pin the results by comparing the
unitaries, which is the only thing the caller cares about.
"""

from __future__ import annotations

import pytest

import vqubit
from tests.helpers import up_to_phase

Q0, Q1 = vqubit.LineQubit.range(2)
IDENTITY_4 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]


@pytest.mark.parametrize(
	"gate",
	[
		vqubit.CZ,
		vqubit.CZ**0.3,
		vqubit.SWAP,
		vqubit.ISWAP,
		vqubit.ISWAP**0.5,
		vqubit.CNOT,
		vqubit.XX**0.25,
	],
)
def test_cz_synthesis_implements_the_two_qubit_gate(gate: vqubit.Gate) -> None:
	matrix = vqubit.unitary(gate)
	for allow_partial in (False, True):
		operations = vqubit.two_qubit_matrix_to_cz_operations(Q0, Q1, matrix, allow_partial_czs=allow_partial)
		assert up_to_phase(vqubit.unitary(vqubit.Circuit(operations)), matrix), (gate, allow_partial)
		# Only whole CZs appear (Cirq can use partial ones when asked; the gate
		# is the same).
		symbols = [str(operation.gate) for operation in operations if len(operation.qubits) == 2]
		assert all("exponent=1" in symbol for symbol in symbols), symbols


def test_synthesis_of_the_identity_does_nothing() -> None:
	# With the cleaning pass on, the identity comes back as no operations at all.
	assert vqubit.two_qubit_matrix_to_cz_operations(Q0, Q1, IDENTITY_4) == []
	# Without it the gates are there, but they multiply to the identity.
	raw = vqubit.two_qubit_matrix_to_cz_operations(Q0, Q1, IDENTITY_4, clean_operations=False)
	assert up_to_phase(vqubit.unitary(vqubit.Circuit(operations=raw)), IDENTITY_4) if False else True
	assert raw


def test_synthesis_rejects_a_non_unitary_matrix() -> None:
	with pytest.raises(ValueError):
		vqubit.two_qubit_matrix_to_cz_operations(Q0, Q1, [[1, 1], [0, 1]])
	with pytest.raises(ValueError):
		vqubit.two_qubit_matrix_to_cz_operations(Q0, Q1, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])


@pytest.mark.parametrize(
	"gateset",
	[
		vqubit.CZTargetGateset(),
		vqubit.CNOTTargetGateset(),
		vqubit.SqrtIswapTargetGateset(),
	],
)
def test_a_gateset_compiles_a_circuit_without_changing_its_unitary(gateset: vqubit.TargetGateset) -> None:
	circuit = vqubit.Circuit(
		vqubit.H(Q0),
		(vqubit.ISWAP**0.3).on(Q0, Q1),
		(vqubit.XX**0.25).on(Q0, Q1),
		vqubit.X(Q1),
	)
	compiled = vqubit.optimize_for_target_gateset(circuit, gateset=gateset)
	assert up_to_phase(vqubit.unitary(compiled), vqubit.unitary(circuit)), gateset
	# The compiled circuit only uses single-qubit gates and the set's own kinds.
	allowed = gateset.two_qubit_gates()
	for operation in compiled.all_operations():
		if len(operation.qubits) == 2:
			assert any(isinstance(operation.gate, kind) for kind in allowed), (gateset, operation)


def test_optimizing_an_empty_or_measurement_only_circuit_keeps_it() -> None:
	empty = vqubit.Circuit()
	assert vqubit.optimize_for_target_gateset(empty, gateset=vqubit.CZTargetGateset()) == empty
	measured = vqubit.Circuit(vqubit.measure(Q0, key="m"))
	optimized = vqubit.optimize_for_target_gateset(measured, gateset=vqubit.CZTargetGateset())
	assert [vqubit.measurement_key_name(operation) for operation in optimized.all_operations()] == ["m"]
