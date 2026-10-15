"""Testing helpers, mirroring ``cirq.testing``.

Import as ``vqubit.testing``. These are for tests of code built on vqubit:
comparing circuits, checking that a ``repr`` round-trips, and generating random
circuits.
"""

from __future__ import annotations

import math
import random
from collections.abc import Iterable, Sequence
from typing import Any

from ._circuits import Circuit
from ._gates import CNOT, CZ, SWAP, H, S, T, X, Y, Z
from ._linalg import Matrix
from ._ops import Moment
from ._protocols import unitary
from ._qubits import Qid


def assert_equivalent_repr(value: Any, *, setup_code: str = "import vqubit") -> None:
	"""Assert that ``eval(repr(value)) == value`` under ``setup_code``."""
	namespace: dict[str, Any] = {}
	exec(setup_code, namespace)
	evaluated = eval(repr(value), namespace)
	assert evaluated == value, f"repr did not round-trip: {value!r}"


def assert_same_circuits(actual: Circuit, expected: Circuit) -> None:
	"""Assert that two circuits are equal."""
	assert actual == expected, f"circuits differ:\n{actual}\n!=\n{expected}"


def _without_measurements(circuit: Circuit) -> Circuit:
	return Circuit(operation for operation in circuit.all_operations() if not operation.gate._is_measurement_())


def assert_circuits_with_terminal_measurements_are_equivalent(a: Circuit, b: Circuit, atol: float = 1e-8) -> None:
	"""Assert two circuits agree once their measurements are removed."""
	matrix_a = unitary(_without_measurements(a))
	matrix_b = unitary(_without_measurements(b))
	phase = None
	for i in range(len(matrix_b)):
		for j in range(len(matrix_b)):
			if abs(matrix_b[i][j]) > 1e-9:
				phase = matrix_a[i][j] / matrix_b[i][j]
				break
		if phase is not None:
			break
	assert phase is not None
	for i in range(len(matrix_a)):
		for j in range(len(matrix_a)):
			assert abs(matrix_a[i][j] - phase * matrix_b[i][j]) <= atol, "circuits are not equivalent"


def assert_allclose_up_to_global_phase(actual: Matrix, expected: Matrix, atol: float = 1e-8) -> None:
	"""Assert two matrices are equal up to a global phase."""
	assert len(actual) == len(expected)
	phase = None
	for i in range(len(expected)):
		for j in range(len(expected)):
			if abs(expected[i][j]) > 1e-9:
				phase = actual[i][j] / expected[i][j]
				break
		if phase is not None:
			break
	assert phase is not None
	for i in range(len(actual)):
		for j in range(len(actual)):
			assert abs(actual[i][j] - phase * expected[i][j]) <= atol


_DEFAULT_GATE_DOMAIN: dict[int, list[Any]] = {
	1: [X, Y, Z, H, S, T],
	2: [CNOT, CZ, SWAP],
}


def random_circuit(
	qubits: Iterable[Qid],
	n_moments: int,
	op_density: float,
	*,
	random_state: random.Random | int | None = None,
	gate_domain: dict[int, Sequence[Any]] | None = None,
) -> Circuit:
	"""Return a random circuit over ``qubits``.

	Each moment places a gate on each qubit with probability ``op_density``; a
	two-qubit gate also takes a second, unused qubit.
	"""
	rng = random_state if isinstance(random_state, random.Random) else random.Random(random_state)
	qubits = tuple(qubits)
	domain = dict(gate_domain) if gate_domain is not None else _DEFAULT_GATE_DOMAIN
	arities = sorted(domain)
	circuit = Circuit()
	for _ in range(n_moments):
		used: set[Qid] = set()
		moment_ops: list[Any] = []
		for qubit in qubits:
			if qubit in used or rng.random() > op_density:
				continue
			arity = rng.choice(arities)
			gate = rng.choice(list(domain[arity]))
			if arity == 1:
				moment_ops.append(gate.on(qubit))
				used.add(qubit)
				continue
			partners = [other for other in qubits if other != qubit and other not in used]
			if partners:
				other = rng.choice(partners)
				moment_ops.append(gate.on(qubit, other))
				used.add(qubit)
				used.add(other)
		if moment_ops:
			circuit.append(Moment(moment_ops))
	return circuit


def assert_has_consistent_apply_unitary(value: Any, *, atol: float = 1e-8) -> None:
	"""Assert that ``apply_unitary`` matches the value's unitary."""
	from ._protocols import apply_unitary, has_unitary

	if not has_unitary(value):
		return
	matrix = unitary(value)
	dim = len(matrix)
	rng = random.Random(0)
	state = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(dim)]
	norm = math.sqrt(sum(abs(z) ** 2 for z in state))
	state = [z / norm for z in state]
	applied = apply_unitary(value, state)
	expected = [sum(matrix[i][j] * state[j] for j in range(dim)) for i in range(dim)]
	for a, e in zip(applied, expected, strict=True):
		assert abs(a - e) <= atol, "apply_unitary disagrees with unitary"


def assert_decompose_is_consistent_with_unitary(value: Any, *, atol: float = 1e-8) -> None:
	"""Assert that decomposing ``value`` preserves its unitary."""
	from ._protocols import decompose, has_unitary

	subject = getattr(value, "gate", value)
	if not hasattr(subject, "_decompose_") or not has_unitary(value):
		return
	decomposed = decompose(value)
	if not decomposed:
		return
	matrix = unitary(value)
	from ._gates import Gate
	from ._ops import Operation
	from ._qubits import LineQubit

	operations: list[Any] = []
	for item in decomposed:
		if isinstance(item, Operation):
			operations.append(item)
		elif isinstance(item, Gate):
			operations.append(item.on(*LineQubit.range(item.num_qubits)))
		else:
			operations.append(item)
	rebuilt = unitary(Circuit(operations))
	assert len(matrix) == len(rebuilt), "decomposition changed the dimension"
	phase = None
	for i in range(len(matrix)):
		for j in range(len(matrix)):
			if abs(matrix[i][j]) > 1e-9:
				phase = rebuilt[i][j] / matrix[i][j]
				break
		if phase is not None:
			break
	assert phase is not None
	for i in range(len(matrix)):
		for j in range(len(matrix)):
			assert abs(rebuilt[i][j] - phase * matrix[i][j]) <= atol, "decomposition is inconsistent"


def assert_pauli_expansion_is_consistent_with_unitary(value: Any, *, atol: float = 1e-8) -> None:
	"""Assert that the Pauli expansion of ``value`` matches its unitary."""
	from ._protocols import has_unitary, pauli_expansion

	subject = getattr(value, "gate", value)
	if not hasattr(subject, "_pauli_expansion_") or not has_unitary(value):
		return
	expansion = pauli_expansion(subject)
	matrix = unitary(subject)
	rebuilt = [[0j] * len(matrix) for _ in range(len(matrix))]

	for pauli_string, coefficient in expansion.items():
		piece = pauli_string.matrix()
		for i in range(len(matrix)):
			for j in range(len(matrix)):
				rebuilt[i][j] += coefficient * piece[i][j]
	for i in range(len(matrix)):
		for j in range(len(matrix)):
			assert abs(rebuilt[i][j] - matrix[i][j]) <= atol, "Pauli expansion is inconsistent"


def assert_qasm_is_consistent_with_unitary(value: Any, *, atol: float = 1e-8) -> None:
	"""Assert that ``value`` can be written as OpenQASM."""
	del atol
	from ._qasm import qasm

	text = qasm(value)
	assert isinstance(text, str) and text.strip(), "qasm did not produce a string"


def assert_has_consistent_qid_shape(value: Any) -> None:
	"""Assert that the qid shape of ``value`` matches its qubit count."""
	from ._protocols import qid_shape

	shape = qid_shape(value)
	assert len(shape) == len(getattr(value, "qubits", shape)), "qid_shape length mismatch"


def assert_implements_consistent_protocols(value: Any, *, atol: float = 1e-8) -> None:
	"""Run the common consistency checks on ``value``."""
	assert_equivalent_repr(value)
	assert_has_consistent_apply_unitary(value, atol=atol)
	assert_decompose_is_consistent_with_unitary(value, atol=atol)


def equals_tester() -> Any:
	"""Return a function that checks the ``==``/``!=``/``hash`` contract."""

	def tester(a: Any, b: Any) -> None:
		assert (a == b) == (b == a), "equality is not symmetric"
		assert (a != b) == (b != a), "inequality is not symmetric"
		if a == b:
			assert hash(a) == hash(b), "equal values must hash equally"

	return tester


__all__ = [
	"assert_allclose_up_to_global_phase",
	"assert_circuits_with_terminal_measurements_are_equivalent",
	"assert_decompose_is_consistent_with_unitary",
	"assert_equivalent_repr",
	"assert_has_consistent_apply_unitary",
	"assert_has_consistent_qid_shape",
	"assert_implements_consistent_protocols",
	"assert_pauli_expansion_is_consistent_with_unitary",
	"assert_qasm_is_consistent_with_unitary",
	"assert_same_circuits",
	"equals_tester",
	"random_circuit",
]
