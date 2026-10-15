"""Testing helpers, mirroring ``cirq.testing``.

Import as ``vqubit.testing``. These are for tests of code built on vqubit:
comparing circuits, checking that a ``repr`` round-trips, and generating random
circuits.
"""

from __future__ import annotations

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


__all__ = [
	"assert_allclose_up_to_global_phase",
	"assert_circuits_with_terminal_measurements_are_equivalent",
	"assert_equivalent_repr",
	"assert_same_circuits",
	"random_circuit",
]
