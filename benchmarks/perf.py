"""Regression gate for the native kernels.

``verify.py`` checks that the numbers are right; this checks that they are still
*fast*. Each case is timed three times and the best run is compared with
``perf_budget.json``, which holds an absolute ceiling per case.

The ceilings are deliberately generous (about twenty times the time measured on
a development machine, because a shared CI runner is several times slower and
noisy): the gate is there to catch a kernel that broke or became an order of
magnitude slower, not a few percent. Re-baseline with
``python benchmarks/perf.py --write`` on a quiet machine, with a release build
(``vcraft build --release``): a debug extension is several times slower.

Run it with ``python benchmarks/perf.py``; it exits non-zero when a case exceeds
its ceiling.
"""

from __future__ import annotations

import itertools
import json
import os
import sys
import time
from collections.abc import Callable
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import vqubit as vq

BUDGET = os.path.join(os.path.dirname(__file__), "perf_budget.json")
REPEATS = 3


def state_vector_gates() -> None:
	"""Sixteen qubits, and a couple of hundred single- and two-qubit gates."""
	qubits = vq.LineQubit.range(16)
	operations: list[Any] = []
	for _ in range(20):
		for qubit in qubits:
			operations.append(vq.H(qubit))
		for first, second in itertools.pairwise(qubits):
			operations.append(vq.CNOT(first, second))
	result = vq.Simulator().simulate(vq.Circuit(operations))
	result.state_vector()


def density_matrix_gates() -> None:
	"""Eight qubits as a density matrix."""
	qubits = vq.LineQubit.range(8)
	operations: list[Any] = []
	for _ in range(20):
		for qubit in qubits:
			operations.append(vq.H(qubit))
		for first, second in itertools.pairwise(qubits):
			operations.append(vq.CNOT(first, second))
	result = vq.DensityMatrixSimulator().simulate(vq.Circuit(operations))
	np_array = result.final_density_matrix
	del np_array


def cliford_tableau() -> None:
	"""A hundred qubits in the stabilizer formalism, five hundred gates."""
	qubits = vq.LineQubit.range(100)
	operations: list[Any] = []
	for index in range(5):
		for qubit in qubits:
			operations.append(vq.H(qubit) if index % 2 else vq.S(qubit))
		for first, second in itertools.pairwise(qubits):
			operations.append(vq.CNOT(first, second))
	vq.CliffordSimulator().simulate(vq.Circuit(operations))


def sampling() -> None:
	"""Sixteen qubits measured two thousand times."""
	qubits = vq.LineQubit.range(16)
	circuit = vq.Circuit(vq.H(qubits[0]), vq.CNOT(qubits[0], qubits[1]), vq.measure(*qubits, key="m"))
	vq.Simulator().run(circuit, repetitions=2000)


CASES: list[tuple[str, Callable[[], None]]] = [
	("state-vector 16q, 220 gates", state_vector_gates),
	("density-matrix 8q, 220 gates", density_matrix_gates),
	("stabilizer 100q, 500 gates", cliford_tableau),
	("samples 16q, 2000 runs", sampling),
]


def measure(case: Callable[[], None]) -> float:
	"""The best of ``REPEATS`` runs, in seconds."""
	best = float("inf")
	for _ in range(REPEATS):
		start = time.perf_counter()
		case()
		best = min(best, time.perf_counter() - start)
	return best


def main(argv: list[str]) -> int:
	if "--write" in argv:
		budget = {name: round(measure(case) * 20, 3) for name, case in CASES}
		with open(BUDGET, "w", encoding="utf-8") as handle:
			json.dump(budget, handle, indent="\t", sort_keys=True)
			handle.write("\n")
		print(f"wrote {len(budget)} ceilings to {BUDGET}")
		return 0

	if not os.path.exists(BUDGET):
		print(f"{BUDGET} is missing; run with --write on a quiet machine")
		return 1
	with open(BUDGET, encoding="utf-8") as handle:
		budget = json.load(handle)

	failures = 0
	print(f"{'case':32s} {'best':>8s} {'ceiling':>8s}  status")
	for name, case in CASES:
		best = measure(case)
		ceiling = float(budget.get(name, 0))
		ok = ceiling and best <= ceiling
		failures += 0 if ok else 1
		print(f"{name:32s} {best:8.3f} {ceiling:8.3f}  {'ok' if ok else 'SLOW'}")
	return 1 if failures else 0


if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))
