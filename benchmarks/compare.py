"""Benchmark vqubit against Cirq on the same tasks.

Run ``python benchmarks/compare.py`` (with Cirq installed) to print a table of
wall time and peak memory for each task, for both packages. The script measures
each task in a fresh subprocess so that peak RSS (``VmHWM``) is not polluted by
earlier tasks; time is the best of several runs inside that subprocess.

Only the standard library is used. If Cirq is not installed the vqubit column
is still produced and the Cirq column is left blank.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import random
import subprocess
import sys
import time
from collections.abc import Callable

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import vqubit as vq

try:
	import cirq

	HAS_CIRQ = True
except ImportError:  # pragma: no cover - depends on the environment
	cirq = None
	HAS_CIRQ = False

# --- task parameters ---------------------------------------------------------

N_MOMENTS = 300
SIM_QUBITS = 12
SIM_QUBITS_LARGE = 20
DENSITY_QUBITS = 7
DENSITY_QUBITS_LARGE = 10
REPETITIONS = 1000
UNITARY_QUBITS = 6
QASM_GATES = 50
PAULI_MULTIPLY_COUNT = 1000
REPEATS = 2


# --- shared data -------------------------------------------------------------


def _random_unitary(dim: int, rng: random.Random) -> list[list[complex]]:
	columns: list[list[complex]] = []
	for _ in range(dim):
		vector = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(dim)]
		for column in columns:
			overlap = sum(column[k].conjugate() * vector[k] for k in range(dim))
			vector = [vector[k] - overlap * column[k] for k in range(dim)]
		norm = sum(abs(z) ** 2 for z in vector) ** 0.5
		columns.append([z / norm for z in vector])
	return [[columns[j][i] for j in range(dim)] for i in range(dim)]


_RNG = random.Random(1234)
_MATRIX_4 = _random_unitary(4, _RNG)
_MATRIX_2 = _random_unitary(2, _RNG)


# --- vqubit tasks ------------------------------------------------------------


def vq_build_circuit() -> None:
	qubits = vq.LineQubit.range(5)
	circuit = vq.Circuit()
	for _ in range(N_MOMENTS):
		circuit.append([vq.H(q) for q in qubits])
		circuit.append([vq.CNOT(qubits[i], qubits[i + 1]) for i in range(4)])


def vq_unitary() -> None:
	qubits = vq.LineQubit.range(UNITARY_QUBITS)
	circuit = vq.Circuit(vq.H(q) for q in qubits)
	for i in range(UNITARY_QUBITS - 1):
		circuit.append(vq.CNOT(qubits[i], qubits[i + 1]))
	vq.unitary(circuit)


def vq_simulate_state_vector() -> None:
	qubits = vq.LineQubit.range(SIM_QUBITS)
	circuit = vq.Circuit(vq.H(qubits[0]))
	for i in range(SIM_QUBITS - 1):
		circuit.append(vq.CNOT(qubits[i], qubits[i + 1]))
	vq.Simulator().simulate(circuit).state_vector()


def vq_sample() -> None:
	qubits = vq.LineQubit.range(SIM_QUBITS)
	circuit = vq.Circuit(vq.H(qubits[0]))
	for i in range(SIM_QUBITS - 1):
		circuit.append(vq.CNOT(qubits[i], qubits[i + 1]))
	circuit.append(vq.measure(*qubits, key="m"))
	vq.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)


def vq_density_matrix_simulate() -> None:
	qubits = vq.LineQubit.range(DENSITY_QUBITS)
	circuit = vq.Circuit(vq.H(qubits[0]))
	for i in range(DENSITY_QUBITS - 1):
		circuit.append(vq.CNOT(qubits[i], qubits[i + 1]))
	circuit.append(vq.depolarize(0.01).on_each(*qubits))
	vq.DensityMatrixSimulator().simulate(circuit).density_matrix()


def vq_simulate_large() -> None:
	qubits = vq.LineQubit.range(SIM_QUBITS_LARGE)
	circuit = vq.Circuit(vq.H(qubits[0]))
	for i in range(SIM_QUBITS_LARGE - 1):
		circuit.append(vq.CNOT(qubits[i], qubits[i + 1]))
	vq.Simulator().simulate(circuit).state_vector()


def vq_density_large() -> None:
	qubits = vq.LineQubit.range(DENSITY_QUBITS_LARGE)
	circuit = vq.Circuit(vq.H(qubits[0]))
	for i in range(DENSITY_QUBITS_LARGE - 1):
		circuit.append(vq.CNOT(qubits[i], qubits[i + 1]))
	circuit.append(vq.depolarize(0.01).on_each(*qubits))
	vq.DensityMatrixSimulator().simulate(circuit).density_matrix()


def vq_kak() -> None:
	vq.kak_decomposition(_MATRIX_4)


def vq_single_qubit_decompose() -> None:
	vq.single_qubit_matrix_to_phased_x_z(_MATRIX_2)


def vq_pauli_multiply() -> None:
	q0, q1 = vq.LineQubit.range(2)
	left = vq.PauliString({q0: vq.Pauli.X, q1: vq.Pauli.Z})
	right = vq.PauliString({q0: vq.Pauli.Y, q1: vq.Pauli.X})
	for _ in range(PAULI_MULTIPLY_COUNT):
		left * right


def vq_qasm() -> None:
	qubits = vq.LineQubit.range(5)
	circuit = vq.Circuit()
	for i in range(QASM_GATES):
		circuit.append(vq.H(qubits[i % 5]))
		if i % 3 == 0:
			circuit.append(vq.CNOT(qubits[i % 5], qubits[(i + 1) % 5]))
	vq.qasm(circuit)


def vq_bloch_vector() -> None:
	qubits = vq.LineQubit.range(SIM_QUBITS)
	circuit = vq.Circuit(vq.H(qubits[0]))
	for i in range(SIM_QUBITS - 1):
		circuit.append(vq.CNOT(qubits[i], qubits[i + 1]))
	vector = vq.Simulator().simulate(circuit).state_vector()
	vq.bloch_vector_from_state_vector(vector, 0)


# --- cirq tasks --------------------------------------------------------------


def cq_build_circuit() -> None:
	qubits = cirq.LineQubit.range(5)
	circuit = cirq.Circuit()
	for _ in range(N_MOMENTS):
		circuit.append([cirq.H(q) for q in qubits])
		circuit.append([cirq.CNOT(qubits[i], qubits[i + 1]) for i in range(4)])


def cq_unitary() -> None:
	qubits = cirq.LineQubit.range(UNITARY_QUBITS)
	circuit = cirq.Circuit(cirq.H(q) for q in qubits)
	for i in range(UNITARY_QUBITS - 1):
		circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
	cirq.unitary(circuit)


def cq_simulate_state_vector() -> None:
	qubits = cirq.LineQubit.range(SIM_QUBITS)
	circuit = cirq.Circuit(cirq.H(qubits[0]))
	for i in range(SIM_QUBITS - 1):
		circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
	cirq.Simulator().simulate(circuit).state_vector()


def cq_sample() -> None:
	qubits = cirq.LineQubit.range(SIM_QUBITS)
	circuit = cirq.Circuit(cirq.H(qubits[0]))
	for i in range(SIM_QUBITS - 1):
		circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
	circuit.append(cirq.measure(*qubits, key="m"))
	cirq.Simulator(seed=0).run(circuit, repetitions=REPETITIONS)


def cq_density_matrix_simulate() -> None:
	qubits = cirq.LineQubit.range(DENSITY_QUBITS)
	circuit = cirq.Circuit(cirq.H(qubits[0]))
	for i in range(DENSITY_QUBITS - 1):
		circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
	circuit.append(cirq.depolarize(0.01).on_each(*qubits))
	_ = cirq.DensityMatrixSimulator().simulate(circuit).final_density_matrix


def cq_simulate_large() -> None:
	qubits = cirq.LineQubit.range(SIM_QUBITS_LARGE)
	circuit = cirq.Circuit(cirq.H(qubits[0]))
	for i in range(SIM_QUBITS_LARGE - 1):
		circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
	cirq.Simulator().simulate(circuit).state_vector()


def cq_density_large() -> None:
	qubits = cirq.LineQubit.range(DENSITY_QUBITS_LARGE)
	circuit = cirq.Circuit(cirq.H(qubits[0]))
	for i in range(DENSITY_QUBITS_LARGE - 1):
		circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
	circuit.append(cirq.depolarize(0.01).on_each(*qubits))
	_ = cirq.DensityMatrixSimulator().simulate(circuit).final_density_matrix


def cq_kak() -> None:
	import numpy as np

	cirq.kak_decomposition(np.array(_MATRIX_4))


def cq_single_qubit_decompose() -> None:
	import numpy as np

	cirq.single_qubit_matrix_to_phased_x_z(np.array(_MATRIX_2))


def cq_pauli_multiply() -> None:
	q0, q1 = cirq.LineQubit.range(2)
	left = cirq.PauliString({q0: cirq.X, q1: cirq.Z})
	right = cirq.PauliString({q0: cirq.Y, q1: cirq.X})
	for _ in range(PAULI_MULTIPLY_COUNT):
		left * right


def cq_qasm() -> None:
	qubits = cirq.LineQubit.range(5)
	circuit = cirq.Circuit()
	for i in range(QASM_GATES):
		circuit.append(cirq.H(qubits[i % 5]))
		if i % 3 == 0:
			circuit.append(cirq.CNOT(qubits[i % 5], qubits[(i + 1) % 5]))
	cirq.qasm(circuit)


def cq_bloch_vector() -> None:
	qubits = cirq.LineQubit.range(SIM_QUBITS)
	circuit = cirq.Circuit(cirq.H(qubits[0]))
	for i in range(SIM_QUBITS - 1):
		circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
	vector = cirq.Simulator().simulate(circuit).state_vector()
	cirq.bloch_vector_from_state_vector(vector, 0)


TASKS: dict[str, tuple[str, Callable[[], None], Callable[[], None] | None]] = {
	"build_circuit": ("Build a 5-qubit, 300-moment circuit", vq_build_circuit, cq_build_circuit if HAS_CIRQ else None),
	"unitary": (f"Unitary of a {UNITARY_QUBITS}-qubit circuit", vq_unitary, cq_unitary if HAS_CIRQ else None),
	"simulate": (
		f"State-vector simulate {SIM_QUBITS} qubits",
		vq_simulate_state_vector,
		cq_simulate_state_vector if HAS_CIRQ else None,
	),
	"sample": (f"Sample {REPETITIONS} shots of {SIM_QUBITS} qubits", vq_sample, cq_sample if HAS_CIRQ else None),
	"density_matrix": (
		f"Density-matrix simulate {DENSITY_QUBITS} qubits",
		vq_density_matrix_simulate,
		cq_density_matrix_simulate if HAS_CIRQ else None,
	),
	"simulate_large": (
		f"State-vector simulate {SIM_QUBITS_LARGE} qubits",
		vq_simulate_large,
		cq_simulate_large if HAS_CIRQ else None,
	),
	"density_large": (
		f"Density-matrix simulate {DENSITY_QUBITS_LARGE} qubits",
		vq_density_large,
		cq_density_large if HAS_CIRQ else None,
	),
	"kak": ("KAK decomposition of a random 4x4 unitary", vq_kak, cq_kak if HAS_CIRQ else None),
	"single_qubit_decompose": (
		"Single-qubit matrix to PhasedXZ",
		vq_single_qubit_decompose,
		cq_single_qubit_decompose if HAS_CIRQ else None,
	),
	"pauli_multiply": (
		f"Multiply Pauli strings {PAULI_MULTIPLY_COUNT} times",
		vq_pauli_multiply,
		cq_pauli_multiply if HAS_CIRQ else None,
	),
	"qasm": (f"OpenQASM export of a {QASM_GATES}-gate circuit", vq_qasm, cq_qasm if HAS_CIRQ else None),
	"bloch_vector": (f"Bloch vector of {SIM_QUBITS} qubits", vq_bloch_vector, cq_bloch_vector if HAS_CIRQ else None),
}


def _rss_kb() -> int:
	try:
		with open("/proc/self/statm") as handle:
			pages = int(handle.read().split()[1])
		return pages * (os.sysconf("SC_PAGE_SIZE") // 1024)
	except OSError:  # pragma: no cover - non-Linux
		import resource

		return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def _peak_rss_kb() -> int:
	try:
		with open("/proc/self/status") as handle:
			for line in handle:
				if line.startswith("VmHWM:"):
					return int(line.split()[1])
	except OSError:  # pragma: no cover - non-Linux
		pass
	return _rss_kb()


def _run_one(name: str, package: str) -> dict[str, float]:
	task = TASKS[name][1] if package == "vqubit" else TASKS[name][2]
	assert task is not None
	task()  # warm up imports and caches
	gc.collect()
	baseline = _rss_kb()
	best = float("inf")
	for _ in range(REPEATS):
		start = time.perf_counter()
		task()
		elapsed = time.perf_counter() - start
		best = min(best, elapsed)
	return {"time_s": best, "baseline_kb": baseline, "peak_kb": _peak_rss_kb()}


def _measure(name: str, package: str) -> dict[str, float]:
	"""Run one benchmark in a fresh subprocess and return its time and memory."""
	result = subprocess.run(
		[sys.executable, os.path.abspath(__file__), "--one", name, package],
		capture_output=True,
		text=True,
		check=True,
	)
	payload = result.stdout.strip().splitlines()[-1]
	return json.loads(payload)


def main() -> None:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--one", nargs=2, metavar=("NAME", "PACKAGE"))
	parser.add_argument("--json", action="store_true", help="print raw JSON instead of a table")
	args = parser.parse_args()

	if args.one:
		print(json.dumps(_run_one(args.one[0], args.one[1])))
		return

	rows = []
	for name, (description, _vq, cq) in TASKS.items():
		vq_result = _measure(name, "vqubit")
		cq_result = _measure(name, "cirq") if cq is not None else None
		rows.append((name, description, vq_result, cq_result))

	if args.json:
		print(json.dumps([{"task": n, "description": d, "vqubit": v, "cirq": c} for n, d, v, c in rows], indent=2))
		return

	header = f"{'task':<24}{'vqubit time':>14}{'cirq time':>14}{'speedup':>10}{'vqubit RAM':>14}{'cirq RAM':>14}"
	print(header)
	print("-" * len(header))
	for name, _description, v, c in rows:
		if c is None:
			print(
				f"{name:<24}{v['time_s'] * 1000:>12.1f}ms{'n/a':>14}{'n/a':>10}{v['peak_kb'] / 1024:>12.1f}MB{'n/a':>14}"
			)
			continue
		speedup = c["time_s"] / v["time_s"]
		print(
			f"{name:<24}{v['time_s'] * 1000:>12.1f}ms{c['time_s'] * 1000:>12.1f}ms{speedup:>9.2f}x"
			f"{v['peak_kb'] / 1024:>12.1f}MB{c['peak_kb'] / 1024:>12.1f}MB"
		)


if __name__ == "__main__":
	main()
