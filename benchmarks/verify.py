"""Cross-validate vqubit's numerical results against Cirq.

Run ``python benchmarks/verify.py`` with Cirq installed (``uv pip install
--group bench``). The script compares the *results* of both libraries, not their
speed: gate unitaries, channels, circuit simulation, sampling, decompositions,
Pauli/Clifford algebra and the public API surface. It exits with a non-zero
status if any check fails, so it can gate a CI job.

Gate unitaries are compared up to a global phase (the usual notion of gate
equivalence); global-phase-only gates are compared exactly. Everything else is
compared with a tight absolute tolerance.
"""

from __future__ import annotations

import argparse
import os
import sys
import types
from collections.abc import Callable
from typing import Any

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

try:
	import cirq
except ImportError:  # pragma: no cover - exercised only without the bench group
	print("Cirq is required for the cross-validation harness.", file=sys.stderr)
	print("Install it with: uv pip install --group bench", file=sys.stderr)
	raise SystemExit(2) from None

import vqubit as vq

ATOL = 1e-7
SAMPLING_TV = 0.06


class Report:
	"""Collect pass/fail checks grouped by section."""

	def __init__(self) -> None:
		self.sections: dict[str, list[tuple[str, bool, str]]] = {}

	def check(self, section: str, name: str, ok: bool, detail: str = "") -> None:
		self.sections.setdefault(section, []).append((name, bool(ok), detail))

	def failures(self) -> int:
		return sum(1 for checks in self.sections.values() for _n, ok, _d in checks if not ok)

	def total(self) -> int:
		return sum(len(checks) for checks in self.sections.values())

	def print(self) -> None:
		for section, checks in self.sections.items():
			failed = [c for c in checks if not c[1]]
			status = "OK" if not failed else f"FAIL ({len(failed)}/{len(checks)})"
			print(f"[{section:>16}] {status}")
			for name, ok, detail in checks:
				if not ok:
					print(f"    - {name}{(' :: ' + detail) if detail else ''}")


def up_to_phase(a: Any, b: Any, atol: float = ATOL) -> bool:
	"""Whether two matrices/vectors are equal up to a global phase."""
	a = np.asarray(a, dtype=complex)
	b = np.asarray(b, dtype=complex)
	if a.shape != b.shape:
		return False
	if a.size == 0:
		return True
	if a.size == 1:
		return bool(np.allclose(a, b, atol=atol))
	index = np.unravel_index(np.argmax(np.abs(b)), b.shape)
	if abs(b[index]) < 1e-12:
		return bool(np.allclose(a, b, atol=atol))
	phase = a[index] / b[index]
	return bool(np.allclose(a, phase * b, atol=atol))


def choi(kraus: Any) -> np.ndarray:
	"""Canonical Choi matrix of a channel given its Kraus operators."""
	ops = [np.asarray(k, dtype=complex) for k in kraus]
	dim = ops[0].shape[0]
	out = np.zeros((dim * dim, dim * dim), dtype=complex)
	for op in ops:
		flat = op.reshape(-1)
		out += np.outer(flat, np.conj(flat))
	return out


# --- gate unitaries ----------------------------------------------------------


def _gate_factories() -> dict[str, Callable[[Any], Any]]:
	V = vq

	def both(name: str) -> Callable[[Any], Any]:
		return lambda m: getattr(m, name)

	return {
		"X": both("X"),
		"Y": both("Y"),
		"Z": both("Z"),
		"H": both("H"),
		"S": both("S"),
		"T": both("T"),
		"I": both("I"),
		"XPow": lambda m: m.XPowGate(exponent=0.3),
		"YPow": lambda m: m.YPowGate(exponent=0.3),
		"ZPow": lambda m: m.ZPowGate(exponent=0.3),
		"HPow": lambda m: m.HPowGate(exponent=0.3),
		"CXPow": lambda m: m.CXPowGate(exponent=0.3),
		"CZPow": lambda m: m.CZPowGate(exponent=0.3),
		"CCXPow": lambda m: m.CCXPowGate(exponent=0.3),
		"CCZPow": lambda m: m.CCZPowGate(exponent=0.3),
		"XXPow": lambda m: m.XXPowGate(exponent=0.3),
		"YYPow": lambda m: m.YYPowGate(exponent=0.3),
		"ZZPow": lambda m: m.ZZPowGate(exponent=0.3),
		"SwapPow": lambda m: m.SwapPowGate(exponent=0.3),
		"ISwapPow": lambda m: m.ISwapPowGate(exponent=0.3),
		"PhasedISwapPow": lambda m: m.PhasedISwapPowGate(exponent=0.3, phase_exponent=0.25),
		"PhasedXPow": lambda m: m.PhasedXPowGate(phase_exponent=0.2, exponent=0.3),
		"PhasedXZ": lambda m: m.PhasedXZGate(x_exponent=0.2, z_exponent=0.3, axis_phase_exponent=0.4),
		"FSim": lambda m: m.FSimGate(0.3, 0.4),
		"PhasedFSim": lambda m: m.PhasedFSimGate(0.3, 0.1, 0.2, 0.3, 0.4),
		"MS": lambda m: m.MSGate(rads=0.3) if m is V else m.ms(0.3),
		"Rx": lambda m: m.rx(0.3),
		"Ry": lambda m: m.ry(0.3),
		"Rz": lambda m: m.rz(0.3),
		"cphase": lambda m: m.cphase(0.3),
		"givens": lambda m: m.givens(0.3),
		"MatrixGate": lambda m: m.MatrixGate(np.array([[0, 1], [1, 0]])),
		"QubitPermutation": lambda m: m.QubitPermutationGate([1, 0, 3, 2]),
		"GlobalPhase": lambda m: m.global_phase(0.3) if m is V else m.global_phase_operation(np.exp(1j * np.pi * 0.3)),
		"PhaseGradient": lambda m: m.PhaseGradientGate(num_qubits=2, exponent=0.5),
		"QFT2": lambda m: m.QuantumFourierTransformGate(2),
		"QFT2nr": lambda m: m.QuantumFourierTransformGate(2, without_reverse=True),
		"QFT3": lambda m: m.QuantumFourierTransformGate(3),
		"QFT3nr": lambda m: m.QuantumFourierTransformGate(3, without_reverse=True),
		"Diagonal": lambda m: m.DiagonalGate([0.5, 1.0, 1.5, 2.0]),
		"TwoQubitDiagonal": lambda m: m.TwoQubitDiagonalGate([0.5, 1.0, 1.5, 2.0]),
		"ThreeQubitDiagonal": lambda m: m.ThreeQubitDiagonalGate([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]),
	}


def check_gates(report: Report) -> None:
	for name, factory in _gate_factories().items():
		try:
			vgate = factory(vq)
		except Exception as exc:
			report.check("gates", name, False, f"vqubit: {type(exc).__name__}: {exc}")
			continue
		try:
			cgate = factory(cirq)
		except Exception as exc:
			report.check("gates", name, False, f"cirq: {type(exc).__name__}: {exc}")
			continue
		try:
			got = vq.unitary(vgate)
			expected = cirq.unitary(cgate)
		except Exception as exc:
			report.check("gates", name, False, f"unitary: {type(exc).__name__}: {exc}")
			continue
		report.check("gates", name, up_to_phase(got, expected))


# --- channels ----------------------------------------------------------------


def _channel_factories() -> dict[str, Callable[[Any], Any]]:
	return {
		"depolarize": lambda m: m.depolarize(0.2),
		"amplitude_damp": lambda m: m.amplitude_damp(0.2),
		"phase_damp": lambda m: m.phase_damp(0.2),
		"bit_flip": lambda m: m.bit_flip(0.2),
		"phase_flip": lambda m: m.phase_flip(0.2),
		"generalized_amplitude_damp": lambda m: m.generalized_amplitude_damp(0.2, 0.3),
		"asymmetric_depolarize": lambda m: m.asymmetric_depolarize(0.1, 0.2, 0.3),
		"reset": lambda m: m.reset(m.LineQubit(0)).gate,
		"kraus": lambda m: m.KrausChannel(
			[np.array([[1, 0], [0, np.sqrt(0.8)]]), np.array([[0, np.sqrt(0.2)], [0, 0]])]
		),
	}


def check_channels(report: Report) -> None:
	for name, factory in _channel_factories().items():
		try:
			got = choi(vq.kraus(factory(vq)))
			expected = choi(cirq.kraus(factory(cirq)))
		except Exception as exc:
			report.check("channels", name, False, f"{type(exc).__name__}: {exc}")
			continue
		report.check("channels", name, np.allclose(got, expected, atol=ATOL))


# --- circuits ----------------------------------------------------------------


def _unitary_circuits() -> dict[str, Callable[[Any], Any]]:
	def ghz(m: Any) -> Any:
		q = [m.LineQubit(i) for i in range(3)]
		return m.Circuit(m.H(q[0]), m.CNOT(q[0], q[1]), m.CNOT(q[1], q[2]), m.rx(0.3)(q[0]), m.rz(0.7)(q[2]))

	def dense(m: Any) -> Any:
		q = [m.LineQubit(i) for i in range(4)]
		return m.Circuit(
			[m.H(x) for x in q],
			[m.CZ(q[i], q[i + 1]) for i in range(3)],
			[m.X(q[0]), m.Y(q[1]), m.S(q[2]), m.T(q[3])],
		)

	return {"ghz": ghz, "dense4": dense}


def _noisy_circuits() -> dict[str, Callable[[Any], Any]]:
	def ghz_noisy(m: Any) -> Any:
		q = [m.LineQubit(i) for i in range(3)]
		return m.Circuit(
			m.H(q[0]),
			m.CNOT(q[0], q[1]),
			m.CNOT(q[1], q[2]),
			m.depolarize(0.05)(q[1]),
			m.amplitude_damp(0.1)(q[2]),
		)

	return {"ghz_noisy": ghz_noisy}


def check_circuits(report: Report) -> None:
	for name, factory in _unitary_circuits().items():
		vcircuit = factory(vq)
		ccircuit = factory(cirq)

		report.check("circuits", f"{name}.unitary", up_to_phase(vq.unitary(vcircuit), cirq.unitary(ccircuit)))

		vstate = vq.Simulator().simulate(vcircuit).state_vector()
		cstate = cirq.Simulator().simulate(ccircuit).state_vector()
		report.check("circuits", f"{name}.state_vector", up_to_phase(vstate, cstate))

		vdm = np.asarray(vq.DensityMatrixSimulator().simulate(vcircuit).density_matrix())
		cdm = np.asarray(cirq.DensityMatrixSimulator().simulate(ccircuit).final_density_matrix)
		report.check("circuits", f"{name}.density_matrix", np.allclose(vdm, cdm, atol=ATOL))

	for name, factory in _noisy_circuits().items():
		vcircuit = factory(vq)
		ccircuit = factory(cirq)
		vdm = np.asarray(vq.DensityMatrixSimulator().simulate(vcircuit).density_matrix())
		cdm = np.asarray(cirq.DensityMatrixSimulator().simulate(ccircuit).final_density_matrix)
		report.check("circuits", f"{name}.density_matrix", np.allclose(vdm, cdm, atol=ATOL))


# --- sampling ----------------------------------------------------------------


def check_sampling(report: Report) -> None:
	def build(m: Any) -> Any:
		q = [m.LineQubit(i) for i in range(3)]
		return m.Circuit(
			m.H(q[0]),
			m.CNOT(q[0], q[1]),
			m.CNOT(q[1], q[2]),
			m.measure(*q, key="m"),
		)

	reps = 4000
	vhist = dict(vq.Simulator(seed=1).run(build(vq), repetitions=reps).histogram(key="m"))
	chist = dict(cirq.Simulator(seed=1).run(build(cirq), repetitions=reps).histogram(key="m"))
	keys = set(vhist) | set(chist)
	tv = 0.5 * sum(abs(vhist.get(k, 0) - chist.get(k, 0)) for k in keys) / reps
	report.check("sampling", "histogram_matches", tv <= SAMPLING_TV, f"total variation {tv:.3f}")

	# The GHZ circuit only ever produces |000> and |111>.
	support_ok = set(vhist) <= {0, 7}
	report.check("sampling", "histogram_support", support_ok, f"keys {sorted(vhist)}")


# --- decompositions ----------------------------------------------------------


def _random_unitary(dim: int, seed: int) -> np.ndarray:
	rng = np.random.default_rng(seed)
	matrix = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
	q, _r = np.linalg.qr(matrix)
	return q


def _apply(ops: Any, qubit: Any) -> list[Any]:
	"""Turn a sequence of gates or operations into operations on ``qubit``."""
	out = []
	for op in ops:
		try:
			out.append(op(qubit))
		except TypeError:
			out.append(op)
	return out


def check_decompositions(report: Report) -> None:
	u4 = _random_unitary(4, 0)
	u2 = _random_unitary(2, 1)

	vk = np.asarray(vq.kak_vector(u4))
	ck = np.asarray(cirq.kak_vector(u4))
	report.check("decompositions", "kak_vector", np.allclose(np.sort(np.abs(vk)), np.sort(np.abs(ck)), atol=ATOL))

	q0 = vq.LineQubit(0)
	vphxz = vq.unitary(vq.Circuit(_apply(vq.single_qubit_matrix_to_phased_x_z(u2), q0)))
	report.check("decompositions", "single_qubit_reconstructs", up_to_phase(vphxz, u2))
	cphxz = cirq.unitary(cirq.Circuit(_apply(cirq.single_qubit_matrix_to_phased_x_z(u2), cirq.LineQubit(0))))
	report.check("decompositions", "single_qubit_matches_cirq", up_to_phase(vphxz, cphxz))

	q0, q1 = vq.LineQubit(0), vq.LineQubit(1)
	vops = vq.two_qubit_matrix_to_cz_operations(q0, q1, u4)
	vcz = vq.unitary(vq.Circuit(vops))
	report.check("decompositions", "two_qubit_cz_reconstructs", up_to_phase(vcz, u4))


# --- pauli and clifford ------------------------------------------------------


def check_algebra(report: Report) -> None:
	vps = vq.PauliString({vq.LineQubit(0): vq.Pauli.X, vq.LineQubit(1): vq.Pauli.Z})
	cps = cirq.PauliString({cirq.LineQubit(0): cirq.X, cirq.LineQubit(1): cirq.Z})
	report.check("algebra", "pauli_string_matrix", up_to_phase(vps.matrix(), cps.matrix()))

	vprod = vq.PauliString({vq.LineQubit(0): vq.Pauli.X}) * vq.PauliString({vq.LineQubit(0): vq.Pauli.Y})
	cprod = cirq.PauliString({cirq.LineQubit(0): cirq.X}) * cirq.PauliString({cirq.LineQubit(0): cirq.Y})
	report.check("algebra", "pauli_string_product", up_to_phase(vprod.matrix(), cprod.matrix()))

	def clifford_circuit(m: Any) -> Any:
		q = [m.LineQubit(i) for i in range(3)]
		return m.Circuit(m.H(q[0]), m.CNOT(q[0], q[1]), m.CNOT(q[1], q[2]), m.S(q[2]))

	vres = vq.CliffordSimulator().simulate(clifford_circuit(vq))
	cres = cirq.CliffordSimulator().simulate(clifford_circuit(cirq))
	vstate = vres.state_vector() if hasattr(vres, "state_vector") else vres.final_state.state_vector()
	cstate = cres.state_vector() if hasattr(cres, "state_vector") else cres.final_state.state_vector()
	report.check("algebra", "clifford_state_vector", up_to_phase(vstate, cstate))

	ctab_api = {name for name in dir(cirq.CliffordTableau) if not name.startswith("_")}
	vtab_api = {name for name in dir(vq.CliffordTableau) if not name.startswith("_")}
	missing = sorted(ctab_api - vtab_api)
	report.check("algebra", "clifford_tableau_api", not missing, f"missing {missing}" if missing else "")


# --- public API surface ------------------------------------------------------


def check_api(report: Report) -> None:
	cirq_names = {
		name for name in dir(cirq) if not name.startswith("_") and not isinstance(getattr(cirq, name), types.ModuleType)
	}
	vq_names = {name for name in dir(vq) if not name.startswith("_")}
	missing = sorted(cirq_names - vq_names)
	report.check("api", "top_level_names", not missing, f"missing {missing}" if missing else "")

	import inspect

	signature_cases = [
		"two_qubit_matrix_to_cz_operations",
		"two_qubit_matrix_to_cz_isometry",
		"two_qubit_matrix_to_diagonal_and_cz_operations",
		"merge_operations",
		"merge_moments",
	]
	for name in signature_cases:
		vparams = list(inspect.signature(getattr(vq, name)).parameters)
		cparams = list(inspect.signature(getattr(cirq, name)).parameters)
		report.check(
			"api",
			f"signature.{name}",
			vparams[: len(cparams)] == cparams,
			f"vqubit {vparams} vs cirq {cparams}",
		)


def main() -> None:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--only", help="run only the named sections (comma separated)")
	args = parser.parse_args()

	sections: dict[str, Callable[[Report], None]] = {
		"gates": check_gates,
		"channels": check_channels,
		"circuits": check_circuits,
		"sampling": check_sampling,
		"decompositions": check_decompositions,
		"algebra": check_algebra,
		"api": check_api,
	}
	if args.only:
		wanted = {s.strip() for s in args.only.split(",")}
		sections = {k: v for k, v in sections.items() if k in wanted}

	report = Report()
	for check in sections.values():
		check(report)
	report.print()

	failures = report.failures()
	print(f"\n{report.total() - failures}/{report.total()} checks passed.")
	raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
	main()
