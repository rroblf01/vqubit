"""Check vqubit's public API contract against Cirq.

Where ``benchmarks/verify.py`` compares *numbers* (unitaries, channels,
decompositions), this script compares the *contract*, i.e. the promises a Cirq
user relies on when they swap the import:

* every public name Cirq exports also exists here;
* raising a gate to a power (``gate ** 0.5``) is possible exactly when Cirq
  allows it, and yields the same unitary;
* the ``isinstance`` lattice over the shared gate classes is the same, so
  ``isinstance(gate, vqubit.CZPowGate)`` means what it means in Cirq;
* gate classes that take no arguments build the same unitary.

Run ``python benchmarks/parity.py`` (or ``pytest tests/test_cirq_parity.py``)
with Cirq installed. A divergence is not necessarily a bug: some are deliberate.
Every deliberate one must be listed, with a reason, in
``benchmarks/parity_tolerated.toml``. The list cannot rot, because an allowance
that no longer silences a real divergence is itself reported as a failure.

Exit status: 0 when every expected divergence is tolerated, 1 when something is
unexpected (a new divergence, or a stale allowance), 2 when Cirq is missing.
"""

from __future__ import annotations

import argparse
import fnmatch
import importlib
import os
import sys
import tomllib
import types
from collections.abc import Callable
from typing import Any, NamedTuple

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

try:
	import cirq
except ImportError:  # pragma: no cover - exercised only without the bench group
	print("Cirq is required for the parity harness.", file=sys.stderr)
	print("Install it with: uv pip install --group bench", file=sys.stderr)
	raise SystemExit(2) from None

import vqubit as vq

ATOL = 1e-7
TOLERATED_PATH = os.path.join(os.path.dirname(__file__), "parity_tolerated.toml")


class Tolerance(NamedTuple):
	"""One deliberate divergence, as written in the TOML file.

	``name`` may be a shell glob (``fnmatch``), e.g. ``*-EigenGate``, so a whole
	class of structural divergences shares one reasoned entry.
	"""

	check: str
	name: str
	reason: str
	plan: str = "permanent"

	@property
	def key(self) -> str:
		return f"{self.check}:{self.name}"


def load_tolerated(path: str) -> dict[str, Tolerance]:
	"""Read the tolerated-divergence list, refusing entries without a reason."""
	with open(path, "rb") as handle:
		raw = tomllib.load(handle)
	tolerated: dict[str, Tolerance] = {}
	for entry in raw.get("tolerated", []):
		missing = {"check", "name", "reason"} - set(entry)
		if missing:
			raise SystemExit(f"{path}: an entry is missing {sorted(missing)}")
		reason = (entry["reason"] or "").strip()
		if not reason:
			raise SystemExit(f"{path}: {entry['check']}:{entry['name']} has an empty reason")
		item = Tolerance(check=entry["check"], name=entry["name"], reason=reason, plan=entry.get("plan", "permanent"))
		if item.key in tolerated:
			raise SystemExit(f"{path}: duplicate entry {item.key}")
		tolerated[item.key] = item
	return tolerated


class Report:
	"""Collect pass/fail checks grouped by section."""

	def __init__(self, tolerated: dict[str, Tolerance]) -> None:
		self.sections: dict[str, list[tuple[str, bool, str]]] = {}
		self.tolerated = tolerated
		self.used: set[str] = set()

	def tolerance(self, section: str, name: str) -> Tolerance | None:
		"""The allowance covering this check, if any (``name`` supports globs)."""
		exact = self.tolerated.get(f"{section}:{name}")
		if exact is not None:
			return exact
		for pattern, item in self.tolerated.items():
			if fnmatch.fnmatchcase(f"{section}:{name}", pattern):
				return item
		return None

	def check(self, section: str, name: str, ok: bool, detail: str = "") -> None:
		if not ok:
			item = self.tolerance(section, name)
			if item is not None:
				self.used.add(item.key)
		self.sections.setdefault(section, []).append((name, bool(ok), detail))

	def failures(self) -> list[tuple[str, str, str]]:
		"""Unexpected failures: (section, name, detail)."""
		out = []
		for section, checks in self.sections.items():
			for name, ok, detail in checks:
				if not ok and self.tolerance(section, name) is None:
					out.append((section, name, detail))
		return out

	def stale(self) -> list[Tolerance]:
		"""Allowances that silenced nothing this run."""
		return [t for key, t in self.tolerated.items() if key not in self.used]

	def print(self) -> None:
		for section, checks in self.sections.items():
			failed = [c for c in checks if not c[1] and self.tolerance(section, c[0]) is None]
			tolerated = [c for c in checks if not c[1] and self.tolerance(section, c[0]) is not None]
			label = f"{len(checks)} checks"
			if tolerated:
				label += f", {len(tolerated)} tolerated"
			status = "OK" if not failed else f"FAIL ({len(failed)})"
			print(f"[{section:>16}] {status} ({label})")
			for name, _ok, detail in failed:
				print(f"    ! {name}{(' :: ' + detail) if detail else ''}")

	def total(self) -> int:
		return sum(len(checks) for checks in self.sections.values())


# --- helpers -----------------------------------------------------------------


def up_to_phase(a: Any, b: Any, atol: float = ATOL) -> bool:
	"""Whether two matrices are equal up to a global phase."""
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


def public_names(module: types.ModuleType) -> set[str]:
	"""Public names of a module, submodules included.

	Submodules count: ``import vqubit.work`` mirroring ``import cirq.work`` is
	part of the contract (roadmap workstream 1.5).
	"""
	return {name for name in dir(module) if not name.startswith("_")}


def gate_instances(module: types.ModuleType) -> dict[str, Any]:
	"""Public names that are gate *instances* (``cirq.X``, ``cirq.CZ``)."""
	return {
		name: getattr(module, name) for name in public_names(module) if isinstance(getattr(module, name), module.Gate)
	}


def gate_classes(module: types.ModuleType) -> dict[str, type]:
	"""Public names that are gate *classes*."""
	out = {}
	for name in public_names(module):
		value = getattr(module, name)
		if isinstance(value, type) and issubclass(value, module.Gate):
			out[name] = value
	return out


def unitary(module: types.ModuleType, gate: Any) -> np.ndarray | None:
	try:
		return np.asarray(module.unitary(gate), dtype=complex)
	except Exception:
		return None


def power(gate: Any, exponent: float = 0.5) -> tuple[bool, Any]:
	"""Whether ``gate ** exponent`` works, and what it returns."""
	try:
		return True, gate**exponent
	except Exception:
		return False, None


def text_diagram(module: types.ModuleType, gate: Any) -> str | None:
	"""Draw ``gate`` on ``LineQubit``\\ s, or ``None`` when it cannot be drawn.

	Cirq 1.7 exposes ``num_qubits`` as a method on its Pauli singletons and as a
	property elsewhere, so both spellings are accepted.
	"""
	try:
		num_qubits = gate.num_qubits
		count = num_qubits() if callable(num_qubits) else int(num_qubits)
		if count > 4:
			return None
		qubits = [module.LineQubit(i) for i in range(count)]
		return module.Circuit(gate.on(*qubits)).to_text_diagram().rstrip()
	except Exception:
		return None


def same_value(gate: Any, power: int) -> bool | None:
	"""Whether ``gate ** power == gate``, or ``None`` when ``**`` is refused."""
	try:
		raised = gate**power
	except Exception:
		return None
	return bool(raised == gate)


# --- checks ------------------------------------------------------------------


def check_namespace(report: Report) -> None:
	"""Every public name of Cirq exists here, and mirrors submodule importability."""
	missing = sorted(public_names(cirq) - public_names(vq))
	for name in missing:
		report.check("namespace", name, False, "missing from vqubit")
	# A name being present in `dir()` does not prove the submodule imports.
	for name in sorted(public_names(cirq) & public_names(vq)):
		if not isinstance(getattr(cirq, name), types.ModuleType):
			continue
		try:
			importlib.import_module(f"vqubit.{name}")
		except Exception as exc:
			report.check("namespace", name, False, f"vqubit.{name} does not import: {type(exc).__name__}: {exc}")


def check_unitary(report: Report) -> None:
	shared = gate_instances(cirq).keys() & gate_instances(vq).keys()
	for name in sorted(shared):
		got = unitary(vq, getattr(vq, name))
		want = unitary(cirq, getattr(cirq, name))
		if want is None or got is None:
			report.check("unitary", name, got is want, f"cirq {'ok' if want is not None else 'raises'}")
			continue
		report.check("unitary", name, up_to_phase(got, want))


def check_power(report: Report) -> None:
	"""``**`` must be available exactly when Cirq allows it, and agree."""
	shared = gate_instances(cirq).keys() & gate_instances(vq).keys()
	for name in sorted(shared):
		cirq_ok, cirq_pow = power(getattr(cirq, name))
		vq_ok, vq_pow = power(getattr(vq, name))
		if not cirq_ok:
			report.check("power", name, not vq_ok, "vqubit accepts ** where Cirq refuses")
			continue
		if not vq_ok:
			report.check("power", name, False, f"Cirq accepts ** (-> {type(cirq_pow).__name__})")
			continue
		got = unitary(vq, vq_pow)
		want = unitary(cirq, cirq_pow)
		report.check(
			"power",
			name,
			got is not None and want is not None and up_to_phase(got, want),
			f"vqubit {type(vq_pow).__name__} vs cirq {type(cirq_pow).__name__}",
		)


def check_rotations(report: Report) -> None:
	"""The ``rx``/``ry``/``rz`` factories must scale their angle under ``**``."""
	for name in ("rx", "ry", "rz"):
		fv, fc = getattr(vq, name, None), getattr(cirq, name, None)
		if fv is None or fc is None:
			report.check("rotation", name, False, "missing factory")
			continue
		cirq_ok, cirq_pow = power(fc(0.3))
		vq_ok, vq_pow = power(fv(0.3))
		if not cirq_ok:
			report.check("rotation", name, not vq_ok)
			continue
		if not vq_ok:
			report.check("rotation", name, False, f"Cirq accepts ** (-> {cirq_pow!r})")
			continue
		got, want = unitary(vq, vq_pow), unitary(cirq, cirq_pow)
		report.check(
			"rotation",
			name,
			got is not None and want is not None and up_to_phase(got, want),
			f"vqubit {type(vq_pow).__name__} vs cirq {type(cirq_pow).__name__}",
		)


def check_isinstance(report: Report) -> None:
	"""The ``isinstance`` lattice over the shared gate classes must agree."""
	shared_instances = gate_instances(cirq).keys() & gate_instances(vq).keys()
	shared_classes = gate_classes(cirq).keys() & gate_classes(vq).keys()
	for name in sorted(shared_instances):
		cirq_gate, vq_gate = getattr(cirq, name), getattr(vq, name)
		for cls in sorted(shared_classes):
			want = isinstance(cirq_gate, getattr(cirq, cls))
			got = isinstance(vq_gate, getattr(vq, cls))
			report.check("isinstance", f"{name}<-{cls}", got == want, f"cirq {want} vs vqubit {got}")


def check_classes(report: Report) -> None:
	"""Gate classes buildable with no arguments must build the same unitary."""
	shared = gate_classes(cirq).keys() & gate_classes(vq).keys()
	for name in sorted(shared):
		try:
			cirq_gate = getattr(cirq, name)()
		except Exception:
			continue
		try:
			vq_gate = getattr(vq, name)()
		except Exception as exc:
			report.check("class", name, False, f"vqubit raises {type(exc).__name__}: {exc}")
			continue
		want, got = unitary(cirq, cirq_gate), unitary(vq, vq_gate)
		if want is None or got is None:
			report.check("class", name, got is want)
			continue
		report.check("class", name, up_to_phase(got, want))


def check_diagram(report: Report) -> None:
	"""The text diagram of a shared gate must read the same as Cirq's.

	Each shared gate instance is drawn twice, once as it is and once raised to
	``** 0.5``, which is how the power gates advertise their exponent
	(``@^0.5``). This is the check that sees the plain gates fold into the power
	gates: ``vq.CZ`` and ``cirq.CZ`` are both a ``CZPowGate`` at exponent 1 and
	both draw a bare ``@``.
	"""
	shared = gate_instances(cirq).keys() & gate_instances(vq).keys()
	for name in sorted(shared):
		cirq_gate, vq_gate = getattr(cirq, name), getattr(vq, name)
		want, got = text_diagram(cirq, cirq_gate), text_diagram(vq, vq_gate)
		report.check("diagram", name, got is not None and got == want, f"cirq {want!r} vs vqubit {got!r}")
		cirq_ok, cirq_pow = power(cirq_gate)
		vq_ok, vq_pow = power(vq_gate)
		if cirq_ok != vq_ok:
			report.check("diagram", f"{name}**0.5", False, "vqubit and Cirq disagree on whether ** is allowed")
		elif not cirq_ok:
			report.check("diagram", f"{name}**0.5", True)
		else:
			want, got = text_diagram(cirq, cirq_pow), text_diagram(vq, vq_pow)
			report.check("diagram", f"{name}**0.5", got is not None and got == want, f"cirq {want!r} vs vqubit {got!r}")


#: Gates whose value equality is a different contract entirely -- the Clifford
#: gates compare their tableau, and vqubit compares their matrices -- so the
#: ``equality`` section skips them instead of pretending to check them. Recorded
#: in ``ROADMAP.md``: until it is closed, ``vq.CZSWAP ** 3 != vq.CZSWAP`` where
#: Cirq says they are equal.
EQUALITY_EXEMPT = frozenset({"CZSWAP", "CXSWAP"})

#: Powers to use for the equality check: two whole turns, one more (which the
#: period folds back) and the inverse (which does not fold for every gate).
EQUALITY_POWERS = (2, 3, -1)


def check_equality(report: Report) -> None:
	"""``gate ** power == gate`` must answer as Cirq answers.

	A power gate *is* its eigenphases, and Cirq compares them modulo their period
	(through ``PeriodicValue``), which is why ``cirq.X ** 3 == cirq.X`` and
	``cirq.CNOT ** -1 == cirq.CNOT`` while ``cirq.ISWAP ** -1 != cirq.ISWAP``.
	When neither library defines the power the check passes: whether ``**`` is
	available at all is the ``power`` section's business.
	"""
	shared = gate_instances(cirq).keys() & gate_instances(vq).keys()
	for name in sorted(shared - EQUALITY_EXEMPT):
		for power in EQUALITY_POWERS:
			want = same_value(getattr(cirq, name), power)
			got = same_value(getattr(vq, name), power)
			if want is None or got is None:
				report.check(
					"equality", f"{name}**{power}", want is got, "one library defines the power, the other refuses it"
				)
				continue
			report.check("equality", f"{name}**{power}", got == want, f"cirq {want} vs vqubit {got}")


SECTIONS: dict[str, Callable[[Report], None]] = {
	"namespace": check_namespace,
	"unitary": check_unitary,
	"power": check_power,
	"rotation": check_rotations,
	"isinstance": check_isinstance,
	"class": check_classes,
	"diagram": check_diagram,
	"equality": check_equality,
}


def build_report(only: str | None = None) -> Report:
	"""Run the requested checks and return the report (used by the tests too)."""
	tolerated = load_tolerated(TOLERATED_PATH)
	report = Report(tolerated)
	sections = SECTIONS
	if only:
		wanted = {s.strip() for s in only.split(",")}
		sections = {k: v for k, v in SECTIONS.items() if k in wanted}
	for check in sections.values():
		check(report)
	return report


def main() -> None:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--only", help="run only the named sections (comma separated)")
	args = parser.parse_args()

	report = build_report(args.only)
	report.print()

	failures = report.failures()
	stale = report.stale()
	for section, name, detail in failures:
		print(f"UNEXPECTED [{section}] {name}{(' :: ' + detail) if detail else ''}")
	for item in stale:
		print(f"STALE      [{item.check}] {item.name} (plan {item.plan}) :: no longer failing, drop it")

	tolerated_count = len(report.used)
	print(
		f"\n{report.total()} checks, {tolerated_count} tolerated divergences, {len(failures)} unexpected, {len(stale)} stale."
	)
	raise SystemExit(1 if failures or stale else 0)


if __name__ == "__main__":
	main()
