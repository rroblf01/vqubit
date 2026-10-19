"""Self-tests for the Cirq parity harness in ``benchmarks/parity.py``.

The harness itself runs in CI as a script (the ``cirq-parity`` job). These tests
cover the parts that are easy to get wrong: that the contract currently holds,
that a documented divergence cannot be added without a reason, and that an
allowance which no longer silences anything is reported as stale. They are
skipped when Cirq is not installed.
"""

from __future__ import annotations

import importlib.util
import pathlib
from typing import Any

import pytest

HARNESS = pathlib.Path(__file__).resolve().parent.parent / "benchmarks" / "parity.py"


@pytest.fixture(scope="module")
def parity() -> Any:
	"""The harness module, loaded by path (``benchmarks/`` is not a package)."""
	pytest.importorskip("cirq")
	spec = importlib.util.spec_from_file_location("vqubit_parity_harness", HARNESS)
	assert spec is not None and spec.loader is not None
	module = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(module)
	return module


def test_the_contract_holds(parity: Any) -> None:
	report = parity.build_report()
	unexpected = report.failures()
	stale = report.stale()
	assert unexpected == [], f"undocumented divergences: {unexpected}"
	assert stale == [], f"allowances that silence nothing: {[t.key for t in stale]}"


def test_every_section_reports(parity: Any) -> None:
	assert set(parity.build_report().sections) == set(parity.SECTIONS)


def test_the_diagram_section_covers_the_plain_gates(parity: Any) -> None:
	"""``vq.CZ`` and its fractional power must be drawn against Cirq's."""
	drawn = {name: ok for name, ok, _detail in parity.build_report().sections["diagram"]}
	for name in ("CZ", "CNOT", "SWAP", "ISWAP", "CCX", "CZ**0.5", "SWAP**0.5", "ISWAP**0.5"):
		assert drawn[name], f"the diagram of {name} is not compared with Cirq's"


def test_an_allowance_without_a_reason_is_refused(parity: Any, tmp_path: Any) -> None:
	path = tmp_path / "parity_tolerated.toml"
	path.write_text('[[tolerated]]\ncheck = "power"\nname = "CZ"\nreason = "   "\n')
	with pytest.raises(SystemExit):
		parity.load_tolerated(str(path))


def test_an_incomplete_allowance_is_refused(parity: Any, tmp_path: Any) -> None:
	path = tmp_path / "parity_tolerated.toml"
	path.write_text('[[tolerated]]\ncheck = "power"\nreason = "because"\n')
	with pytest.raises(SystemExit):
		parity.load_tolerated(str(path))


def test_an_allowance_that_silences_nothing_is_stale(parity: Any) -> None:
	item = parity.Tolerance(check="power", name="NOPE", reason="invented")
	report = parity.Report({item.key: item})
	report.check("power", "CZ", True)
	assert report.stale() == [item]
	assert report.failures() == []


def test_a_glob_allowance_covers_its_family(parity: Any) -> None:
	item = parity.Tolerance(check="isinstance", name="*-EigenGate", reason="hierarchy")
	report = parity.Report({item.key: item})
	report.check("isinstance", "CZ<-EigenGate", False)
	report.check("isinstance", "SWAP<-SwapPowGate", False)
	assert report.stale() == []
	assert [name for _section, name, _detail in report.failures()] == ["SWAP<-SwapPowGate"]
