"""The stubtest allowlist is a short, pinned list of vcraft gaps.

`mypy.stubtest` compares the stub in the wheel with the extension itself. Six
things it reports come from the code generator or from the C API and cannot be
fixed here (`UPSTREAM.md` §1 and §2), so they are allowlisted — and this test
pins the list, so adding an entry to silence a *new* error takes a deliberate
change here as well.
"""

from __future__ import annotations

from pathlib import Path

ALLOWLIST = Path(__file__).resolve().parent.parent / "tools" / "stubtest_allowlist.txt"

#: Exactly what may be silenced: the constructors, which the extension
#: implements in `tp_new` (so the generated `__init__` has no counterpart at run
#: time), and the three types, which mypy 1.19 wants marked `@disjoint_base` (a
#: PEP 800 decorator that does not exist on Python 3.11, where the stub also has
#: to work).
EXPECTED = {
	"vqubit_native.NativeState.__init__",
	"vqubit_native.NativeDensity.__init__",
	"vqubit_native.NativeTableau.__init__",
	"vqubit_native.NativeState",
	"vqubit_native.NativeDensity",
	"vqubit_native.NativeTableau",
}


def _entries() -> list[str]:
	return [
		line.strip()
		for line in ALLOWLIST.read_text(encoding="utf-8").splitlines()
		if line.strip() and not line.strip().startswith("#")
	]


def test_the_allowlist_is_exactly_the_documented_gaps() -> None:
	assert sorted(_entries()) == sorted(EXPECTED)


def test_no_entry_is_repeated() -> None:
	entries = _entries()
	assert len(entries) == len(set(entries))


def test_the_allowlist_explains_itself() -> None:
	text = ALLOWLIST.read_text(encoding="utf-8")
	assert "UPSTREAM.md" in text
	assert "Any new" in text or "Anything new" in text
