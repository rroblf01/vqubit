"""mkdocs hook: substitute the project version into the pages.

The version lives in ``vcraft.toml`` (and is re-exported as
``vqubit.__version__``). Reading it here, at build time, keeps the published
documentation from showing a stale number when it is built for a release.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

_ROOT = Path(__file__).resolve().parent.parent
_VERSION_RE = re.compile(r'^version\s*=\s*"([^"]+)"', re.MULTILINE)


def _version() -> str:
	text = (_ROOT / "vcraft.toml").read_text(encoding="utf-8")
	match = _VERSION_RE.search(text)
	return match.group(1) if match else "0.0.0"


def on_page_markdown(markdown: str, **kwargs: Any) -> str:
	del kwargs
	return markdown.replace("{{ vqubit_version }}", _version())
