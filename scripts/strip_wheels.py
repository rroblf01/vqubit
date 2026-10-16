"""Strip the native library inside built wheels.

``vcraft build --strip`` is a no-op in vcraft 1.4.0 (the flag parses and
``vcraft info`` reports it, but the extension in the wheel keeps its symbol
table), so the release workflow runs this after the build.

For every wheel given on the command line it extracts the extension, runs the
platform ``strip`` on it, and rewrites the wheel with a recomputed ``RECORD``.
It is best effort: if ``strip`` is missing or fails, the wheel is left exactly
as it was, so a platform without a usable ``strip`` still publishes.
"""

from __future__ import annotations

import base64
import hashlib
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

# The extension module of the package: `.so` on Linux and macOS (CPython uses
# `.so` for Mach-O extension modules too) and `.pyd` on Windows.
_NATIVE_SUFFIXES = (".so", ".pyd")
_RECORD_SUFFIX = ".dist-info/RECORD"


def _strip(data: bytes) -> bytes | None:
	strip = shutil.which("strip")
	if strip is None:
		print("strip: not found, wheels left untouched", file=sys.stderr)
		return None
	with tempfile.TemporaryDirectory() as tmp:
		lib = Path(tmp) / "lib.bin"
		lib.write_bytes(data)
		result = subprocess.run([strip, str(lib)], capture_output=True)
		if result.returncode != 0:
			message = result.stderr.decode(errors="replace").strip()
			print(f"strip: failed ({message or result.returncode})", file=sys.stderr)
			return None
		return lib.read_bytes()


def _record_line(name: str, data: bytes | None) -> str:
	if data is None:
		return f"{name},,"
	digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).decode().rstrip("=")
	return f"{name},sha256={digest},{len(data)}"


def strip_wheel(path: Path) -> bool:
	with zipfile.ZipFile(path) as archive:
		names = [item.filename for item in archive.infolist()]
		contents = {name: archive.read(name) for name in names}

	native = [name for name in names if name.endswith(_NATIVE_SUFFIXES)]
	if len(native) != 1:
		print(f"{path.name}: expected one native library, found {native}", file=sys.stderr)
		return False

	target = native[0]
	stripped = _strip(contents[target])
	if stripped is None or stripped == contents[target]:
		return False

	contents[target] = stripped
	record = next((name for name in names if name.endswith(_RECORD_SUFFIX)), None)
	if record is not None:
		lines = [_record_line(name, None if name == record else contents[name]) for name in names]
		contents[record] = ("\n".join(lines) + "\n").encode()

	temporary = path.with_name(path.name + ".tmp")
	with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as archive:
		for name in names:
			archive.writestr(name, contents[name])
	temporary.replace(path)
	return True


def main(argv: list[str]) -> int:
	wheels = [Path(argument) for argument in argv[1:]]
	stripped = sum(1 for wheel in wheels if strip_wheel(wheel))
	print(f"stripped {stripped} of {len(wheels)} wheel(s)")
	return 0


if __name__ == "__main__":
	raise SystemExit(main(sys.argv))
