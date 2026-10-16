"""Repair the source distribution that ``vcraft sdist`` builds.

``vcraft sdist`` (1.4.0) packs ``src/<module>.v`` and the Python package but
leaves out the module's sibling ``.v`` files, so a wheel rebuilt from the
published sdist exposes only ``version()``. This adds every ``src/**/*.v`` file
the project tracks that the archive is missing and rewrites it in place.
"""

from __future__ import annotations

import sys
import tarfile
from pathlib import Path

# `src/_vcraft_generated.v` is written by the build, never shipped.
_SKIP = {"_vcraft_generated.v"}
_TAR_SUFFIXES = (".tar.gz", ".tgz")


def _missing(root: Path, members: list[tarfile.TarInfo]) -> tuple[str, list[tuple[str, Path]]]:
	prefix = members[0].name.split("/", 1)[0] + "/"
	present = {member.name for member in members}
	missing: list[tuple[str, Path]] = []
	for source in sorted((root / "src").rglob("*.v")):
		if source.name in _SKIP:
			continue
		arcname = prefix + "src/" + source.relative_to(root / "src").as_posix()
		if arcname not in present:
			missing.append((arcname, source))
	return prefix, missing


def fix_sdist(path: Path, root: Path) -> bool:
	temporary = path.with_name(path.name + ".tmp")
	with tarfile.open(path) as archive:
		members = archive.getmembers()
		_, missing = _missing(root, members)
		if not missing:
			return False
		with tarfile.open(temporary, "w:gz") as repaired:
			for member in members:
				repaired.addfile(member, archive.extractfile(member))
			for arcname, source in missing:
				repaired.add(source, arcname=arcname)
	temporary.replace(path)
	return True


def main(argv: list[str]) -> int:
	root = Path.cwd()
	archives = [
		Path(argument) for argument in argv[1:] if argument.endswith(_TAR_SUFFIXES) and not argument.endswith(".tmp")
	]
	repaired = sum(1 for archive in archives if fix_sdist(archive, root))
	print(f"repaired {repaired} of {len(archives)} sdist(s)")
	return 0


if __name__ == "__main__":
	raise SystemExit(main(sys.argv))
