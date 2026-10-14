"""Run this after `vcraft develop` to check the extension is importable."""

from __future__ import annotations

import vqubit


def main() -> None:
	print(f"vqubit {vqubit.__version__}")


if __name__ == "__main__":
	main()
