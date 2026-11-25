"""Locate the repository a command is acting on: the first directory at or above
the current one that holds data/data_host.yaml. `HERRING_REPO` overrides it."""
from __future__ import annotations

import os
from pathlib import Path

MARKER = Path("data") / "data_host.yaml"


def find_repo(start: Path | None = None) -> Path:
    env = os.environ.get("HERRING_REPO")
    if env:
        return Path(env).resolve()
    here = (start or Path.cwd()).resolve()
    for d in (here, *here.parents):
        if (d / MARKER).is_file():
            return d
    package_repo = Path(__file__).resolve().parents[2]
    if (package_repo / MARKER).is_file():
        return package_repo
    raise FileNotFoundError(
        f"No project-herring-copy repository at or above {here}; run from a checkout or set HERRING_REPO."
    )
