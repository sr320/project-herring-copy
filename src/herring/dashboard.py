"""Build site/data.json, the data behind the GitHub Pages dashboard (site/index.html).

Modelled on project-white-whale's dashboard. Nothing on the page is typed in;
it is assembled from:

* data/mirror_projects.yaml                the BioProjects in scope, in mirror order
* data/phases.yaml                         where each phase of PLAN.md stands
* metadata/ena-runs-*.tsv                  the committed ENA inventory: species, country
                                           and collection date of every run
* ENA's read-run report (live)             FASTQ files and bytes each project will need
* raw/<ACCESSION>/mirror_manifest.tsv      what is mirrored on gannet (data_host.yaml url_base)

The Pages workflow (.github/workflows/pages.yml) runs `herring build-dashboard`
on every push to main and every few hours, so mirror progress shows up without
a commit.
"""
from __future__ import annotations

import csv
import json
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import yaml

from herring.mirror import ena_md5, load_mirror_config

INVENTORY_GLOB = "ena-runs-*.tsv"
SPECIES = {"Clupea harengus": "Atlantic", "Clupea pallasii": "Pacific"}
EARLIEST_YEAR = 1970   # dates before this are counted, not plotted (one run says 1905)


def load_phases(repo: Path) -> list[dict]:
    path = Path(repo) / "data" / "phases.yaml"
    return (yaml.safe_load(path.read_text()) or {}).get("phases", []) if path.exists() else []


def load_inventory(repo: Path) -> list[dict]:
    """Every run in the newest committed inventory of each species (files sort by date)."""
    newest: dict[str, Path] = {}
    for p in sorted((Path(repo) / "metadata").glob(INVENTORY_GLOB)):
        newest[p.name.rsplit("-", 3)[0]] = p
    rows = []
    for p in newest.values():
        with p.open(newline="") as fh:
            rows += list(csv.DictReader(fh, delimiter="\t"))
    return rows


def _in_scope(row: dict, cfg: dict) -> bool:
    return (row.get("library_strategy") in cfg["library_strategies"]
            and row.get("instrument_platform") in cfg["platforms"]
            and row.get("library_source") in cfg["library_sources"]
            and row.get("scientific_name") in SPECIES)


def _year(date: str) -> int | None:
    try:
        return int((date or "")[:4])
    except ValueError:
        return None


def project_inventory(rows: list[dict], cfg: dict) -> dict[str, dict]:
    """Per BioProject: runs per species, countries, collection years, bases (in-scope runs only)."""
    accs = {p["accession"] for p in cfg["projects"]}
    out: dict[str, dict] = {}
    for r in rows:
        acc = r.get("study_accession")
        if acc not in accs or not _in_scope(r, cfg):
            continue
        p = out.setdefault(acc, {"species": Counter(), "countries": Counter(), "years": Counter(),
                                 "undated": 0, "bases": 0})
        p["species"][SPECIES[r["scientific_name"]]] += 1
        country = (r.get("country") or "").split(":")[0].strip()
        if country:
            p["countries"][country] += 1
        y = _year(r.get("collection_date"))
        if y is None:
            p["undated"] += 1
        else:
            p["years"][y] += 1
        p["bases"] += int(r.get("base_count") or 0)
    return out


def sampling_through_time(rows: list[dict], cfg: dict) -> dict:
    """Runs per collection year and species, over every in-scope project."""
    accs = {p["accession"] for p in cfg["projects"]}
    counts = {sp: Counter() for sp in SPECIES.values()}
    early = Counter()
    undated = Counter()
    for r in rows:
        if r.get("study_accession") not in accs or not _in_scope(r, cfg):
            continue
        sp, y = SPECIES[r["scientific_name"]], _year(r.get("collection_date"))
        if y is None:
            undated[sp] += 1
        elif y < EARLIEST_YEAR:
            early[sp] += 1
        else:
            counts[sp][y] += 1
    years = sorted({y for c in counts.values() for y in c})
    span = list(range(years[0], years[-1] + 1)) if years else []
    return {"years": span,
            "series": [{"species": sp, "counts": [counts[sp].get(y, 0) for y in span],
                        "early": early[sp], "undated": undated[sp]} for sp in SPECIES.values()]}


def mirror_summary(files: list[dict] | None, ledger: dict[str, dict] | None) -> dict | None:
    """What of a project's wanted files is on gannet against ENA's current MD5. None if unknown."""
    if ledger is None:
        return None
    if files is None:
        rows = list(ledger.values())
    else:
        rows = [ledger[f["file"]] for f in files
                if f["file"] in ledger and ena_md5(ledger[f["file"]]) == f["md5"]]
    sources = Counter(r.get("source") or "" for r in rows)
    return {"files": len(rows),
            "gb": round(sum(int(r.get("bytes") or 0) for r in rows) / 1e9, 1),
            "runs": len({r.get("run_accession") for r in rows}),
            "from_ncbi": sources.get("ncbi_sra", 0),
            "from_ena": len(rows) - sources.get("ncbi_sra", 0),
            "last_mirrored": max((r.get("mirrored_on") or "" for r in rows), default="")}


def _git(repo: Path, *args: str) -> str:
    proc = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True)
    return proc.stdout.strip() if proc.returncode == 0 else ""


def build(repo: Path, *, wanted: dict[str, list[dict] | None], ledgers: dict[str, dict] | None,
          repo_url: str, url_base: str | None, now: datetime | None = None) -> dict:
    """The dashboard's data.

    `wanted`: per accession, the files the mirror wants (from ENA), or None if ENA
    could not be read. `ledgers`: per accession, the gannet ledger ({} when there is
    none yet); None for the whole dict when gannet could not be read.
    """
    repo = Path(repo)
    cfg = load_mirror_config(repo)
    rows = load_inventory(repo)
    inv = project_inventory(rows, cfg)
    projects = []
    for p in cfg["projects"]:
        acc = p["accession"]
        files = wanted.get(acc)
        i = inv.get(acc) or {"species": Counter(), "countries": Counter(), "years": Counter(),
                             "undated": 0, "bases": 0}
        m = mirror_summary(files, (ledgers or {}).get(acc)) if ledgers is not None else None
        planned_gb = round(sum(int(f["bytes"] or 0) for f in files) / 1e9, 1) if files is not None else None
        years = sorted(i["years"])
        if files is not None and m is not None and files and m["files"] >= len(files):
            state = "complete"
        elif m and m["files"]:
            state = "mirroring"
        else:
            state = "waiting"
        projects.append({
            "accession": acc, "note": p.get("note", ""),
            "species": dict(i["species"]),
            "runs": sum(i["species"].values()),
            "gbases": round(i["bases"] / 1e9),
            "countries": [c for c, _ in i["countries"].most_common()],
            "years": [years[0], years[-1]] if years else None,
            "undated": i["undated"],
            "files_planned": len(files) if files is not None else None,
            "gb_planned": planned_gb,
            "mirror": m, "state": state,
            "url": f"https://www.ebi.ac.uk/ena/browser/view/{acc}",
            "mirror_url": f"{url_base.rstrip('/')}/raw/{acc}/" if url_base else None,
        })

    phases = load_phases(repo)
    known = [x for x in projects if x["gb_planned"] is not None]
    totals = {
        "projects": len(projects),
        "runs": sum(x["runs"] for x in projects),
        "runs_by_species": dict(sum((Counter(x["species"]) for x in projects), Counter())),
        "gb_planned": round(sum(x["gb_planned"] for x in known), 1),
        "planned_complete": len(known) == len(projects),
        "files_planned": sum(x["files_planned"] or 0 for x in projects),
        "mirrored_gb": round(sum((x["mirror"] or {}).get("gb", 0) for x in projects), 1),
        "mirrored_files": sum((x["mirror"] or {}).get("files", 0) for x in projects),
        "projects_complete": sum(x["state"] == "complete" for x in projects),
        "countries": len({c for x in projects for c in x["countries"]}),
        "phases_done": sum(ph.get("status") == "done" for ph in phases),
    }
    return {
        "generated_at": (now or datetime.now(timezone.utc)).isoformat(timespec="seconds"),
        "commit": _git(repo, "rev-parse", "--short", "HEAD"),
        "repo_url": repo_url,
        "ena_available": all(wanted.get(p["accession"]) is not None for p in cfg["projects"]),
        "mirror_available": ledgers is not None,
        "phases": phases,
        "projects": projects,
        "time": sampling_through_time(rows, cfg),
        "totals": totals,
    }


def write(data: dict, out: Path) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, indent=1) + "\n")
    return out
