"""Offline tests for `herring build-dashboard` (src/herring/dashboard.py)."""
from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pytest
from click.testing import CliRunner

from herring.cli import main
from herring.dashboard import build, mirror_summary, sampling_through_time

REPO_ROOT = Path(__file__).resolve().parents[1]
HEADER = ["run_accession", "study_accession", "scientific_name", "library_strategy", "library_source",
          "library_layout", "instrument_platform", "base_count", "country", "collection_date", "study_title"]


def _inv(path: Path, rows: list[list[str]]) -> None:
    path.write_text("\t".join(HEADER) + "\n" + "".join("\t".join(r) + "\n" for r in rows))


@pytest.fixture
def repo(tmp_path):
    r = tmp_path / "repo"
    (r / "data").mkdir(parents=True)
    (r / "metadata").mkdir()
    shutil.copy(REPO_ROOT / "data" / "data_host.yaml", r / "data")
    shutil.copy(REPO_ROOT / "data" / "phases.yaml", r / "data")
    (r / "data" / "mirror_projects.yaml").write_text(
        "library_strategies: [WGS]\nplatforms: [ILLUMINA]\nlibrary_sources: [GENOMIC]\ntaxids: [7950, 30724]\n"
        "projects:\n  - accession: PRJA\n    note: a\n  - accession: PRJB\n    note: b\n")
    wgs = ["WGS", "GENOMIC", "PAIRED", "ILLUMINA"]
    # An older inventory that must be ignored in favour of the newer one.
    _inv(r / "metadata" / "ena-runs-clupea-harengus-2026-01-01.tsv",
         [["OLD1", "PRJA", "Clupea harengus", *wgs, "100", "Norway", "2001-01-01", "t"]])
    _inv(r / "metadata" / "ena-runs-clupea-harengus-2026-11-25.tsv", [
        ["R1", "PRJA", "Clupea harengus", *wgs, "1000000000", "Sweden: Baltic", "2010-05-01", "t"],
        ["R2", "PRJA", "Clupea harengus", *wgs, "1000000000", "Sweden", "2012", "t"],
        ["R3", "PRJA", "Clupea harengus", "RNA-Seq", "TRANSCRIPTOMIC", "PAIRED", "ILLUMINA", "5", "", "", "t"],
        ["R4", "PRJA", "Clupea harengus", "WGS", "GENOMIC", "SINGLE", "PACBIO_SMRT", "5", "", "", "t"],
        ["R5", "PRJA", "Clupea harengus", *wgs, "1", "", "1905", "t"],
        ["R6", "PRJX", "Clupea harengus", *wgs, "1", "", "2010", "t"],
    ])
    _inv(r / "metadata" / "ena-runs-clupea-pallasii-2026-11-25.tsv", [
        ["P1", "PRJB", "Clupea pallasii", *wgs, "2000000000", "USA: Alaska", "2010", "t"],
        ["P2", "PRJB", "Clupea pallasii", *wgs, "2000000000", "", "", "t"],
    ])
    return r


def _files(acc, runs):
    return [{"run_accession": run, "file": f"{run}_{m}.fastq.gz", "md5": f"md5-{run}-{m}", "bytes": "1000000000"}
            for run in runs for m in (1, 2)]


def _row(run, m, md5, source="ncbi_sra", date="2026-11-25"):
    return {"run_accession": run, "file": f"{run}_{m}.fastq.gz", "md5": "local", "ena_md5": md5,
            "bytes": "1000000000", "source": source, "mirrored_on": date}


NOW = datetime(2026, 11, 25, tzinfo=timezone.utc)


def test_build_counts_only_in_scope_runs_from_newest_inventory(repo):
    data = build(repo, wanted={"PRJA": _files("PRJA", ["R1", "R2", "R5"]), "PRJB": _files("PRJB", ["P1", "P2"])},
                 ledgers={"PRJA": {}, "PRJB": {}}, repo_url="u", url_base="https://host/x", now=NOW)
    a, b = data["projects"]
    assert a["runs"] == 3 and a["species"] == {"Atlantic": 3}
    assert a["countries"] == ["Sweden"] and a["years"] == [1905, 2012]
    assert b["undated"] == 1 and b["countries"] == ["USA"]
    assert data["totals"]["runs_by_species"] == {"Atlantic": 3, "Pacific": 2}
    assert data["totals"]["gb_planned"] == 10.0 and data["totals"]["files_planned"] == 10
    assert data["ena_available"] and data["mirror_available"]
    assert a["mirror_url"] == "https://host/x/raw/PRJA/"
    assert [p["id"] for p in data["phases"]][:2] == ["0", "0b"]


def test_mirror_progress_and_state(repo):
    files = _files("PRJA", ["R1", "R2"])
    ledger = {r["file"]: r for r in [
        _row("R1", 1, "md5-R1-1"), _row("R1", 2, "md5-R1-2", source="fastq_ftp", date="2026-11-26"),
        _row("R2", 1, "stale-md5"),          # ENA changed it since: not counted
    ]}
    m = mirror_summary(files, ledger)
    assert m == {"files": 2, "gb": 2.0, "runs": 1, "from_ncbi": 1, "from_ena": 1, "last_mirrored": "2026-11-26"}
    data = build(repo, wanted={"PRJA": files, "PRJB": _files("PRJB", ["P1"])},
                 ledgers={"PRJA": ledger, "PRJB": {f["file"]: _row("P1", f["file"][-10], f["md5"]) for f in _files("PRJB", ["P1"])}},
                 repo_url="u", url_base=None, now=NOW)
    assert [p["state"] for p in data["projects"]] == ["mirroring", "complete"]
    assert data["totals"]["projects_complete"] == 1 and data["totals"]["mirrored_gb"] == 4.0


def test_unreadable_sources_are_flagged_not_guessed(repo):
    data = build(repo, wanted={"PRJA": None, "PRJB": None}, ledgers=None, repo_url="u", url_base=None, now=NOW)
    assert not data["ena_available"] and not data["mirror_available"]
    assert all(p["mirror"] is None and p["gb_planned"] is None for p in data["projects"])
    assert data["totals"]["planned_complete"] is False


def test_sampling_through_time_fills_years_and_sets_aside_early_dates(repo):
    from herring.dashboard import load_inventory
    from herring.mirror import load_mirror_config
    t = sampling_through_time(load_inventory(repo), load_mirror_config(repo))
    assert t["years"] == [2010, 2011, 2012]
    atl, pac = t["series"]
    assert atl["counts"] == [1, 0, 1] and atl["early"] == 1 and atl["undated"] == 0
    assert pac["counts"] == [1, 0, 0] and pac["undated"] == 1


def test_cli_offline(repo, tmp_path, monkeypatch):
    monkeypatch.setenv("HERRING_REPO", str(repo))
    out = tmp_path / "site" / "data.json"
    res = CliRunner().invoke(main, ["build-dashboard", "--no-ena", "--no-mirror", "--out", str(out)])
    assert res.exit_code == 0, res.output
    data = json.loads(out.read_text())
    assert data["totals"]["projects"] == 2 and not data["mirror_available"]
