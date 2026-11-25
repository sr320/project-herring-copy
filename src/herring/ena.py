"""ENA's read-run report for a BioProject, and the FASTQ files of each run.

The file list a mirror works from (names, MD5s, read and base counts) always
comes from here. Nothing here downloads sequence data.

Vendored from project-white-whale's src/whale/ena.py (VENDORED_FROM_WHITE_WHALE.md):
only the run report and `run_fastq_files` are kept; the sample-sheet builder
is not. The report also asks for instrument_platform and library_source, and
the expected taxa are the two herring species.
"""
from __future__ import annotations

import re
import urllib.parse

from herring.net import get_text

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

RUN_FIELDS = [
    "run_accession", "sample_accession", "sample_alias", "experiment_accession",
    "experiment_title", "library_name", "library_strategy", "library_layout",
    "library_selection", "library_source", "instrument_platform", "instrument_model",
    "scientific_name", "tax_id", "read_count", "base_count", "fastq_ftp", "fastq_md5",
    "fastq_bytes", "submitted_ftp", "submitted_md5", "submitted_bytes", "submitted_format",
]

# Clupea harengus, Clupea pallasii.
EXPECTED_TAXIDS = ["7950", "30724"]

_FASTQ_NAME = re.compile(r"\.f(ast)?q(\.gz|\.bz2)?$", re.IGNORECASE)
_MATE_SUFFIX = re.compile(r"_R?([12])(_\d+)?\.f(ast)?q(\.gz|\.bz2)?$", re.IGNORECASE)

TIMEOUT = 120


class ENAError(RuntimeError):
    """A query returned nothing usable, or the project is not what was expected."""


def _get(url: str) -> str:
    return get_text(url, timeout=TIMEOUT)


def _run_fields(extra: list[str] | None = None) -> list[str]:
    return list(dict.fromkeys(RUN_FIELDS + (extra or [])))


def run_report_url(bioproject: str, extra_fields: list[str] | None = None) -> str:
    query = urllib.parse.urlencode(
        {"accession": bioproject, "result": "read_run",
         "fields": ",".join(_run_fields(extra_fields)), "format": "tsv"})
    return f"{ENA_PORTAL}?{query}"


def fetch_run_report(bioproject: str, extra_fields: list[str] | None = None) -> list[dict]:
    """One record per sequencing run, from ENA's read_run report."""
    text = _get(run_report_url(bioproject, extra_fields))
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise ENAError(
            f"ENA returned no runs for {bioproject}. Either the accession is wrong, or the "
            "project has no public reads (some BioProjects register metadata only).")
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"), strict=False)) for ln in lines[1:]]


def _split(value) -> list[str]:
    return [v for v in (value or "").split(";") if v]


def run_fastq_files(run: dict) -> tuple[str | None, list[tuple[str, str, str]], list[str]]:
    """The FASTQ files for one run: (source, [(url, md5, bytes)], dropped urls).

    ENA's generated FASTQ (`fastq_ftp`) is preferred; when there is none, the
    submitter's files (`submitted_ftp`) are used if every one is a FASTQ.
    Anything else returns source None. A paired run with a third, unpaired
    file has it dropped and returned. Mates are ordered by their _1/_2 suffix.
    """
    def files(prefix):
        urls, md5s, sizes = (_split(run.get(f"{prefix}_{k}")) for k in ("ftp", "md5", "bytes"))
        return [(u, md5s[i] if i < len(md5s) else "", sizes[i] if i < len(sizes) else "")
                for i, u in enumerate(urls)]

    source, found = "fastq_ftp", files("fastq")
    if not found:
        formats = {f.upper() for f in _split(run.get("submitted_format"))}
        submitted = files("submitted")
        if submitted and formats <= {"FASTQ"} and all(_FASTQ_NAME.search(u) for u, _, _ in submitted):
            source, found = "submitted_ftp", submitted
        else:
            return None, [], []

    dropped: list[str] = []
    if run.get("library_layout") == "PAIRED" and len(found) == 3:
        unpaired = [f for f in found if not _MATE_SUFFIX.search(f[0])]
        if len(unpaired) == 1:
            dropped = [unpaired[0][0]]
            found = [f for f in found if f is not unpaired[0]]
    mates = [_MATE_SUFFIX.search(f[0]) for f in found]
    if len(found) > 1 and all(mates) and len({m.group(1) for m in mates}) == len(found):
        found = [f for _, f in sorted(zip((m.group(1) for m in mates), found, strict=True))]
    return source, found, dropped
