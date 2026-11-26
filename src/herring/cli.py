"""`herring` command line. The mirror command is vendored from project-white-whale
(VENDORED_FROM_WHITE_WHALE.md)."""
from __future__ import annotations

import sys
from pathlib import Path

import click

from herring import __version__
from herring.paths import find_repo


def _repo() -> Path:
    try:
        return find_repo()
    except FileNotFoundError as exc:
        raise click.ClickException(str(exc)) from None


def _table(rows: list[dict], columns: list[str]) -> str:
    if not rows:
        return "(none)"
    widths = {c: max(len(c), *(len(str(r.get(c, ""))) for r in rows)) for c in columns}
    lines = ["  ".join(c.ljust(widths[c]) for c in columns),
             "  ".join("-" * widths[c] for c in columns)]
    lines += ["  ".join(str(r.get(c, "")).ljust(widths[c]) for c in columns) for r in rows]
    return "\n".join(lines)


@click.group()
@click.version_option(__version__, prog_name="herring")
def main():
    """Herring copy number variation from public short-read data."""


@main.command("mirror-fastq")
@click.option("--accession", "accessions", multiple=True,
              help="BioProject to mirror (repeatable; default: every project in data/mirror_projects.yaml, in that order).")
@click.option("--dest", default=None, type=click.Path(file_okay=False),
              help="Mirror into this local directory instead of the data host's upload_target.")
@click.option("--staging", default=None, type=click.Path(file_okay=False),
              help="Where files are downloaded and checked before upload "
                   "(default: $HERRING_MIRROR_STAGING, else $HERRING_RAW/.mirror-staging).")
@click.option("--source", type=click.Choice(["ncbi", "ena"]), default="ncbi", show_default=True,
              help="ncbi: .sra from NCBI's S3 bucket converted with fasterq-dump, ENA as fallback per run; "
                   "ena: ENA's FASTQ only.")
@click.option("--jobs", default=4, show_default=True, help="Runs fetched at once.")
@click.option("--threads", default=4, show_default=True, help="fasterq-dump / pigz threads per run (ncbi).")
@click.option("--tmp", default=None, type=click.Path(file_okay=False),
              help="fasterq-dump scratch, ideally fast local disk (default: <staging>/tmp).")
@click.option("--limit-rate", "rate_mb", default=None, type=float, help="Total download rate cap, MB/s.")
@click.option("--max-gb", default=None, type=float,
              help="Stop starting new files once this much has been queued in this run.")
@click.option("--library-strategy", "library_strategies", multiple=True,
              help="Mirror runs with this ENA library_strategy (repeatable; default: those in "
                   "data/mirror_projects.yaml). Projects listed in metadata/samples.tsv use its include column.")
@click.option("--no-verify-remote", is_flag=True, help="Skip re-computing each MD5 on the data host.")
@click.option("--no-verify-fastq", is_flag=True,
              help="Skip decompressing ENA downloads to check whole records and ENA's read count.")
@click.option("--plan", "plan_only", is_flag=True, help="Show what is mirrored and what is left; download nothing.")
def mirror_fastq_cmd(accessions, dest, staging, source, jobs, threads, tmp, rate_mb, max_gb,
                     library_strategies, no_verify_remote, no_verify_fastq, plan_only):
    """Prefetch public herring FASTQ to the data host (gannet), verified and resumable.

    By default each run is taken from NCBI as .sra (MD5-checked), converted with
    fasterq-dump and checked against ENA's read and base counts; ENA's FASTQ
    (MD5-checked) is the fallback. Safe to rerun, e.g. nightly from cron on
    raven: files already mirrored are skipped, partial downloads resume, and a
    second copy started while one is running exits at once.
    """
    from herring.mirror import (
        LocalStore,
        MirrorError,
        SSHStore,
        configured_accessions,
        default_staging,
        load_data_host,
        mirror,
        plan_mirror,
    )

    repo = _repo()
    accs = list(accessions) or configured_accessions(repo)
    try:
        store = LocalStore(Path(dest)) if dest else SSHStore(load_data_host(repo)["upload_target"])
        if plan_only:
            plans = plan_mirror(repo, accs, store, library_strategies=list(library_strategies) or None)
        else:
            res = mirror(repo, accs, store, staging=Path(staging) if staging else default_staging(repo),
                         jobs=jobs, rate_mb=rate_mb, max_gb=max_gb, verify_remote=not no_verify_remote, verify_fastq=not no_verify_fastq,
                         source=source, threads=threads, tmp=Path(tmp) if tmp else None,
                         library_strategies=list(library_strategies) or None, log=click.echo)
            if res["locked"]:
                click.echo("Another mirror-fastq is running with this staging directory; nothing to do.")
                return
            plans = res["plans"]
    except MirrorError as exc:
        raise click.ClickException(str(exc)) from None

    just_done = set() if plan_only else set(res["mirrored"])
    rows = []
    for p in plans:
        todo = [f for f in p["todo"] if f"{p['accession']}/{f['file']}" not in just_done]
        if p.get("error"):
            rows.append({"accession": p["accession"], "files": "", "GB": "", "mirrored": "", "left GB": "",
                         "state": f"error: {p['error'][:60]}"})
            continue
        total = sum(int(f["bytes"] or 0) for f in p["files"]) / 1e9
        left = sum(int(f["bytes"] or 0) for f in todo) / 1e9
        done = len(p["files"]) - len(todo)
        rows.append({"accession": p["accession"], "files": len(p["files"]), "GB": f"{total:.1f}",
                     "mirrored": done, "left GB": f"{left:.1f}",
                     "state": "complete" if p["files"] and not todo else ("no files" if not p["files"] else "")})
    if plan_only:
        click.echo(f"Mirror at {store.label}\n")
    else:
        click.echo(f"\nMirrored {len(res['mirrored'])} file(s), {res['bytes'] / 1e9:.1f} GB this run "
                   f"({res['sources']['ncbi_sra']} converted from NCBI, {res['sources']['ena']} from ENA).\n")
    click.echo(_table(rows, ["accession", "files", "GB", "mirrored", "left GB", "state"]))
    for p in plans:
        for n in p.get("notes", []):
            click.echo(click.style(f"  {p['accession']}: {n}", fg="yellow"))
    if not plan_only and res["failed"]:
        click.echo(click.style(f"\n{len(res['failed'])} failed:", fg="red"))
        for what, why in res["failed"]:
            click.echo(click.style(f"  {what}: {why}", fg="red"))
        sys.exit(1)


@main.command("build-dashboard")
@click.option("--out", default=None, type=click.Path(dir_okay=False), help="Default: site/data.json.")
@click.option("--no-ena", is_flag=True, help="Do not ask ENA for each project's FASTQ sizes.")
@click.option("--no-mirror", is_flag=True, help="Do not read the gannet mirror ledgers.")
@click.option("--remote", default="origin", show_default=True)
def build_dashboard_cmd(out, no_ena, no_mirror, remote):
    """Write site/data.json for the GitHub Pages dashboard (site/index.html)."""
    import urllib.error

    from herring.dashboard import build, write
    from herring.ena import ENAError
    from herring.mirror import (
        LEDGER_NAME,
        MirrorError,
        configured_accessions,
        load_data_host,
        read_ledger_text,
        wanted_files,
    )
    from herring.net import get_text

    repo = _repo()
    accs = configured_accessions(repo)
    wanted: dict = {}
    for acc in accs:
        if no_ena:
            wanted[acc] = None
            continue
        try:
            wanted[acc] = wanted_files(repo, acc)[0]
        except (ENAError, MirrorError, urllib.error.URLError, OSError) as exc:
            click.echo(click.style(f"ENA unreadable for {acc}: {exc}", fg="yellow"), err=True)
            wanted[acc] = None

    base = load_data_host(repo).get("url_base")
    ledgers = None
    if base and not no_mirror:
        ledgers = {}
        for acc in accs:
            try:
                ledgers[acc] = read_ledger_text(get_text(f"{base.rstrip('/')}/raw/{acc}/{LEDGER_NAME}"))
            except urllib.error.HTTPError as exc:
                if exc.code != 404:
                    click.echo(click.style(f"mirror ledger unreadable for {acc}: {exc}", fg="yellow"), err=True)
                    ledgers = None
                    break
                ledgers[acc] = {}          # nothing mirrored yet
            except (urllib.error.URLError, OSError) as exc:
                click.echo(click.style(f"gannet unreachable: {exc}", fg="yellow"), err=True)
                ledgers = None
                break

    data = build(repo, wanted=wanted, ledgers=ledgers, repo_url=_repo_url(repo, remote), url_base=base)
    path = write(data, Path(out) if out else repo / "site" / "data.json")
    t = data["totals"]
    click.echo(f"Wrote {path}: {t['projects']} projects, {t['runs']} runs, {t['gb_planned']} GB planned, "
               f"{t['mirrored_gb']} GB on gannet ({t['projects_complete']} complete)")


def _repo_url(repo: Path, remote: str) -> str:
    """The GitHub URL of the repository: from Actions' environment, else the remote."""
    import os
    import re
    import subprocess

    if os.environ.get("GITHUB_REPOSITORY"):
        return f"{os.environ.get('GITHUB_SERVER_URL', 'https://github.com')}/{os.environ['GITHUB_REPOSITORY']}"
    origin = subprocess.run(["git", "remote", "get-url", remote], cwd=repo, capture_output=True,
                            text=True).stdout.strip()
    m = re.search(r"github\.com[:/](.+?)(\.git)?$", origin)
    return f"https://github.com/{m.group(1)}" if m else ""
