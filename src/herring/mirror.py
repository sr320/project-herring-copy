"""Mirror public herring FASTQ to the lab data host (gannet) ahead of any compute.

Downloading is the slowest step and needs no compute, so it runs first, a
little at a time, from raven (PLAN.md 3.2). The file list (names, MD5s, read
and base counts) always comes from ENA's read-run report. Each run is then
fetched from one of two sources:

* NCBI (default): the run's .sra from NCBI's AWS open-data bucket, found with
  NCBI's locator service and checked against its MD5, converted here with
  fasterq-dump, gzipped under ENA's file names, and checked against ENA's
  read and base counts. project-white-whale measured ~25 MB/s from raven
  against ENA's ~1 MB/s.
* ENA (fallback, or ``--source ena``): each FASTQ over HTTPS, checked against
  ENA's MD5. Used for a run NCBI has no full-quality copy of, one whose NCBI
  checks fail, or one only partly missing.

Downloads go to a staging directory and resume a partial file. Each finished
file is copied to the data host (data/data_host.yaml) with rsync, its MD5
optionally re-checked there, and recorded. On the data host:

* ``raw/<ACCESSION>/<file>``               every FASTQ, named as ENA names it
* ``raw/<ACCESSION>/mirror_manifest.tsv``  one row per verified file: run, mate,
                                           source URL, MD5 and bytes of the file
                                           as stored, ENA's MD5, the .sra MD5
                                           and read count when converted, date

A converted FASTQ holds the same reads as ENA's but is not byte-identical, so
its MD5 differs from ENA's: check staged copies against the ledger's ``md5``.
ENA's MD5 is kept as ``ena_md5`` so a run ENA reprocesses (new MD5) is fetched
again.

Which runs: the BioProjects, in order, are listed in data/mirror_projects.yaml.
Once metadata/samples.tsv lists runs for a project, exactly the runs it marks
``include`` are mirrored. Until then, every run that matches the config's
library strategies, platforms, library sources and taxa is taken. A file
already in the ledger against ENA's current MD5, and present on the host at
its recorded size, is skipped, so the command can be rerun (from cron) until
everything is mirrored.

Vendored from project-white-whale's src/whale/mirror.py
(VENDORED_FROM_WHITE_WHALE.md). Only the run selection differs.
"""
from __future__ import annotations

import csv
import fcntl
import hashlib
import io
import json
import os
import shlex
import shutil
import subprocess
import threading
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from http.client import HTTPException
from pathlib import Path

import yaml

from herring.checksums import md5_file
from herring.ena import ENAError, fetch_run_report, run_fastq_files
from herring.net import https_url

LEDGER_NAME = "mirror_manifest.tsv"
LEDGER_COLUMNS = ["run_accession", "mate", "file", "url", "md5", "bytes", "source", "ena_md5", "sra_md5",
                  "read_count", "mirrored_on"]
NCBI_LOCATE = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"
TIMEOUT = 120
CHUNK = 1 << 20
RETRIES = 5
BACKOFF = [10, 30, 60, 120, 300]
LEDGER_LOCK = ".ledger.lock"     # a directory beside the ledger: mkdir is atomic everywhere
LEDGER_LOCK_WAIT = 300           # seconds to wait for another writer's lock

_sleep = time.sleep


class MirrorError(RuntimeError):
    """A file could not be mirrored and verified."""


# --- what to mirror ------------------------------------------------------------

MIRROR_CONFIG = Path("data") / "mirror_projects.yaml"
SAMPLE_SHEET = Path("metadata") / "samples.tsv"
_YES = {"yes", "y", "true", "1", "include"}


def load_mirror_config(repo: Path) -> dict:
    """data/mirror_projects.yaml: the BioProjects to mirror, in order, and the run filters."""
    return yaml.safe_load((Path(repo) / MIRROR_CONFIG).read_text())


def excluded_runs(repo: Path) -> dict[str, str]:
    """Runs never to mirror, with the reason (`exclude_runs` in data/mirror_projects.yaml)."""
    return {str(k): str(v) for k, v in (load_mirror_config(repo).get("exclude_runs") or {}).items()}


def included_runs(repo: Path, accession: str) -> set[str] | None:
    """Runs metadata/samples.tsv marks `include` for this BioProject; None if it lists none of its runs.

    The sheet needs the columns run_accession, bioproject and include.
    """
    sheet = Path(repo) / SAMPLE_SHEET
    if not sheet.exists():
        return None
    with sheet.open(newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("bioproject") == accession]
    if not rows:
        return None
    return {r["run_accession"] for r in rows if (r.get("include") or "").strip().lower() in _YES}


def wanted_files(repo: Path, accession: str, library_strategies: list[str] | None = None,
                 fetch_runs: Callable[[str], list[dict]] | None = None) -> tuple[list[dict], list[str]]:
    """The FASTQ files to mirror for one accession, with ENA's current MD5s, and notes on what was left out."""
    runs = (fetch_runs or fetch_run_report)(accession)
    notes = []
    excluded = excluded_runs(repo)
    dropped = sorted(r["run_accession"] for r in runs if r["run_accession"] in excluded)
    if dropped:
        notes.append(f"not mirrored (excluded in {MIRROR_CONFIG}): {dropped}")
        runs = [r for r in runs if r["run_accession"] not in excluded]
    keep_runs = included_runs(repo, accession)
    if keep_runs is not None:
        kept = [r for r in runs if r["run_accession"] in keep_runs]
        gone = keep_runs - {r["run_accession"] for r in kept}
        if gone:
            notes.append(f"{len(gone)} included run(s) no longer in ENA's report: {sorted(gone)[:5]}")
        if len(runs) > len(kept):
            notes.append(f"not mirrored (not marked include in {SAMPLE_SHEET}): {len(runs) - len(kept)} run(s)")
    else:
        cfg = load_mirror_config(repo)
        filters = [("library_strategy", list(library_strategies or cfg["library_strategies"])),
                   ("instrument_platform", list(cfg["platforms"])),
                   ("library_source", list(cfg["library_sources"])),
                   ("tax_id", [str(t) for t in cfg["taxids"]])]
        kept = []
        left_out: dict[str, dict[str, int]] = {}
        for r in runs:
            for field, allowed in filters:
                if r.get(field) not in allowed:
                    why = left_out.setdefault(field, {})
                    why[r.get(field) or "?"] = why.get(r.get(field) or "?", 0) + 1
                    break
            else:
                kept.append(r)
        for field, counts in left_out.items():
            notes.append(f"not mirrored ({field.replace('_', ' ')}): {dict(sorted(counts.items()))}")

    files, no_fastq = [], []
    for r in sorted(kept, key=lambda x: x["run_accession"]):
        source, found, _dropped = run_fastq_files(r)
        if source is None:
            no_fastq.append(r["run_accession"])
            continue
        for idx, (url, md5, size) in enumerate(found):
            files.append({"run_accession": r["run_accession"], "mate": str(idx + 1),
                          "file": url.rsplit("/", 1)[-1], "url": f"ftp://{url}", "md5": md5,
                          "bytes": size, "source": source, "read_count": r.get("read_count") or "",
                          "base_count": r.get("base_count") or "",
                          "library_layout": r.get("library_layout") or ""})
    if no_fastq:
        notes.append(f"{len(no_fastq)} run(s) have no usable FASTQ: {no_fastq[:5]}")
    names = [f["file"] for f in files]
    clashes = sorted({n for n in names if names.count(n) > 1})
    if clashes:
        raise MirrorError(f"{accession}: several runs have files named {clashes[:5]}; "
                          "they cannot share raw/<ACCESSION>/")
    return files, notes


# --- where it goes -------------------------------------------------------------

def read_ledger_text(text: str) -> dict[str, dict]:
    if not text.strip():
        return {}
    return {r["file"]: r for r in csv.DictReader(io.StringIO(text), delimiter="\t")}


def _lock_held(path: str, wait: int) -> MirrorError:
    return MirrorError(f"the ledger lock {path} was held for {wait}s. If no mirror-fastq is "
                       "writing this accession (on any machine), remove that directory and re-run.")


def record_in_ledger(store, acc: str, row: dict, staging: Path, *, wait: int = LEDGER_LOCK_WAIT) -> None:
    """Add one row to the accession's ledger on the store, keeping rows others wrote meanwhile.

    The ledger is rewritten whole, so two writers (two machines, or two staging
    directories) could each drop the other's rows. The store's ledger lock is
    held across re-reading the current ledger, adding the row and writing it.
    """
    store.lock_ledger(acc, wait)
    try:
        rows = store.read_ledger(acc)
        rows[row["file"]] = row
        store.write_ledger(acc, rows, staging)
    finally:
        store.unlock_ledger(acc)


def _ledger_text(rows: dict[str, dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=LEDGER_COLUMNS, delimiter="\t", lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    w.writerows(rows[k] for k in sorted(rows))
    return buf.getvalue()


class LocalStore:
    """A directory on this machine (a mounted share, or a test)."""

    def __init__(self, root: Path):
        self.root = Path(root)
        self.label = str(self.root)

    def _dir(self, acc: str) -> Path:
        return self.root / "raw" / acc

    def listing(self, acc: str) -> dict[str, int]:
        d = self._dir(acc)
        return {p.name: p.stat().st_size for p in d.iterdir() if p.is_file()} if d.is_dir() else {}

    def read_ledger(self, acc: str) -> dict[str, dict]:
        p = self._dir(acc) / LEDGER_NAME
        return read_ledger_text(p.read_text()) if p.exists() else {}

    def lock_ledger(self, acc: str, wait: int) -> None:
        lock = self._dir(acc) / LEDGER_LOCK
        lock.parent.mkdir(parents=True, exist_ok=True)
        deadline = time.monotonic() + wait
        while True:
            try:
                lock.mkdir()
                return
            except FileExistsError:
                if time.monotonic() >= deadline:
                    raise _lock_held(str(lock), wait) from None
                _sleep(1)

    def unlock_ledger(self, acc: str) -> None:
        (self._dir(acc) / LEDGER_LOCK).rmdir()

    def write_ledger(self, acc: str, rows: dict[str, dict], staging: Path) -> None:
        d = self._dir(acc)
        d.mkdir(parents=True, exist_ok=True)
        tmp = d / (LEDGER_NAME + ".tmp")
        tmp.write_text(_ledger_text(rows))
        tmp.replace(d / LEDGER_NAME)

    def put(self, src: Path, acc: str, name: str) -> None:
        d = self._dir(acc)
        d.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(d / name))

    def md5(self, acc: str, name: str) -> str:
        return md5_file(self._dir(acc) / name)


class SSHStore:
    """``user@host:/path`` reached over SSH (key auth) and rsync, as gannet is from raven."""

    def __init__(self, target: str):
        host, sep, path = target.partition(":")
        if not sep or not host or not path:
            raise MirrorError(f"upload target must look like user@host:/path, not {target!r}")
        self.host, self.path, self.label = host, path.rstrip("/"), target

    def _dir(self, acc: str) -> str:
        return f"{self.path}/raw/{acc}"

    def _ssh(self, command: str) -> str:
        proc = subprocess.run(["ssh", "-o", "BatchMode=yes", self.host, command],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise MirrorError(f"ssh {self.host} failed: {proc.stderr.strip() or proc.returncode}")
        return proc.stdout

    def _rsync(self, src: Path, dest: str) -> None:
        proc = subprocess.run(
            ["rsync", "-tp", "--partial", "--chmod=D755,F644", "-e", "ssh -o BatchMode=yes",
             str(src), f"{self.host}:{dest}"], capture_output=True, text=True)
        if proc.returncode != 0:
            raise MirrorError(f"rsync to {self.host}:{dest} failed: {proc.stderr.strip()}")

    def listing(self, acc: str) -> dict[str, int]:
        d = shlex.quote(self._dir(acc))
        out = self._ssh(f"if [ -d {d} ]; then find {d} -maxdepth 1 -type f -printf '%f\\t%s\\n'; fi")
        return {name: int(size) for name, size in (ln.split("\t") for ln in out.splitlines() if ln)}

    def read_ledger(self, acc: str) -> dict[str, dict]:
        p = shlex.quote(f"{self._dir(acc)}/{LEDGER_NAME}")
        return read_ledger_text(self._ssh(f"if [ -f {p} ]; then cat {p}; fi"))

    def lock_ledger(self, acc: str, wait: int) -> None:
        d, lock = shlex.quote(self._dir(acc)), shlex.quote(f"{self._dir(acc)}/{LEDGER_LOCK}")
        # One round trip: wait on the host itself; exit 75 if still held after `wait` s.
        proc = subprocess.run(
            ["ssh", "-o", "BatchMode=yes", self.host,
             f"mkdir -p {d} || exit 1; i=0; until mkdir {lock} 2>/dev/null; do i=$((i+1)); "
             f"[ $i -ge {int(wait)} ] && exit 75; sleep 1; done"],
            capture_output=True, text=True)
        if proc.returncode == 75:
            raise _lock_held(f"{self.host}:{self._dir(acc)}/{LEDGER_LOCK}", wait)
        if proc.returncode != 0:
            raise MirrorError(f"ssh {self.host} failed: {proc.stderr.strip() or proc.returncode}")

    def unlock_ledger(self, acc: str) -> None:
        self._ssh(f"rmdir {shlex.quote(f'{self._dir(acc)}/{LEDGER_LOCK}')}")

    def write_ledger(self, acc: str, rows: dict[str, dict], staging: Path) -> None:
        local = Path(staging) / acc / LEDGER_NAME
        local.parent.mkdir(parents=True, exist_ok=True)
        local.write_text(_ledger_text(rows))
        self._ssh(f"mkdir -p {shlex.quote(self._dir(acc))}")
        self._rsync(local, f"{self._dir(acc)}/{LEDGER_NAME}")

    def put(self, src: Path, acc: str, name: str) -> None:
        self._ssh(f"mkdir -p {shlex.quote(self._dir(acc))}")
        self._rsync(src, f"{self._dir(acc)}/{name}")
        src.unlink()

    def md5(self, acc: str, name: str) -> str:
        return self._ssh(f"md5sum {shlex.quote(f'{self._dir(acc)}/{name}')}").split()[0]


def load_data_host(repo: Path) -> dict:
    return yaml.safe_load((Path(repo) / "data" / "data_host.yaml").read_text())


def default_staging(repo: Path) -> Path:
    if os.environ.get("HERRING_MIRROR_STAGING"):
        return Path(os.environ["HERRING_MIRROR_STAGING"])
    raw = os.environ.get("HERRING_RAW")
    return (Path(raw) if raw else Path(repo) / "data" / "raw") / ".mirror-staging"


# --- downloading ---------------------------------------------------------------

class Throttle:
    """A total download rate shared by every worker thread (bytes/s; None = unlimited)."""

    def __init__(self, rate: float | None):
        self.rate = rate
        self._lock = threading.Lock()
        self._next = time.monotonic()

    def take(self, n: int) -> None:
        if not self.rate:
            return
        with self._lock:
            now = time.monotonic()
            self._next = max(self._next, now) + n / self.rate
            wait = self._next - now
        if wait > 0:
            _sleep(wait)


def _urlopen(req):
    return urllib.request.urlopen(req, timeout=TIMEOUT)  # noqa: S310 - ENA URLs from ENA's report


def download(url: str, dest: Path, *, expected_bytes: int | None = None,
             throttle: Throttle | None = None, opener=None) -> str:
    """Fetch url to dest via dest.part, resuming a partial file; return the MD5 of what was written.

    Network errors are retried with back-off, each time resuming where the
    partial file ends, and so is a response that ends early without an error
    (when the size is known). A server that ignores the Range header restarts
    the file.
    """
    part = dest.with_name(dest.name + ".part")
    if part.exists() and expected_bytes is not None and part.stat().st_size > expected_bytes:
        part.unlink()
    for attempt in range(RETRIES + 1):
        h = hashlib.md5()
        have = 0
        if part.exists():
            with open(part, "rb") as fh:
                for chunk in iter(lambda: fh.read(CHUNK), b""):
                    h.update(chunk)
                    have += len(chunk)
        if expected_bytes is not None and have == expected_bytes:
            break
        req = urllib.request.Request(url, headers={"Range": f"bytes={have}-"} if have else {})
        try:
            with (opener or _urlopen)(req) as resp:
                status = getattr(resp, "status", 200)
                if have and status != 206:
                    have, h = 0, hashlib.md5()
                with open(part, "ab" if have else "wb") as out:
                    for chunk in iter(lambda: resp.read(CHUNK), b""):
                        out.write(chunk)
                        h.update(chunk)
                        have += len(chunk)
                        if throttle:
                            throttle.take(len(chunk))
            if expected_bytes is not None and have < expected_bytes:
                # Retried below like a reset, resuming from here rather than from zero.
                raise ConnectionError(f"response ended at {have:,} of {expected_bytes:,} bytes")
            break
        except (urllib.error.URLError, HTTPException, OSError) as exc:
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (403, 404, 410):
                raise MirrorError(f"{url}: HTTP {exc.code}") from None
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 416:
                part.unlink(missing_ok=True)
            if attempt == RETRIES:
                raise MirrorError(f"{url}: gave up after {RETRIES + 1} attempts ({exc})") from None
            _sleep(BACKOFF[min(attempt, len(BACKOFF) - 1)])
    part.replace(dest)
    return h.hexdigest()


# --- the mirror ----------------------------------------------------------------

class _Lock:
    """One mirror at a time per staging directory, so overlapping cron runs exit quietly."""

    def __init__(self, path: Path):
        self.path = path

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # Append mode, so a run that finds the lock taken leaves the holder's PID in the file.
        self.fh = open(self.path, "a")  # noqa: SIM115 - held for the run
        try:
            fcntl.flock(self.fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.fh.close()
            return False
        self.fh.truncate(0)
        self.fh.write(f"{os.getpid()}\n")
        self.fh.flush()
        return True

    def __exit__(self, *exc):
        if not self.fh.closed:
            fcntl.flock(self.fh, fcntl.LOCK_UN)
            self.fh.close()


def ena_md5(done: dict) -> str:
    """The ENA MD5 a ledger row was mirrored against (rows converted from NCBI hold their own MD5 in md5)."""
    return done.get("ena_md5") or done.get("md5", "")


def plan_mirror(repo: Path, accessions: list[str], store, *, library_strategies=None,
                fetch_runs=None) -> list[dict]:
    """Per accession: the files wanted, which are already mirrored and verified, and what is left."""
    plans = []
    for acc in accessions:
        try:
            files, notes = wanted_files(repo, acc, library_strategies, fetch_runs)
        except (ENAError, MirrorError, urllib.error.URLError, OSError) as exc:
            plans.append({"accession": acc, "error": str(exc), "files": [], "todo": [], "notes": [],
                          "ledger": {}})
            continue
        ledger = store.read_ledger(acc)
        present = store.listing(acc)
        todo = []
        for f in files:
            done = ledger.get(f["file"])
            if done and ena_md5(done) == f["md5"] and f["file"] in present \
                    and str(present[f["file"]]) == str(done.get("bytes") or present[f["file"]]):
                continue
            if done and ena_md5(done) != f["md5"]:
                notes.append(f"{f['file']}: ENA's MD5 changed since it was mirrored; fetching again")
            todo.append(f)
        plans.append({"accession": acc, "files": files, "todo": todo, "notes": notes, "ledger": ledger})
    return plans


# --- checking a gzipped FASTQ's structure -----------------------------------------

def gz_line_count(path: Path, pigz: str | None = None) -> int:
    """Lines in a gzipped file (a last line without a newline counts)."""
    cmd = [pigz, "-dc", str(path)] if pigz else ["gzip", "-dc", str(path)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    lines, last = 0, b"\n"
    for chunk in iter(lambda: proc.stdout.read(CHUNK * 4), b""):
        lines += chunk.count(b"\n")
        last = chunk[-1:]
    proc.stdout.close()
    if proc.wait() != 0:
        raise MirrorError(f"{path.name}: {cmd[0]} -dc exited {proc.returncode}; the file is not valid gzip")
    return lines + (last != b"\n")


def expected_records(f: dict, n_files: int) -> int | None:
    """FASTQ records this file must hold by ENA's read_count, or None when that cannot be said.

    A paired run's two files and a single-end run's one file each hold read_count
    records; any other arrangement (one interleaved file, say) is not checked.
    """
    if not f.get("read_count"):
        return None
    paired = (f.get("library_layout") or "").upper() == "PAIRED"
    if (paired and n_files == 2) or (not paired and n_files == 1):
        return int(f["read_count"])
    return None


def check_fastq_structure(path: Path, f: dict, n_files: int, pigz: str | None = None) -> None:
    """Whole 4-line records, as many as ENA's read_count says. ENA's MD5 alone cannot catch a
    deposited file that is itself malformed (a stray line break inside a read, say)."""
    want = expected_records(f, n_files)
    lines = gz_line_count(path, pigz)
    if lines % 4:
        raise MirrorError(f"{f['file']}: {lines} lines is not whole FASTQ records")
    if want is not None and lines // 4 != want:
        raise MirrorError(f"{f['file']}: {lines // 4} records, ENA reports {want} reads")


# --- NCBI: .sra from NCBI's AWS open-data bucket, converted on this machine -----------

class NcbiUnavailable(MirrorError):
    """The run cannot be taken from NCBI (not on S3, check failed); ENA is used instead."""


def _get_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=TIMEOUT) as fh:  # noqa: S310 - fixed NCBI host
        return json.loads(fh.read().decode())


def locate_sra(run: str) -> dict:
    """The run's full-quality .sra on S3, with NCBI's MD5 and size, from NCBI's locator service."""
    try:
        doc = _get_json(f"{NCBI_LOCATE}?acc={run}&location=s3.us-east-1&filetype=run")
    except (urllib.error.URLError, HTTPException, OSError, ValueError) as exc:
        raise NcbiUnavailable(f"NCBI locator failed ({exc})") from None
    for res in doc.get("result", []):
        for f in res.get("files", []):
            if f.get("type") != "sra" or not f.get("md5"):
                continue
            for loc in f.get("locations", []):
                if loc.get("service") == "s3" and loc.get("link"):
                    return {"url": loc["link"], "md5": f["md5"], "bytes": int(f["size"])}
    raise NcbiUnavailable("no full-quality .sra on NCBI's S3 bucket")


def find_tools() -> dict:
    """fasterq-dump ($HERRING_SRATOOLKIT/bin, else PATH) and a gzip compressor (pigz if present)."""
    kit = os.environ.get("HERRING_SRATOOLKIT")
    fasterq = str(Path(kit) / "bin" / "fasterq-dump") if kit else shutil.which("fasterq-dump")
    if not fasterq or not Path(fasterq).exists():
        raise MirrorError("fasterq-dump not found: set HERRING_SRATOOLKIT to an sratoolkit directory "
                          "(on raven /home/shared/sratoolkit.3.0.2-ubuntu64), or use --source ena")
    return {"fasterq_dump": fasterq, "pigz": os.environ.get("HERRING_PIGZ") or shutil.which("pigz")}


def _count_fastq(src: Path, dest: Path | None, compress_cmd: list[str] | None) -> tuple[int, int]:
    """(reads, bases) of a FASTQ, gzipping it to dest on the way when dest is given."""
    proc = out = None
    if dest is not None:
        out = open(dest, "wb")  # noqa: SIM115 - closed below
        proc = subprocess.Popen(compress_cmd, stdin=subprocess.PIPE, stdout=out)
    lines = bases = 0
    try:
        with open(src, "rb") as fh:
            while True:
                batch = fh.readlines(1 << 22)
                if not batch:
                    break
                for k in range((1 - lines) % 4, len(batch), 4):
                    bases += len(batch[k].rstrip(b"\r\n"))
                lines += len(batch)
                if proc:
                    proc.stdin.write(b"".join(batch))
    finally:
        if proc:
            proc.stdin.close()
            rc = proc.wait()
            out.close()
            if rc != 0:
                raise NcbiUnavailable(f"{' '.join(compress_cmd)} exited {rc}")
    if lines % 4:
        raise NcbiUnavailable(f"{src.name}: {lines} lines is not whole FASTQ records")
    return lines // 4, bases


def convert_sra(sra: Path, outdir: Path, files: list[dict], *, tools: dict, threads: int,
                tmp: Path) -> list[dict]:
    """fasterq-dump --split-3 the .sra, gzip each mate under ENA's file name, and check the counts.

    Reads and bases must equal ENA's read_count / base_count for the run (an
    unpaired leftover file counts towards them and is then discarded, as the
    ENA path drops it). Returns [{file, path, reads}] in mate order.
    """
    run = sra.stem
    tmp.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run([tools["fasterq_dump"], "--split-3", "--threads", str(threads), "--temp", str(tmp),
                           "--outdir", str(outdir), "--force", str(sra)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise NcbiUnavailable(f"fasterq-dump exited {proc.returncode}: {proc.stderr.strip()[-300:]}")
    paired = len(files) == 2
    mates = [outdir / f"{run}_{i + 1}.fastq" for i in range(2)] if paired else [outdir / f"{run}.fastq"]
    if not all(m.exists() for m in mates):
        raise NcbiUnavailable(f"fasterq-dump did not write {[m.name for m in mates]}")
    compress = ([tools["pigz"], "-p", str(max(1, threads // len(mates))), "-c"] if tools.get("pigz")
                else ["gzip", "-c"])
    with ThreadPoolExecutor(max_workers=len(mates)) as pool:
        counts = list(pool.map(lambda mf: _count_fastq(mf[0], outdir / mf[1]["file"], compress),
                               zip(mates, files, strict=True)))
    unpaired = outdir / f"{run}.fastq"
    extra = _count_fastq(unpaired, None, None) if paired and unpaired.exists() else (0, 0)
    for m in [*mates, unpaired]:
        m.unlink(missing_ok=True)

    spots = counts[0][0] + extra[0]
    bases = sum(b for _, b in counts) + extra[1]
    want_reads, want_bases = files[0].get("read_count"), files[0].get("base_count")
    if paired and counts[0][0] != counts[1][0]:
        raise NcbiUnavailable(f"mates differ: {counts[0][0]} vs {counts[1][0]} reads")
    if want_reads and spots != int(want_reads):
        raise NcbiUnavailable(f"{spots} spots, ENA reports {want_reads}")
    if want_bases and bases != int(want_bases):
        raise NcbiUnavailable(f"{bases} bases, ENA reports {want_bases}")
    return [{"file": f["file"], "path": outdir / f["file"], "reads": c[0]}
            for f, c in zip(files, counts, strict=True)]


# --- the mirror ----------------------------------------------------------------

def mirror(repo: Path, accessions: list[str], store, *, staging: Path, jobs: int = 4,
           rate_mb: float | None = None, max_gb: float | None = None, verify_remote: bool = True,
           source: str = "ncbi", threads: int = 4, tmp: Path | None = None, verify_fastq: bool = True,
           library_strategies=None, fetch_runs=None, opener=None,
           log=print, today: str | None = None) -> dict:
    """Download, verify and upload what the plan says is missing; record each file in the ledger.

    Work is done per run. With source "ncbi", a run whose files are all missing
    is taken as one .sra from NCBI's S3 bucket (MD5-checked), converted with
    fasterq-dump and checked against ENA's read and base counts; if any of that
    fails, or only some of a run's files are missing, its files come from ENA.

    Returns {"locked": bool, "plans", "mirrored": [...], "failed": [(what, reason)], "bytes": int,
    "sources": {"ncbi_sra": n, "ena": n}}.
    """
    staging = Path(staging)
    tmp = Path(tmp) if tmp else staging / "tmp"
    today = today or date.today().isoformat()
    result = {"locked": False, "plans": [], "mirrored": [], "failed": [], "bytes": 0,
              "sources": {"ncbi_sra": 0, "ena": 0}}
    tools = find_tools() if source == "ncbi" else None
    with _Lock(staging / "mirror.lock") as got:
        if not got:
            result["locked"] = True
            return result
        plans = plan_mirror(repo, accessions, store, library_strategies=library_strategies,
                            fetch_runs=fetch_runs)
        result["plans"] = plans
        for p in plans:
            if p.get("error"):
                result["failed"].append((p["accession"], p["error"]))

        tasks, budget = [], (max_gb * 1e9 if max_gb else None)
        for p in plans:
            by_run: dict[str, list[dict]] = {}
            for f in p["files"]:
                by_run.setdefault(f["run_accession"], []).append(f)
            todo_runs: dict[str, list[dict]] = {}
            for f in p["todo"]:
                todo_runs.setdefault(f["run_accession"], []).append(f)
            for run, fs in todo_runs.items():
                if budget is not None and budget <= 0:
                    break
                tasks.append((p["accession"], run, fs, len(fs) == len(by_run[run]), len(by_run[run])))
                if budget is not None:
                    budget -= sum(int(f["bytes"] or 0) for f in fs)

        ledger_lock = threading.Lock()
        throttle = Throttle(rate_mb * 1e6 if rate_mb else None)

        def publish(acc: str, row: dict, local: Path) -> int:
            size = local.stat().st_size
            store.put(local, acc, row["file"])
            if verify_remote:
                remote = store.md5(acc, row["file"])
                if remote != row["md5"]:
                    raise MirrorError(f"MD5 on {store.label} is {remote}, not {row['md5']}")
            with ledger_lock:
                record_in_ledger(store, acc, {**row, "bytes": str(size), "mirrored_on": today}, staging)
            return size

        pigz = os.environ.get("HERRING_PIGZ") or shutil.which("pigz")

        def fetch_ena(acc: str, f: dict, n_files: int) -> Path:
            """Download one ENA file to staging, MD5-checked against ENA and structure-checked."""
            if not f["md5"]:
                raise MirrorError(f"{f['file']}: ENA lists no MD5 for it, so it cannot be verified")
            local = staging / acc / f["file"]
            local.parent.mkdir(parents=True, exist_ok=True)
            expected = int(f["bytes"]) if f["bytes"] else None
            got_md5 = download(https_url(f["url"]), local, expected_bytes=expected, throttle=throttle,
                               opener=opener)
            if got_md5 != f["md5"]:
                local.unlink()
                got_md5 = download(https_url(f["url"]), local, expected_bytes=expected,
                                   throttle=throttle, opener=opener)
                if got_md5 != f["md5"]:
                    local.unlink()
                    raise MirrorError(f"{f['file']}: MD5 {got_md5} does not match ENA's {f['md5']} (twice)")
            if verify_fastq:
                try:
                    check_fastq_structure(local, f, n_files, pigz)
                except MirrorError:
                    local.unlink(missing_ok=True)
                    raise
            return local

        def from_ena(acc: str, fs: list[dict], n_files: int) -> list[tuple[str, int]]:
            """All of a run's missing files from ENA; none is published unless every one passes."""
            locals_ = []
            try:
                for f in fs:
                    locals_.append((f, fetch_ena(acc, f, n_files)))
            except MirrorError:
                for _, path in locals_:
                    path.unlink(missing_ok=True)
                raise
            return [(f["file"], publish(acc, {**f, "ena_md5": f["md5"]}, path)) for f, path in locals_]

        def from_ncbi(acc: str, run: str, fs: list[dict]) -> list[tuple[str, int]]:
            info = locate_sra(run)
            work = staging / acc / run
            work.mkdir(parents=True, exist_ok=True)
            sra = work / f"{run}.sra"
            try:
                for attempt in (1, 2):
                    try:
                        got_md5 = download(info["url"], sra, expected_bytes=info["bytes"], throttle=throttle,
                                           opener=opener)
                    except MirrorError as exc:
                        raise NcbiUnavailable(str(exc)) from None
                    if got_md5 == info["md5"]:
                        break
                    sra.unlink()
                    if attempt == 2:
                        raise NcbiUnavailable(f".sra MD5 {got_md5} does not match NCBI's {info['md5']} (twice)")
                outputs = convert_sra(sra, work, fs, tools=tools, threads=threads, tmp=tmp / run)
                sra.unlink()
                done = []
                for f, out in zip(fs, outputs, strict=True):
                    row = {**f, "url": info["url"], "md5": md5_file(out["path"]), "source": "ncbi_sra",
                           "ena_md5": f["md5"], "sra_md5": info["md5"], "read_count": str(out["reads"])}
                    done.append((f["file"], publish(acc, row, out["path"])))
                return done
            finally:
                shutil.rmtree(work, ignore_errors=True)
                shutil.rmtree(tmp / run, ignore_errors=True)

        def one_run(acc: str, run: str, fs: list[dict], whole: bool, n_files: int) -> list[tuple[str, int, str]]:
            ncbi_problem = None
            if source == "ncbi" and whole:
                try:
                    return [(name, size, "ncbi_sra") for name, size in from_ncbi(acc, run, fs)]
                except NcbiUnavailable as exc:
                    ncbi_problem = str(exc)
                    log(f"note {acc}/{run}: NCBI not used ({exc}); fetching from ENA")
            try:
                return [(name, size, "ena") for name, size in from_ena(acc, fs, n_files)]
            except MirrorError as exc:
                if ncbi_problem:
                    raise MirrorError(f"{exc}; NCBI's copy also failed ({ncbi_problem}). If the deposited "
                                      f"data is malformed, add {run} to exclude_runs in {MIRROR_CONFIG}") from None
                raise

        with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
            futures = {pool.submit(one_run, *t): t for t in tasks}
            for fut in as_completed(futures):
                acc, run = futures[fut][:2]
                try:
                    done = fut.result()
                except (MirrorError, OSError, subprocess.SubprocessError) as exc:
                    result["failed"].append((f"{acc}/{run}", str(exc)))
                    log(f"FAIL {acc}/{run}: {exc}")
                    continue
                for name, size, src in done:
                    result["mirrored"].append(f"{acc}/{name}")
                    result["bytes"] += size
                    result["sources"][src] += 1
                    log(f"ok   {acc}/{name} ({size / 1e9:.2f} GB, {'NCBI .sra' if src == 'ncbi_sra' else 'ENA'})")
    return result


def configured_accessions(repo: Path) -> list[str]:
    return [p["accession"] for p in load_mirror_config(repo)["projects"]]
