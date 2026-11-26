# Mirroring herring FASTQ to gannet

`herring mirror-fastq` runs on raven, a little every night, until every
included run is on gannet. Compute on Hyak then starts from that copy
(PLAN.md 3.2). The code is vendored from project-white-whale
(`VENDORED_FROM_WHITE_WHALE.md`); its `docs/fastq_mirror.md` has the full
account of the checks and the ledger, which are unchanged here.

## What is mirrored

The BioProjects in `data/mirror_projects.yaml`, in that order: the truth-set
fish, the Pacific temporal series, then the Atlantic projects.

- Before `metadata/samples.tsv` lists a project: every run that is WGS,
  Illumina, genomic, and taxid 7950 or 30724. Other runs are reported, not
  mirrored.
- Once `metadata/samples.tsv` has rows for a project (columns `run_accession`,
  `bioproject`, `include`): exactly the runs marked `include` = yes.

On 2026-11-25 the plan came to 3,500 files and 6.7 TB:

| BioProject | Files | GB |
|---|---|---|
| PRJNA1023520 | 24 | 198 |
| PRJNA972458 | 1,784 | 927 |
| PRJNA1356517 | 240 | 250 |
| PRJNA642736 | 784 | 2,373 |
| PRJNA1137584 | 172 | 463 |
| PRJEB25669 | 222 | 260 |
| PRJNA338612 | 66 | 230 |
| PRJNA930418 | 88 | 704 |
| PRJNA356817 | 52 | 764 |
| PRJEB32358 | 26 | 366 |
| PRJNA179107, PRJNA179109, PRJNA179110 | 42 | 201 |

## How each run is fetched

1. The `.sra` is located in NCBI's AWS bucket with NCBI's locator, downloaded
   (resuming a partial file) and checked against NCBI's MD5.
2. `fasterq-dump --split-3`, then pigz, under ENA's file names
   (`<RUN>_1.fastq.gz`, `<RUN>_2.fastq.gz`).
3. Both mates must have the same number of reads, and reads and bases must
   equal ENA's `read_count` and `base_count`.
4. The file is copied to gannet with rsync, its MD5 recomputed there, and a row
   added to `raw/<BIOPROJECT>/mirror_manifest.tsv`.

A run falls back to ENA's FASTQ when NCBI has no full-quality copy or a check
fails. ENA files are checked against ENA's MD5 and then decompressed to check
that each holds whole 4-line records, as many as ENA's read count. A run's files
are published together or not at all. A converted FASTQ is not byte-identical
to ENA's, so check copies against the ledger's `md5`, not ENA's.

### Runs that fail both ways

Some deposited runs are malformed at both archives, so neither copy passes.
The first case was PRJNA179109 (2012): in SRR611629, read 2 has stray bytes at
read ends (one of them a line break, which splits a record), and read 1 has
placeholder qualities. The mirror reports such a run as failed and says so.
Add it to `exclude_runs` in `data/mirror_projects.yaml`, with the reason, so
the nightly mirror stops fetching it; `--plan` lists excluded runs.

## Once, before the first run

On gannet, create the project directory (from raven, which has key access):

```bash
ssh sr320@gannet.fish.washington.edu 'mkdir -p /volume1/v1_web/project-herring-copy/raw && chmod 755 /volume1/v1_web/project-herring-copy /volume1/v1_web/project-herring-copy/raw'
```

On raven, in a checkout of this repo (Python 3.10 or newer):

```bash
python3 -m venv .venv && .venv/bin/pip install -e .
```

```bash
export HERRING_SRATOOLKIT=/home/shared/sratoolkit.3.0.2-ubuntu64
export HERRING_PIGZ=/home/shared/pigz-2.8/bin/pigz
export HERRING_MIRROR_STAGING=/home/shared/8TB_HDD_03/$USER/herring/mirror-staging
.venv/bin/herring mirror-fastq --plan
```

The tool paths are the ones white-whale recorded on raven on 2026-10-23.
White-whale's own mirror may be running nightly on the same disks: use a
separate staging and `--tmp` directory, and check `df -h` before raising
`--jobs`.

## First real run

Mirror the smallest project end to end before starting cron:

```bash
.venv/bin/herring mirror-fastq --accession PRJNA179109 --jobs 2 --threads 4 --tmp /home/shared/16TB_HDD_01/$USER/herring/mirror-tmp
```

The deepest runs are ~43 Gb of bases; `--tmp` needs roughly 3–4 times the
`.sra` size per run in flight.

## Nightly from cron

```cron
# 20:00 every night: up to 500 GB, 3 runs at a time
0 20 * * * cd $HOME/project-herring-copy && HERRING_SRATOOLKIT=/home/shared/sratoolkit.3.0.2-ubuntu64 HERRING_PIGZ=/home/shared/pigz-2.8/bin/pigz HERRING_MIRROR_STAGING=/home/shared/8TB_HDD_03/$USER/herring/mirror-staging .venv/bin/herring mirror-fastq --jobs 3 --threads 4 --tmp /home/shared/16TB_HDD_01/$USER/herring/mirror-tmp --max-gb 500 >> $HOME/herring-mirror.log 2>&1
```

At 500 GB a night the whole set takes about two weeks. A second copy started
while one is running exits at once, and the exit code is 1 when any run failed.

| Option | Meaning |
|---|---|
| `--plan` | Show what is mirrored and what is left; download nothing |
| `--accession PRJ…` | Only these (repeatable) |
| `--source ncbi\|ena` | `ncbi` (default) converts NCBI's `.sra`, ENA as fallback per run; `ena` uses ENA's FASTQ only |
| `--jobs N` | Runs at once (default 4) |
| `--threads N` | fasterq-dump / pigz threads per run (default 4) |
| `--tmp DIR` | fasterq-dump scratch (default `<staging>/tmp`) |
| `--limit-rate MB/s` | Total download cap across all jobs |
| `--max-gb G` | Stop starting new runs once G GB (ENA's sizes) are queued in this run |
| `--no-verify-remote` | Skip the MD5 on gannet |
| `--dest DIR` | Mirror into a local directory instead (testing) |

## Not here yet

Staging from the mirror onto klone (white-whale's `fetch-fastq`) is not
vendored; it comes with the Nextflow workflow in the pilot phase.
