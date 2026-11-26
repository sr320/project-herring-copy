# Code adapted from project-white-whale

This project vendors code from
[project-white-whale](https://github.com/RobertsLab/project-white-whale)
(MIT, Copyright (c) 2026 Roberts Lab, University of Washington) rather than
depending on it. Each entry names the source file and the white-whale commit
that last changed it, so a fix in either project can be carried across by hand.
The checkout copied from was at `05d3f6b`.

| Here | White Whale source | Commit | Changes |
|---|---|---|---|
| `src/herring/mirror.py` (`herring mirror-fastq`) | `src/whale/mirror.py` | `4c374e1` | Run selection only: BioProjects and filters (library strategy, platform, library source, taxa) from `data/mirror_projects.yaml`; `metadata/samples.tsv` `include` column in place of a registered study's `fastq_manifest.tsv`; no backlog or registry. Download, conversion, checks, stores, ledger and locks unchanged. `HERRING_*` variables |
| `src/herring/ena.py` | `src/whale/ena.py` | `4793677` | Only `fetch_run_report` and `run_fastq_files` kept; report also asks for `instrument_platform` and `library_source`; taxa 7950 and 30724 |
| `src/herring/net.py` | `src/whale/net.py` | `4793677` | None |
| `src/herring/checksums.py` | `src/whale/checksums.py` | `4793677` | None |
| `src/herring/cli.py` (`mirror-fastq`, `_table`) | `src/whale/cli.py` | `05d3f6b` | That command only; default accessions from `data/mirror_projects.yaml` |
| `tests/test_mirror.py`, `tests/test_net.py` | same paths | `4062784` | WGS fixtures; selection tests rewritten for the config and sample sheet |
| `src/herring/dashboard.py`, `site/index.html`, `.github/workflows/pages.yml` (`herring build-dashboard`) | `src/whale/dashboard.py`, `site/index.html`, `.github/workflows/pages.yml` | `05d3f6b` | Rewritten for this project: phases, mirror progress per BioProject, sampling through time; page styles, tiles, tooltip and Pages workflow kept |
| `data/data_host.yaml` | same path | `05d3f6b` | This project's directory on gannet; layout for alignments, depth and calls |
| `docs/fastq_mirror.md` | same path | `05d3f6b` | Shortened; herring run selection |

Not vendored: `whale fetch-fastq` (staging from the mirror onto the machine
that runs a workflow), the Nextflow profiles, claims and the registry.

When vendoring more, add a row with the white-whale commit
(`git -C <white-whale checkout> log -1 --format=%h -- <file>`).
