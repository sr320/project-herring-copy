# project-herring-copy

Copy number variation in Atlantic and Pacific herring from public short-read data.

- [PLAN.md](PLAN.md): aims, datasets, methods, compute and storage.
- [docs/fastq_mirror.md](docs/fastq_mirror.md): mirroring the reads to gannet from raven.
- [Dashboard](https://sr320.github.io/project-herring-copy/): mirror progress, phases and sampling through time, rebuilt by `herring build-dashboard` (`.github/workflows/pages.yml`) on every push and every three hours. Preview locally with `herring build-dashboard` then `python3 -m http.server --directory site`.
- [VENDORED_FROM_WHITE_WHALE.md](VENDORED_FROM_WHITE_WHALE.md): code taken from project-white-whale.

```bash
pip install -e ".[dev]"
```

```bash
pytest
```
