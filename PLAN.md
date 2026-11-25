# Herring copy number variation from public short-read data: project plan

Drafted 2026-11-25. Dataset numbers come from an ENA query run that day
(`metadata/ena-runs-*.tsv`); depths are estimates from total bases / ~0.8 Gb genome
and have not yet been checked sample by sample.

## 1. Aim and questions

Use every suitable public whole-genome short-read dataset for Atlantic herring
(*Clupea harengus*) and Pacific herring (*C. pallasii*) to build one consistently
processed copy number variation (CNV) resource, then ask:

1. **Landscape.** How much of the herring genome varies in copy number, which genes
   and gene families are involved, and how does CNV diversity compare to SNP diversity?
2. **Space.** Does copy number differ among populations along known environmental
   axes (Baltic salinity gradient, spawning season, latitude/temperature, NE vs NW
   Atlantic, Pacific coast from Washington to Alaska) more than neutral structure predicts?
3. **Time.** Has copy number shifted within populations across decades, particularly
   across documented collapses and warming, beyond what drift predicts?
4. **Parallelism.** Are the same genes or families copy-number-variable, or
   environment-associated, in both species?
5. **Resilience.** Do CNV loci that differ in space also change in time (standing
   variation being used), and are they enriched for functions such as immunity,
   reproduction timing, osmoregulation, or detoxification?

## 2. Data available

### Core WGS datasets (short read, Illumina)

| BioProject | Species | Runs | Total | Approx. design | Space / time |
|---|---|---|---|---|---|
| PRJNA642736 | Atlantic (+100 runs labelled Pacific) | 392 | 4.0 Tb | Mix of pools and individuals | NE Atlantic, Baltic, Arctic, Pacific; 1979–2013 |
| PRJNA930418 | Atlantic | 44 | 1.9 Tb | ~43 Gb/run, likely pools | NW Atlantic (Canada, USA); dates to verify |
| PRJNA356817 | Atlantic | 26 | 1.5 Tb | Deep, likely pools | Metadata sparse |
| PRJNA1137584 | Atlantic | 86 | 0.95 Tb | ~14x individuals | Baltic (Sweden), 2020–2022 |
| PRJEB32358 | Atlantic | 13 | 0.64 Tb | Deep | Baltic |
| PRJEB25669 | Atlantic | 111 | 0.44 Tb | ~5x individuals | Faroe Islands, 2015–2017 |
| PRJNA338612 | Atlantic (+6 Pacific) | 33 | 0.35 Tb | ~13x | North Atlantic, 2014 |
| PRJNA1023520 | Atlantic | 12 Illumina + 14 HiFi | 0.7 Tb | Same fish as the haplotype assemblies | Ireland, Norway, Sweden; 2019–2021 |
| PRJNA179107/109/110 | Atlantic | ~24 | 0.4 Tb | Early pools | Metadata sparse |
| PRJNA972458 | Pacific | 892 | 1.6 Tb | **~2x low-coverage individuals** | Canada, USA; **1991–2017** |
| PRJNA1356517 | Pacific | 120 | 0.54 Tb | ~5–6x individuals | USA; **1996–2023** |
| PRJEB52723, PRJEB77597, PRJEB54658 | Atlantic (ancient) | 344 | 2.7 Tb | Ancient DNA | Baltic and North Sea, Viking Age onward |

Roughly 15 Tb of modern WGS bases in total: 6.7 TB as gzipped FASTQ by ENA's
byte counts (`herring mirror-fastq --plan`, 2026-11-25), an estimated 3–4 TB as CRAM.

**Excluded:** RAD-seq (PRJNA669068, PRJNA669540, PRJNA508972, PRJNA717668), which
cannot give reliable depth-based copy number. **Deferred:** ancient DNA (damage and
fragmentation bias depth; treat as a stretch goal for large, high-copy loci only).
**Secondary:** RNA-seq (PRJNA1112117, PRJNA1277801, PRJNA862443) for testing whether
copy number tracks expression.

### References

| Assembly | Species | Level | Use |
|---|---|---|---|
| Ch_v3.1 (GCF_040183275.1) | Atlantic | Chromosome, RefSeq annotated | Primary reference for both species |
| Clupal_KotzSound (GCA_056790095.1) | Pacific | Chromosome | Secondary reference for Pacific samples |
| 32 haplotype assemblies (GCA_0463342xx–0463348xx; CS, BS, NSSH fish) | Atlantic | Contig | **Truth set** for CNV validation |
| Ch_v2.0.2 (GCF_900700415.2) | Atlantic | Chromosome | Coordinate lookup for published loci only |

### Consequence for design

The data are three different kinds of sample: population pools at 30–50x,
individuals at 5–15x, and individuals at ~2x. The only quantity all three support is
**population-mean copy number per window or gene**. That is the primary unit of
analysis. Individual CNV genotypes are a second tier, limited to samples at ≥10x
(all sizes ≥1 kb) and ≥5x (events ≥10 kb). The ~2x Pacific samples are merged
in silico by population × year into pseudo-pools.

## 3. Methods

### 3.1 Metadata curation (do this first)

- Pull BioSample records for every run; resolve pool vs individual, pool size,
  population name, lat/long, collection date, spawning season, library prep
  (PCR-free or not), instrument, read length.
- Fill gaps from the source papers' supplementary tables. Several projects have
  empty country/date fields, and PRJNA930418 reports a 1905 collection year that
  needs checking.
- Resolve why 100 runs in PRJNA642736 carry the Pacific herring taxon ID.
- Output: one `metadata/samples.tsv` that drives the workflow, with an
  include/exclude column and reason.

### 3.2 Mirroring reads to gannet (ahead of any compute)

Same approach as project-white-whale (`docs/fastq_mirror.md` there), which
measured ENA at 0.04–1 MB/s per connection from raven against ~25 MB/s from
NCBI's AWS open-data bucket.

- Runs on **raven**, nightly from cron, capped per night (GB and concurrent
  runs), until every included run is on gannet. Needs no Hyak time.
- Per run: locate the `.sra` in `sra-pub-run-odp` with NCBI's locator, download
  and check its MD5, `fasterq-dump --split-3`, compress with pigz under ENA's
  file names, check read and base counts against ENA. Fall back to ENA's FASTQ
  when NCBI has no full-quality copy or a check fails.
- On gannet: `raw/<BIOPROJECT>/<RUN>_{1,2}.fastq.gz` plus a per-project
  `mirror_manifest.tsv` ledger (source, MD5 and bytes of the file on gannet,
  ENA MD5, read count, date). A run already in the ledger is skipped, so the
  job is safe to rerun.
- Order: pilot samples, then the truth-set fish (PRJNA1023520), then the
  Pacific temporal series, then the rest. Ancient DNA is not mirrored unless it
  comes into scope.
- Size: 6.7 TB for the modern WGS (3,500 files in 13 BioProjects). At ~25 MB/s
  that is about 3 days of transfer in total; conversion and compression on raven
  are the slower part.
- Tool: `herring mirror-fastq`, vendored from white-whale
  (`docs/fastq_mirror.md`, `VENDORED_FROM_WHITE_WHALE.md`). Projects and order
  are in `data/mirror_projects.yaml`.

### 3.3 Alignment

- Reads are staged from the gannet mirror, never straight from ENA or NCBI:
  `fastp` → `bwa-mem2` → `samtools markdup` → CRAM.
- All samples to Ch_v3.1 for cross-species comparability. Pacific samples also to
  the Kotzebue Sound assembly for within-Pacific analyses; a whole-genome
  alignment (minimap2/Cactus) links coordinates.
- Per-sample QC: mapping rate, insert size, duplication, GC-bias curve, depth
  evenness, sex (depth over the male-specific region), contamination/species check.

### 3.4 Copy number estimation

- **Depth track:** `mosdepth` in 1 kb windows and per gene, computed twice (all
  reads, and MAPQ ≥ 20). GC correction per sample; normalisation to the median of
  single-copy BUSCO genes; mappability mask from `genmap`.
- **Population tier (all data):** mean normalised copy number per population ×
  time point, with between-sample or bootstrap variance.
- **Individual tier (≥10x):** `CNVpytor` for depth, `delly` and `manta` (or
  `smoove`) for breakpoint evidence, merged with `SURVIVOR`; genotyped across
  samples. Only depth calls for 5–10x samples, ≥10 kb.
- **Gene family tier:** for multicopy families collapsed in the reference
  (immune genes especially), family-level copy number from total depth across all
  paralogs, so the answer does not depend on multi-mapping assignment.

### 3.5 Validation

- **Truth set:** call CNVs directly from the 16 fish with haplotype-resolved
  assemblies, then compare with short-read calls from the same fish's Illumina
  data (PRJNA1023520). Downsample those to 15x, 5x, 2x to get sensitivity and
  false discovery by event size at each depth tier. This sets the size thresholds
  used everywhere else.
- **Positive controls:** the male-specific region should separate sexes by depth;
  previously published CNV loci (to compile in the literature review) should replicate.
- **Negative control:** BUSCO single-copy genes should sit at copy number 2.

### 3.6 Controlling batch effects

Study is confounded with geography, which is the main threat to every spatial
conclusion. Controls:

- Within-study contrasts are primary; cross-study comparisons are secondary.
- Populations sampled by more than one study serve as anchors to estimate batch offsets.
- Library type and instrument enter models as covariates.
- A spatial CNV candidate needs either breakpoint support or replication in a
  second study.
- Temporal comparisons are only made within a single project (PRJNA972458 and
  PRJNA1356517 for Pacific; PRJNA642736 for Atlantic).

### 3.7 Population analyses

- **Neutral baseline:** SNP allele frequencies from the same alignments (ANGSD
  genotype likelihoods for individuals, `poolfstat` for pools) to give FST, a
  population covariance matrix, and Ne estimates.
- **Space:** V_ST per locus against the SNP FST distribution; PCA on copy number;
  genotype–environment association (BayPass or RDA with the SNP covariance
  matrix) against salinity, temperature, latitude, and spawning season.
- **Time:** change in mean copy number between time points within a population,
  tested against drift simulated from the SNP-based Ne; consistency of direction
  across replicate populations.
- **Parallelism:** overlap of CNV genes and of environment-associated CNV genes
  between species, against permutation expectation.
- **Function:** GO/family enrichment; overlap with the known adaptive loci and
  the four large inversions; copy number vs expression in the RNA-seq sets.

## 4. Compute and storage

Follows project-white-whale (`docs/compute_profiles.md` and `data/data_host.yaml`
there): gannet holds everything, raven feeds it, Hyak does the bulk compute on
the checkpoint partition.

### gannet: the single home for data

`/volume1/v1_web/project-herring-copy`, served at
`https://gannet.fish.washington.edu/v1_web/project-herring-copy/`.

| Path | Contents |
|---|---|
| `raw/<BIOPROJECT>/` | Mirrored FASTQ and `mirror_manifest.tsv` |
| `reference/<ASSEMBLY>/` | Genomes, annotation, bwa-mem2 indexes, mappability and BUSCO masks |
| `alignments/<BIOPROJECT>/<sample>/` | CRAM, index, per-sample QC |
| `depth/<BIOPROJECT>/<sample>/` | mosdepth window and gene tracks |
| `calls/` | CNV/SV calls, SNP frequency files, copy number matrices |

- Nothing large is committed. The repo holds a manifest per output set (URL,
  SHA-256, size), as white-whale does for methylomes.
- Expected footprint: 6.7 TB FASTQ plus 3–4 TB CRAM plus small derived files,
  about 11 TB. Gannet had 171 TB free when white-whale checked on 2026-10-23.
- Hyak copies are working copies only and are cleared once results are on
  gannet and verified against the manifest.
- Raven has key access to gannet. Klone needs its own key to upload; reading
  from gannet works over HTTPS without one.

### Hyak klone

One Nextflow workflow (as in white-whale) with `desktop`, `raven` and `klone`
profiles; no config names a machine-specific path.

**Default: `ckpt-all`, account `srlab`, every task capped at 4 h 30 min and
retried six times.** Work is shaped to fit that window rather than moved off it:

- FASTQ is split into fixed-size chunks before alignment, so the ~43 Gb pool
  runs align as many short tasks; chunks are merged before duplicate marking.
- Depth, CNVpytor, SV calling and ANGSD run per sample or per chromosome.
- Every task is idempotent, so preemption costs one task, not a sample.

| Work | Where |
|---|---|
| Trimming, chunked alignment, merge, markdup, mosdepth, per-sample and per-chromosome calling | `ckpt-all` |
| Interactive tests | `ckpt-all` via `salloc` |
| Nextflow head job (2 CPUs, 8 GB, days) | Non-preemptible, first with room |
| Staging FASTQ from gannet; uploading results to gannet | Non-preemptible, first with room |
| Any task that still cannot finish in 4.5 h after splitting | Non-preemptible, as a per-process override in a run-only config |

Non-preemptible order, checked with `hyakalloc`: `coenv compute` (40 CPUs,
188 GB), then `srlab cpu-g2-mem2x` (32 CPUs, 503 GB), then `coenv cpu-g2`
(32 CPUs, 200 GB).

Rules carried over from white-whale:

- Work directory, `NXF_TEMP` and logs on `/gscratch/srlab/$USER/herring`, not on
  scrubbed, whose shared quota fills and killed two white-whale runs on
  2026-11-25. Staged FASTQ, which is only read, can sit on
  `/gscratch/scrubbed/$USER/herring/raw`.
- Read-sized intermediates on the node's `/scr`. Do not export `TMPDIR` in the
  head job.
- Stage by BioProject in batches sized to the space available, run, upload,
  verify, clear, then stage the next batch.
- Hyak needs interactive login with 2FA, so head jobs are submitted by a person.

### raven and desktops

| Work | Where |
|---|---|
| Nightly FASTQ mirror to gannet | raven (cron) |
| Pilot runs and workflow debugging | raven (`raven` profile, Singularity) |
| Assembly-based truth set, whole-genome alignment, joint genotyping, BayPass | raven |
| Metadata curation, statistics on copy number matrices, figures, writing | desktop |

On raven, scratch goes on the NVMe RAID0 (`/home/shared/16TB_HDD_01`) and
anything kept goes on a bulk disk until it is on gannet.

Unmeasured so far: gannet → klone transfer rate, and alignment cost per Gb
(a guess of 10^4 core-hours overall). The pilot measures both.

## 5. Phases

| Phase | Work | Deliverable |
|---|---|---|
| 0. Inventory | Metadata curation; literature review of herring CNV, structural variation, and each source dataset | `metadata/samples.tsv`; reference list; list of published CNV loci |
| 0b. Mirror | Create the gannet directory; install the mirror tool (in the repo, tested offline) on raven; nightly mirror from raven, starting as soon as the run list is fixed and continuing in the background | All included FASTQ on gannet with ledgers |
| 1. Pilot | ~40 samples covering every data type and study, aligned on raven from the mirror; build and test the workflow; one small study end to end on `ckpt-all` | Working workflow with klone profile; measured cost per sample; QC report |
| 2. Validation | Truth set from haplotype assemblies; downsampling experiment | Sensitivity/FDR by depth and size; fixed thresholds |
| 3. Scale-out | Full alignment and depth on `ckpt-all`, one BioProject batch at a time | CRAMs, QC and depth tracks on gannet with manifests |
| 4. CNV catalogue | Population and individual tiers; gene family copy numbers | Copy number matrices; annotated CNV catalogue (answers Q1) |
| 5. Space | V_ST, GEA, within- then across-study | Spatial candidates with batch checks (Q2) |
| 6. Time | Pacific temporal series first, then Atlantic | Temporal candidates vs drift (Q3) |
| 7. Synthesis | Parallelism, function, expression | Q4, Q5; manuscript figures |
| Stretch | Ancient DNA at large high-copy loci; pangenome graph genotyping | Deep-time trajectories |

Phases 1 and 2 gate everything else: if the downsampling experiment shows that
2x pseudo-pools cannot recover known CNVs, the Pacific temporal analysis is
restricted to the ~5x dataset and to large events.

## 6. Risks

- **Batch confounded with geography.** Addressed in 3.6; some cross-study
  spatial comparisons may remain unresolvable.
- **Low coverage in the largest Pacific dataset.** Population-level only; power
  set by phase 2.
- **Reference bias for Pacific reads on an Atlantic reference.** Divergent or
  absent sequence looks like a deletion. Mitigated by dual mapping and by
  requiring within-species polymorphism.
- **Pool-seq depth reflects unequal DNA contribution per fish.** Affects variance
  more than mean; pool size enters as an effective sample size.
- **Incomplete public metadata.** Samples without recoverable location or date
  are used for the catalogue but not for spatial or temporal tests.
- **Storage on Hyak.** Gannet has room for everything; the constraint is the
  lab fileset and scrubbed on klone, handled by batching per BioProject.
- **Preemption.** Chunked, idempotent tasks keep the loss per preemption small;
  tasks that repeatedly hit the 4.5 h limit are split further before being moved.

## 7. Decisions needed

1. Whether ancient DNA is in scope from the start or stays a stretch goal
   (recommended: stretch; it would add about 2.7 Tb of reads to the mirror).
2. Whether to contact the data generators early about metadata gaps and overlap
   with their own plans.
