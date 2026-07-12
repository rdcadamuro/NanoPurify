# NanoPurify

[![DOI](https://zenodo.org/badge/1298314004.svg)](https://zenodo.org/badge/latestdoi/1298314004)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**A consensus-binning and chimerism-aware decontamination pipeline for recovering
metagenome-assembled genomes (MAGs) from Oxford Nanopore metagenomic sequencing,
ending in an ENA-submission-ready metadata table.**

Built and organized by **PhD Rafael Dorighello Cadamuro**.

NanoPurify orchestrates an ensemble of eight independent binners behind DAS Tool,
then runs a GUNC-driven decision loop (Deepurify → MAGpurify2) to clean up
chimeric bins before a final CheckM2 quality gate. Every decision made about
every bin — kept as-is, cleaned, or discarded, and why — is logged to a
per-sample audit table.

This repository contains only the **orchestration code** (a Snakefile and three
Python scripts). It calls, via separate conda environments, a set of third-party
tools — no third-party source code is vendored here. See
[Software and citations](#software-and-citations) for full credit.

## Table of contents

- [Pipeline overview](#pipeline-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running the pipeline](#running-the-pipeline)
- [Outputs](#outputs)
- [Design notes](#design-notes)
- [Software and citations](#software-and-citations)
- [Citing NanoPurify](#citing-nanopurify)
- [License](#license)

## Pipeline overview

```
reads.fastq(.gz)
      │
      ▼
   Flye (--nano-hq/--nano-raw, --meta)  ──►  assembly.fasta + assembly_graph.gfa
      │
      ▼
   Medaka (polishing)  ──►  consensus.fasta
      │
      ▼
   minimap2 + samtools + jgi_summarize_bam_contig_depths
      │   (one shared BAM/depth file, reused by every binner below)
      ▼
   ┌──────────────────────────────────────────────────────────────┐
   │  MetaBAT2 · MaxBin2 · CONCOCT · VAMB · TaxVAMB (+ mmseqs2)    │
   │  SemiBin2 · MetaCoAG (uses the Flye assembly graph) · LRBinner│
   └──────────────────────────────────────────────────────────────┘
      │   (each binner's output standardized to a contig2bin.tsv)
      ▼
   DAS Tool  ──►  consolidates whichever binners actually produced bins
      │            (fault-tolerant: binners that fail or return nothing
      │             are skipped rather than halting the run)
      ▼
   CheckM2 (pass 1)  ──►  initial filter: completeness ≥50%, contamination ≤10%
      │
      ▼
   barrnap  ──►  16S/23S/5S rRNA presence (ENA MIMAG-checklist field)
      │
      ▼
   For every bin that passed CheckM2 (pass 1):
      GUNC (GTDB r95 diamond database)  ──►  clade separation score (CSS)
         │
         ├─ CSS ≤ 0.45 (clean)  ─────────────────────────────────►  keep
         │
         └─ CSS > 0.45 (chimeric)
                │
                ▼
             Deepurify  ──►  GUNC re-check
                │
                ├─ passed  ─────────────────────────────────────►  keep
                │
                └─ still chimeric
                       ▼
                    MAGpurify2 (applied on top of the Deepurify output)
                       │
                       ▼
                    GUNC re-check
                       │
                       ├─ passed  ────────────────────────────►  keep
                       └─ still chimeric  ────────────────────►  DISCARD
      │
      ▼
   CheckM2 (pass 2, final)  ──►  re-checks 50%/10% after cleaning
      │                          (contig removal can shift the numbers)
      ▼
   Approved MAG  ──►  row in the ENA `genome_upload` input TSV
```

Every bin's full trajectory — CSS at each stage, which cleaner (if any) was
applied, completeness/contamination before and after, final verdict — is
written to `audit_table.tsv`, one row per bin.

## Requirements

- Linux, conda/mamba
- Snakemake ≥9 (only needed in its own environment, see below)
- ~20 separate conda environments, one per tool (keeps dependency conflicts
  isolated — several of these tools pin incompatible versions of numpy/PyTorch)
- A GUNC diamond database (`gunc download_db <path> -db gtdb`)
- A CheckM2 diamond database (`checkm2 database --download`)
- A Deepurify model database (see [Deepurify's repository](https://github.com/zoubohao/Deepurify) for the download link)
- An mmseqs2 taxonomy database for TaxVAMB (`mmseqs databases GTDB <path> <tmp>` — use the **full GTDB**, not a small test database, for production runs)

## Installation

Each tool runs in its own conda environment (created with `mamba create -p <path> -c bioconda -c conda-forge <package>`); `config.yaml` just needs the path to each. Exact package/version pins used during development:

| Environment | Package(s) |
|---|---|
| `flye` | flye=2.9.6, medaka=2.1.1, minimap2=2.30 |
| `binning` | metabat2=2.18, maxbin2=2.2.7, concoct=1.1.0 |
| `vamb` | vamb=5.0.4 (provides both `vamb bin default` and `vamb bin taxvamb`) |
| `semibin` | semibin=2.3.0 |
| `metacoag` | metacoag=1.2.2 |
| `lrbinner` | cloned from [anuradhawick/LRBinner](https://github.com/anuradhawick/LRBinner) + its documented conda dependencies (no versioned release) |
| `dastool` | das_tool=1.1.7 |
| `checkm2` | checkm2=1.1.0 |
| `gunc` | gunc=1.0.6 |
| `deepurify` | Deepurify=2.4.3 (pip; installed alongside prodigal, hmmer, galah, concoct, metabat2, semibin per its own requirements) |
| `magpurify` | magpurify=2.1.2 (pip package name for [MAGpurify2](https://github.com/apcamargo/magpurify2)) |
| `mmseqs2` | mmseqs2=18.8cc5c |
| `barrnap` | barrnap=1.10.6 |
| `genome_uploader` | genome-uploader=3.0.2 |
| (runner) | snakemake-minimal≥9 |

## Configuration

Edit `config.yaml`:

- `samples`: map of sample id → path to the nanopore FASTQ(.gz)
- `outdir`: where results are written
- `envs`: absolute path to each conda environment above
- `db`: paths to the GUNC, CheckM2, Deepurify, and mmseqs2/GTDB databases
- `thresholds`: `min_completeness` (default 50), `max_contamination` (default 10), `gunc_css_max` (default 0.45)
- `dastool_score_threshold`: DAS Tool's bin-score cutoff (default 0.5 — the DAS Tool default; only lower this to debug a tiny/toy dataset)

## Running the pipeline

```bash
cd NanoPurify

# 1. assembly -> polishing -> shared coverage -> 8 binners -> DAS Tool -> CheckM2 (pass 1) -> barrnap
snakemake -j <threads> -p --keep-going

# 2. decontamination loop: GUNC -> Deepurify -> GUNC -> MAGpurify2 -> GUNC -> CheckM2 (pass 2)
python scripts/gunc_decontam_loop.py \
  --sample-dir <outdir>/<sample> \
  --gunc-env <path> --gunc-db <gunc_db_gtdb95.dmnd> \
  --deepurify-env <path> --deepurify-db <Deepurify-DB> \
  --magpurify-env <path> \
  --checkm2-env <path> --checkm2-db <checkm2_db.dmnd>

# 3. build the ENA `genome_upload` input table for the approved MAGs
python scripts/build_ena_tsv.py \
  --sample-dir <outdir>/<sample> \
  --sample-accession <SAMEA.../SRS...>
```

Step 2 is a standalone Python script rather than a Snakemake rule because its
control flow (try cleaner 1 → recheck → try cleaner 2 → recheck → discard) is
inherently sequential/conditional and doesn't map cleanly onto Snakemake's
static DAG model — see [Design notes](#design-notes).

## Outputs

Per sample, under `<outdir>/<sample>/`:

- `dastool/dastool_DASTool_bins/*.fa` — consensus bins from DAS Tool
- `checkm2_1/quality_report.tsv` — initial completeness/contamination
- `rrna/rrna_presence.tsv` — 16S/23S/5S presence per bin
- `gunc_loop/<bin_id>/` — GUNC/Deepurify/MAGpurify2 intermediate results per bin
- `checkm2_2/quality_report.tsv` — final completeness/contamination (post-cleaning)
- `mags_approved/*.fa` — final approved MAGs (medium-quality-or-better)
- `audit_table.tsv` — one row per bin: every CSS/CheckM2 value at every stage and the final verdict
- `ena_genome_uploader_input.tsv` — ready for [`genome_upload`](https://github.com/EBI-Metagenomics/genome_uploader) (the `NCBI_lineage` column is left as `PENDING_GTDBTK_ONLINE`; fill it in after running GTDB-Tk separately)

## Design notes

- **Eight binners, not the usual three or four.** Alongside the standard
  MetaBAT2/MaxBin2/CONCOCT trio, NanoPurify includes TaxVAMB (taxonomy-informed
  binning, reported to recover more near-complete bins specifically on
  real long-read data), MetaCoAG (the only binner here that uses the Flye
  assembly graph directly rather than flat contigs), and LRBinner (built
  specifically for long-read binning). DAS Tool is designed to consume an
  arbitrary number of binning sets, so this is a legitimate use of the tool,
  not an unusual one.
- **Fault-tolerant DAS Tool step.** Individual binners routinely fail or
  return zero bins on low-coverage or highly fragmented data (this happened
  repeatedly during testing on a deliberately tiny dataset). Rather than
  letting one failed binner halt the whole pipeline, the binner list handed
  to DAS Tool is built dynamically from whichever binners actually produced
  output.
- **GUNC, not CheckM2, drives chimerism detection.** CheckM2 measures
  completeness/contamination via marker-gene copy number; it is not designed
  to detect cross-species mixing within a single bin. GUNC's clade separation
  score is the right signal for that, which is why the decontamination loop
  is keyed on GUNC's pass/fail rather than CheckM2's.
- **Two CheckM2 passes.** Running GUNC/Deepurify/MAGpurify2 is expensive;
  pass 1 filters out bins that were never going to meet the quality bar
  before spending that compute. Pass 2 exists because contig removal during
  cleaning can push a bin's completeness below the threshold it passed
  before cleaning.
- **Snakemake for the DAG, Python for the decision loop.** The
  assembly→binning→consensus stage is a clean, mostly-parallel DAG (a good
  fit for Snakemake). The decontamination loop's "try cleaner 1, recheck,
  try cleaner 2, recheck, discard" logic is inherently sequential and
  data-dependent; forcing it into Snakemake's declarative rule model would
  have meant Snakemake checkpoints with awkward branching, so it was kept as
  a plain, readable Python script instead.

## Software and citations

NanoPurify does not reimplement any of these tools — it only calls them.
Please cite the original works, not this repository, when reporting
scientific results produced with these tools.

| Stage | Tool | Version used | Citation | DOI |
|---|---|---|---|---|
| Assembly | [Flye](https://github.com/mikolmogorov/Flye) | 2.9.6 | Kolmogorov et al. 2019, *Nat. Biotechnol.* 37:540–546 | [10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8) |
| Assembly (metagenome mode) | metaFlye | 2.9.6 | Kolmogorov et al. 2020, *Nat. Methods* 17:1103–1110 | [10.1038/s41592-020-00971-x](https://doi.org/10.1038/s41592-020-00971-x) |
| Polishing | [Medaka](https://github.com/nanoporetech/medaka) | 2.1.1 | Oxford Nanopore Technologies (no peer-reviewed publication; cite the GitHub repository) | — |
| Mapping | [minimap2](https://github.com/lh3/minimap2) | 2.30 | Li 2018, *Bioinformatics* 34:3094–3100 | [10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191) |
| Binning | [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) | 2.18 | Kang et al. 2019, *PeerJ* 7:e7359 | [10.7717/peerj.7359](https://doi.org/10.7717/peerj.7359) |
| Binning | [MaxBin2](https://sourceforge.net/projects/maxbin2/) | 2.2.7 | Wu et al. 2016, *Bioinformatics* 32(4):605–607 | [10.1093/bioinformatics/btv638](https://doi.org/10.1093/bioinformatics/btv638) |
| Binning | [CONCOCT](https://github.com/BinPro/CONCOCT) | 1.1.0 | Alneberg et al. 2014, *Nat. Methods* 11:1144–1146 | [10.1038/nmeth.3103](https://doi.org/10.1038/nmeth.3103) |
| Binning | [VAMB](https://github.com/RasmussenLab/vamb) | 5.0.4 | Nissen et al. 2021, *Nat. Biotechnol.* 39:555–560 | [10.1038/s41587-020-00777-4](https://doi.org/10.1038/s41587-020-00777-4) |
| Binning (taxonomy-informed) | TaxVAMB | 5.0.4 | Kutuzova et al. 2024, bioRxiv preprint (as of writing) | [10.1101/2024.10.25.620172](https://doi.org/10.1101/2024.10.25.620172) |
| Binning | [SemiBin2](https://github.com/BigDataBiology/SemiBin) | 2.3.0 | Pan et al. 2023, *Bioinformatics* 39(Suppl_1):i21–i29 | [10.1093/bioinformatics/btad209](https://doi.org/10.1093/bioinformatics/btad209) |
| Binning (graph-aware) | [MetaCoAG](https://github.com/metagentools/MetaCoAG) | 1.2.2 | Mallawaarachchi & Lin 2022, *J. Comput. Biol.* 29(12):1357–1376 | [10.1089/cmb.2022.0262](https://doi.org/10.1089/cmb.2022.0262) |
| Binning (long-read-specific) | [LRBinner](https://github.com/anuradhawick/LRBinner) | — (git HEAD) | Wickramarachchi & Lin 2022, *Algorithms Mol. Biol.* 17:14 | [10.1186/s13015-022-00221-z](https://doi.org/10.1186/s13015-022-00221-z) |
| Bin consensus/refinement | [DAS Tool](https://github.com/cmks/DAS_Tool) | 1.1.7 | Sieber et al. 2018, *Nat. Microbiol.* 3:836–843 | [10.1038/s41564-018-0171-1](https://doi.org/10.1038/s41564-018-0171-1) |
| Chimerism detection | [GUNC](https://github.com/grp-bork/gunc) | 1.0.6 | Orakov et al. 2021, *Genome Biol.* 22:178 | [10.1186/s13059-021-02393-0](https://doi.org/10.1186/s13059-021-02393-0) |
| Quality assessment | [CheckM2](https://github.com/chklovski/CheckM2) | 1.1.0 | Chklovski et al. 2023, *Nat. Methods* 20:1203–1212 | [10.1038/s41592-023-01940-w](https://doi.org/10.1038/s41592-023-01940-w) |
| Decontamination | [Deepurify](https://github.com/zoubohao/Deepurify) | 2.4.3 | Zou et al. 2024, *Nat. Mach. Intell.* 6:1245–1255 | [10.1038/s42256-024-00908-5](https://doi.org/10.1038/s42256-024-00908-5) |
| Decontamination | [MAGpurify2](https://github.com/apcamargo/magpurify2) | 2.1.2 | Camargo et al. 2023, *ISME J.* 17:354–370 (tool first described in this study) | [10.1038/s41396-022-01345-1](https://doi.org/10.1038/s41396-022-01345-1) |
| rRNA detection | [barrnap](https://github.com/tseemann/barrnap) | 1.10.6 | Torsten Seemann (no peer-reviewed publication; cite the GitHub repository) | — |
| Taxonomy (for TaxVAMB) | [mmseqs2](https://github.com/soedinglab/MMseqs2) | 18.8cc5c | Steinegger & Söding 2017, *Nat. Biotechnol.* 35(11):1026–1028 | [10.1038/nbt.3988](https://doi.org/10.1038/nbt.3988) |
| Taxonomy assignment | mmseqs2 taxonomy | 18.8cc5c | Mirdita et al. 2021, *Bioinformatics* 37(18):3029–3031 | [10.1093/bioinformatics/btab184](https://doi.org/10.1093/bioinformatics/btab184) |
| Reference taxonomy | [GTDB](https://gtdb.ecogenomic.org/) | r95 (GUNC db) | Parks et al. 2022, *Nucleic Acids Res.* 50(D1):D785–D794 | [10.1093/nar/gkab776](https://doi.org/10.1093/nar/gkab776) |
| ENA submission | [genome_uploader](https://github.com/EBI-Metagenomics/genome_uploader) | 3.0.2 | EBI-Metagenomics (no dedicated tool publication; cite the GitHub repository) | — |
| Orchestration | [Snakemake](https://github.com/snakemake/snakemake) | ≥9 | Mölder et al. 2021, *F1000Research* 10:33 | [10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2) |

## Citing NanoPurify

If NanoPurify (the orchestration pipeline itself, as opposed to the individual
tools it calls — see the table above) was useful in your work, please cite it
via its archived Zenodo record. Each tagged GitHub release is automatically
archived on Zenodo with its own version-specific DOI; the badge at the top of
this README always resolves to the latest one.

```bibtex
@software{cadamuro_nanopurify,
  author  = {Cadamuro, Rafael Dorighello},
  title   = {{NanoPurify: a consensus-binning and chimerism-aware
              MAG recovery pipeline for Oxford Nanopore metagenomes}},
  url     = {https://github.com/rdcadamuro/NanoPurify},
  doi     = {10.5281/zenodo.PENDING},
  version = {v1.0}
}
```

Replace `10.5281/zenodo.PENDING` with the concept DOI (or the specific
version DOI) shown on the Zenodo record once the v1.0 release has been
archived. `CITATION.cff` in the repository root carries the same information
in a machine-readable format and powers GitHub's "Cite this repository"
button.

## License

The orchestration code in this repository (the `Snakefile` and the scripts in
`scripts/`) is released under the MIT License — see [LICENSE](LICENSE). Each
tool called by the pipeline retains its own license; check the tool's own
repository before redistributing it.
