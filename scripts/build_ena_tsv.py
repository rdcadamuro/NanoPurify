#!/usr/bin/env python3
"""
Gera a TSV de entrada para o `genome_upload` (EBI-Metagenomics genome_uploader),
a partir dos resultados do pipeline: audit_table.tsv (status final por bin),
checkm2_2/quality_report.tsv, rrna/rrna_presence.tsv, mapping/depth.txt e
mags_approved/*.fa.

So inclui bins com status_final comecando por 'mantido_' (ou seja: aprovados no
CheckM2 final, medium-quality ou melhor -- o minimo que o ENA aceita).

Colunas exigidas pelo genome_uploader (ordem oficial):
genome_name, genome_path, accessions, assembly_software, binning_software,
binning_parameters, stats_generation_software, completeness, contamination,
genome_coverage, metagenome, co-assembly, broad_environment, local_environment,
environmental_medium, rRNA_presence, NCBI_lineage
"""
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

COLUMNS = [
    "genome_name", "genome_path", "accessions", "assembly_software",
    "binning_software", "binning_parameters", "stats_generation_software",
    "completeness", "contamination", "genome_coverage", "metagenome",
    "co-assembly", "broad_environment", "local_environment",
    "environmental_medium", "rRNA_presence", "NCBI_lineage",
]

BINNING_SOFTWARE = ("MetaBAT2_v2.18;MaxBin2_v2.2.7;CONCOCT_v1.1.0;VAMB_v5.0.4;"
                     "TaxVAMB_v5.0.4;SemiBin2_v2.3.0;MetaCoAG_v1.2.2;"
                     "LRBinner;DAS_Tool_v1.1.7")
BINNING_PARAMETERS = "default parameters; consenso via DAS_Tool (score_threshold=0.5)"
ASSEMBLY_SOFTWARE = "Flye_v2.9.6;Medaka_v2.1.1"
STATS_SOFTWARE = "CheckM2_v1.1.0"


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def bin_coverage(depth_tsv: Path, contig_ids: set[str]) -> float | None:
    if not depth_tsv.exists():
        return None
    depths = []
    with depth_tsv.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("contigName") in contig_ids:
                try:
                    depths.append(float(row["totalAvgDepth"]))
                except (KeyError, ValueError):
                    pass
    return round(statistics.mean(depths), 2) if depths else None


def contigs_in_fasta(fasta_path: Path) -> set[str]:
    ids = set()
    with fasta_path.open() as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split()[0].strip())
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", required=True, help="OUT/<sample>")
    ap.add_argument("--sample-accession", default="", help="Run/sample accession fonte (ex: SAMEAxxxx)")
    ap.add_argument("--metagenome-type", default="biodigester/broiler litter metagenome")
    ap.add_argument("--broad-environment", default="livestock farm")
    ap.add_argument("--local-environment", default="anaerobic biodigester / broiler litter")
    ap.add_argument("--environmental-medium", default="manure slurry / poultry litter")
    ap.add_argument("--co-assembly", default="False")
    args = ap.parse_args()

    sample_dir = Path(args.sample_dir)
    audit = read_tsv(sample_dir / "audit_table.tsv")
    rrna_rows = {r["bin_id"]: r for r in read_tsv(sample_dir / "rrna" / "rrna_presence.tsv")}
    depth_tsv = sample_dir / "mapping" / "depth.txt"
    approved_dir = sample_dir / "mags_approved"

    out_rows = []
    for row in audit:
        if not row["status_final"].startswith("mantido"):
            continue
        bin_id = row["bin_id"]
        fasta_path = approved_dir / f"{bin_id}.fa"
        if not fasta_path.exists():
            continue

        rrna = rrna_rows.get(bin_id, {})
        has_rrna = any(int(rrna.get(k, 0) or 0) > 0 for k in ("has_16S", "has_23S", "has_5S"))

        contig_ids = contigs_in_fasta(fasta_path)
        coverage = bin_coverage(depth_tsv, contig_ids)

        out_rows.append({
            "genome_name": bin_id,
            "genome_path": str(fasta_path),
            "accessions": args.sample_accession,
            "assembly_software": ASSEMBLY_SOFTWARE,
            "binning_software": BINNING_SOFTWARE,
            "binning_parameters": BINNING_PARAMETERS,
            "stats_generation_software": STATS_SOFTWARE,
            "completeness": row["checkm2_2_completude"],
            "contamination": row["checkm2_2_contaminacao"],
            "genome_coverage": coverage if coverage is not None else "",
            "metagenome": args.metagenome_type,
            "co-assembly": args.co_assembly,
            "broad_environment": args.broad_environment,
            "local_environment": args.local_environment,
            "environmental_medium": args.environmental_medium,
            "rRNA_presence": "True" if has_rrna else "False",
            "NCBI_lineage": "PENDING_GTDBTK_ONLINE",
        })

    out_path = sample_dir / "ena_genome_uploader_input.tsv"
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"{len(out_rows)} MAGs elegiveis (medium-quality+) -> {out_path}")
    print("Coluna NCBI_lineage fica 'PENDING_GTDBTK_ONLINE' ate voce rodar o GTDB-Tk "
          "online e atualizar manualmente antes de submeter.")


if __name__ == "__main__":
    main()
