#!/usr/bin/env python3
"""
Loop de decisao GUNC -> Deepurify -> GUNC -> MAGpurify2 -> GUNC -> CheckM2 (final).

So processa bins que ja passaram no CheckM2 #1 (completude >=50%, contaminacao <=10%),
produzidos pelo Snakefile em <outdir>/<sample>/checkm2_1/quality_report.tsv e
<outdir>/<sample>/dastool/dastool_DASTool_bins/*.fa

Para cada bin:
  1. GUNC roda -> CSS (clade_separation_score) no maxCSS_level
     - CSS <= threshold (default 0.45): passa direto pro CheckM2 #2, sem limpeza
     - CSS > threshold: roda Deepurify (limpador 1)
         -> GUNC de novo
            - passa: vai pro CheckM2 #2 (log: limpo_deepurify)
            - nao passa: roda MAGpurify2 (limpador 2, a partir da saida do Deepurify)
                -> GUNC de novo
                   - passa: vai pro CheckM2 #2 (log: limpo_deepurify_magpurify2)
                   - nao passa: DESCARTA (log: descartado_quimerismo_persistente)
  2. CheckM2 #2 roda nos bins que sobraram -> confere de novo >=50%/<=10%
     - passa: MAG aprovado
     - nao passa: DESCARTA (log: descartado_checkm2_final)

Gera <outdir>/<sample>/audit_table.tsv com uma linha por bin.
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from pathlib import Path


def run(cmd: list[str], **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def conda_run(env_path: str, args: list[str], **kw):
    return run(["conda", "run", "--no-capture-output", "-p", env_path, *args], **kw)


def read_checkm2_report(report_tsv: Path) -> dict[str, dict[str, float]]:
    out = {}
    with report_tsv.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["Name"]] = {
                "completeness": float(row["Completeness"]),
                "contamination": float(row["Contamination"]),
            }
    return out


def run_gunc(env_gunc: str, db_path: str, bin_fasta: Path, outdir: Path, threads: int) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    conda_run(env_gunc, [
        "gunc", "run", "-i", str(bin_fasta), "-r", db_path,
        "-o", str(outdir), "-t", str(threads),
    ])
    hits = list(outdir.glob("*.maxCSS_level.tsv")) + list(outdir.glob("GUNC*.tsv"))
    if not hits:
        raise FileNotFoundError(f"GUNC nao produziu maxCSS_level.tsv em {outdir}")
    with hits[0].open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        return {"css": 1.0, "pass_gunc": False, "rrs": None}
    row = rows[0]
    return {
        "css": float(row.get("clade_separation_score", 1.0)),
        "pass_gunc": row.get("pass.GUNC", "False").strip() in ("True", "true", "1"),
        "rrs": float(row["reference_representation_score"]) if row.get("reference_representation_score") else None,
    }


def run_deepurify(env_deepurify: str, db_path: str, bin_fasta: Path, workdir: Path, threads: int) -> Path | None:
    in_dir = workdir / "in"
    out_dir = workdir / "out"
    in_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(bin_fasta, in_dir / bin_fasta.name)
    conda_run(env_deepurify, [
        "deepurify", "clean", "-i", str(in_dir), "-o", str(out_dir),
        "--bin_suffix", bin_fasta.suffix.lstrip("."), "-db", db_path,
        "--gpu_num", "0", "--num_process", str(threads),
    ])
    cleaned = list(out_dir.rglob(f"*{bin_fasta.suffix}"))
    return cleaned[0] if cleaned else None


def run_magpurify2(env_magpurify: str, bin_fasta: Path, workdir: Path, depth_tsv: Path | None) -> Path | None:
    work = workdir / "work"
    work.mkdir(parents=True, exist_ok=True)
    out_fna = workdir / f"cleaned_{bin_fasta.name}"
    conda_run(env_magpurify, ["magpurify", "tetra-freq", str(bin_fasta), str(work)])
    conda_run(env_magpurify, ["magpurify", "gc-content", str(bin_fasta), str(work)])
    if depth_tsv is not None and depth_tsv.exists():
        try:
            conda_run(env_magpurify, ["magpurify", "coverage", str(bin_fasta), str(work),
                                       "--depth-file", str(depth_tsv)])
        except subprocess.CalledProcessError:
            pass
    conda_run(env_magpurify, ["magpurify", "clean-bin", str(bin_fasta), str(work), str(out_fna)])
    return out_fna if out_fna.exists() else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", required=True, help="OUT/<sample>")
    ap.add_argument("--gunc-env", required=True)
    ap.add_argument("--gunc-db", required=True)
    ap.add_argument("--deepurify-env", required=True)
    ap.add_argument("--deepurify-db", required=True)
    ap.add_argument("--magpurify-env", required=True)
    ap.add_argument("--checkm2-env", required=True)
    ap.add_argument("--checkm2-db", required=True)
    ap.add_argument("--css-threshold", type=float, default=0.45)
    ap.add_argument("--min-completeness", type=float, default=50.0)
    ap.add_argument("--max-contamination", type=float, default=10.0)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    sample_dir = Path(args.sample_dir)
    bins_dir = sample_dir / "dastool" / "dastool_DASTool_bins"
    checkm2_1_report = sample_dir / "checkm2_1" / "quality_report.tsv"
    depth_tsv = sample_dir / "mapping" / "depth.txt"

    loop_dir = sample_dir / "gunc_loop"
    loop_dir.mkdir(parents=True, exist_ok=True)

    checkm2_1 = read_checkm2_report(checkm2_1_report)
    audit_rows = []
    survivors: dict[str, Path] = {}

    for bin_fasta in sorted(bins_dir.glob("*.fa")):
        bin_id = bin_fasta.stem
        metrics_1 = checkm2_1.get(bin_id)
        row = {
            "bin_id": bin_id,
            "checkm2_1_completude": metrics_1["completeness"] if metrics_1 else "",
            "checkm2_1_contaminacao": metrics_1["contamination"] if metrics_1 else "",
        }
        if not metrics_1 or metrics_1["completeness"] < args.min_completeness or metrics_1["contamination"] > args.max_contamination:
            row["status_final"] = "descartado_checkm2_inicial"
            audit_rows.append(row)
            continue

        work = loop_dir / bin_id
        current_fasta = bin_fasta
        ferramenta_usada = "nenhuma"

        g1 = run_gunc(args.gunc_env, args.gunc_db, current_fasta, work / "gunc_1", args.threads)
        row["css_inicial"] = g1["css"]

        if g1["css"] > args.css_threshold:
            cleaned = run_deepurify(args.deepurify_env, args.deepurify_db, current_fasta, work / "deepurify", args.threads)
            if cleaned is not None:
                current_fasta = cleaned
                ferramenta_usada = "deepurify"
            g2 = run_gunc(args.gunc_env, args.gunc_db, current_fasta, work / "gunc_2", args.threads)
            row["css_apos_deepurify"] = g2["css"]

            if g2["css"] > args.css_threshold:
                cleaned2 = run_magpurify2(args.magpurify_env, current_fasta, work / "magpurify2", depth_tsv)
                if cleaned2 is not None:
                    current_fasta = cleaned2
                    ferramenta_usada = "deepurify_magpurify2"
                g3 = run_gunc(args.gunc_env, args.gunc_db, current_fasta, work / "gunc_3", args.threads)
                row["css_final"] = g3["css"]
                if g3["css"] > args.css_threshold:
                    row["ferramenta_usada"] = ferramenta_usada
                    row["status_final"] = "descartado_quimerismo_persistente"
                    audit_rows.append(row)
                    continue
            else:
                row["css_final"] = g2["css"]
        else:
            row["css_final"] = g1["css"]

        row["ferramenta_usada"] = ferramenta_usada
        survivors[bin_id] = current_fasta
        audit_rows.append(row)

    # CheckM2 #2 (final) nos sobreviventes
    # input_bins fica FORA do dir de saida do CheckM2: --force limpa o -o inteiro antes de rodar
    checkm2_2_dir = sample_dir / "checkm2_2"
    checkm2_2_input_dir = sample_dir / "checkm2_2_input_bins"
    checkm2_2_input_dir.mkdir(parents=True, exist_ok=True)
    for bin_id, fasta_path in survivors.items():
        dest = checkm2_2_input_dir / f"{bin_id}.fa"
        shutil.copy(fasta_path, dest)

    approved_dir = sample_dir / "mags_approved"
    approved_dir.mkdir(parents=True, exist_ok=True)

    if survivors:
        conda_run(args.checkm2_env, [
            "checkm2", "predict", "-i", str(checkm2_2_input_dir), "-x", "fa",
            "-o", str(checkm2_2_dir), "-t", str(args.threads), "--force",
        ], env={**os.environ, "CHECKM2DB": args.checkm2_db})
        checkm2_2 = read_checkm2_report(checkm2_2_dir / "quality_report.tsv")
    else:
        checkm2_2 = {}

    for row in audit_rows:
        bin_id = row["bin_id"]
        if bin_id not in survivors:
            continue
        m2 = checkm2_2.get(bin_id)
        if not m2:
            row["status_final"] = "descartado_checkm2_final_erro"
            continue
        row["checkm2_2_completude"] = m2["completeness"]
        row["checkm2_2_contaminacao"] = m2["contamination"]
        if m2["completeness"] >= args.min_completeness and m2["contamination"] <= args.max_contamination:
            row["status_final"] = f"mantido_{row['ferramenta_usada']}" if row["ferramenta_usada"] != "nenhuma" else "mantido_limpo"
            shutil.copy(survivors[bin_id], approved_dir / f"{bin_id}.fa")
        else:
            row["status_final"] = "descartado_checkm2_final_pos_limpeza"

    fieldnames = ["bin_id", "checkm2_1_completude", "checkm2_1_contaminacao",
                  "css_inicial", "css_apos_deepurify", "css_final", "ferramenta_usada",
                  "checkm2_2_completude", "checkm2_2_contaminacao", "status_final"]
    audit_path = sample_dir / "audit_table.tsv"
    with audit_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(audit_rows)

    n_approved = sum(1 for r in audit_rows if r["status_final"].startswith("mantido"))
    print(f"MAGs aprovados: {n_approved} / {len(audit_rows)} bins avaliados -> {approved_dir}")
    print(f"Tabela de auditoria: {audit_path}")


if __name__ == "__main__":
    main()
