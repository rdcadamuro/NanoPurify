"""Converte a saida de `mmseqs easy-taxonomy --tax-lineage 1` (formato *_lca.tsv)
para o formato de 2 colunas que o VAMB/TaxVAMB espera: 'contigs\\tpredictions',
onde predictions e uma string com os 7 ranks canonicos separados por ';'
(Domain;Phylum;Class;Order;Family;Genus;Species), sem prefixos de rank.
"""
import csv

RANK_PREFIXES = ["d_", "p_", "c_", "o_", "f_", "g_", "s_"]

def lineage_to_canonical(lineage_field: str) -> str:
    if not lineage_field:
        return ""
    names = []
    for token in lineage_field.split(";"):
        token = token.strip()
        for pref in RANK_PREFIXES:
            if token.startswith(pref):
                name = token[len(pref):].strip()
                if name:
                    names.append(name)
                break
    return ";".join(names)

def main(input_tsv: str, output_tsv: str):
    with open(input_tsv, newline="") as fin, open(output_tsv, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["contigs", "predictions"])
        reader = csv.reader(fin, delimiter="\t")
        for row in reader:
            if len(row) < 5:
                continue
            contig = row[0]
            lineage_field = row[4] if len(row) > 4 else ""
            canonical = lineage_to_canonical(lineage_field)
            if canonical:
                writer.writerow([contig, canonical])

if __name__ == "__main__":
    main(snakemake.input.tax, snakemake.output.vamb_tax)
