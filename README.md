# nanopore-mag-pipeline

Pipeline de recuperação de MAGs (metagenome-assembled genomes) a partir de
leituras Oxford Nanopore, com binning por consenso (8 binners), um loop de
descontaminação guiado por quimerismo (GUNC), e geração automática da tabela
de metadados exigida pelo ENA para submissão de MAGs.

Este repositório contém apenas o **código de orquestração** (Snakemake +
scripts Python). Ele chama, via ambientes conda separados, um conjunto de
ferramentas de terceiros — nenhum código-fonte delas está incluído aqui.
Veja [Ferramentas utilizadas](#ferramentas-utilizadas) para créditos completos.

## Arquitetura

```
reads.fastq(.gz)
      │
      ▼
   Flye (--nano-hq/--nano-raw, --meta)  ──► assembly.fasta + assembly_graph.gfa
      │
      ▼
   Medaka (polimento)  ──► consensus.fasta
      │
      ▼
   minimap2 + samtools + jgi_summarize_bam_contig_depths
      │  (BAM/depth compartilhado por todos os binners abaixo)
      ▼
   ┌─────────────────────────────────────────────────────────────┐
   │  MetaBAT2 · MaxBin2 · CONCOCT · VAMB · TaxVAMB (+mmseqs2)    │
   │  SemiBin2 · MetaCoAG (usa o grafo do Flye) · LRBinner        │
   └─────────────────────────────────────────────────────────────┘
      │  (cada binner convertido em contig2bin.tsv)
      ▼
   DAS Tool  ──► consolida os binners que produziram bins
      │           (tolerante: pula binners sem resultado)
      ▼
   CheckM2 #1  ──► filtro inicial (completude ≥50%, contaminação ≤10%)
      │
      ▼
   barrnap  ──► presença de rRNA (16S/23S/5S, campo exigido pelo ENA)
      │
      ▼
   Para cada bin aprovado no CheckM2 #1:
      GUNC (banco GTDB r95) ──► clade separation score (CSS)
         │
         ├─ CSS ≤ 0.45 (limpo) ──────────────────────────────► segue
         │
         └─ CSS > 0.45 (quimérico)
                │
                ▼
             Deepurify (limpador 1) ──► GUNC de novo
                │
                ├─ passou ──────────────────────────────────► segue
                │
                └─ ainda quimérico
                       ▼
                    MAGpurify2 (limpador 2) ──► GUNC de novo
                       │
                       ├─ passou ───────────────────────────► segue
                       └─ ainda quimérico ──► DESCARTADO
      │
      ▼
   CheckM2 #2 (final, pós-limpeza)  ──► reconfere 50%/10%
      │
      ▼
   MAG aprovado ──► tabela pronta para `genome_upload` (EBI-Metagenomics)
```

Toda decisão fica registrada em `audit_table.tsv` (uma linha por bin, com
os scores de CSS/CheckM2 em cada etapa e o motivo de aprovação/descarte).

## Uso

```bash
cd mag_pipeline
# 1. edite config.yaml: samples, outdir, e o banco de taxonomia do mmseqs2
snakemake -j <threads> -p --keep-going

# 2. loop de descontaminação (GUNC/Deepurify/MAGpurify2/CheckM2 final)
python scripts/gunc_decontam_loop.py \
  --sample-dir <outdir>/<amostra> \
  --gunc-env <env> --gunc-db <gunc_db_gtdb95.dmnd> \
  --deepurify-env <env> --deepurify-db <Deepurify-DB> \
  --magpurify-env <env> \
  --checkm2-env <env> --checkm2-db <checkm2_db.dmnd>

# 3. tabela final para submissão ENA
python scripts/build_ena_tsv.py --sample-dir <outdir>/<amostra> \
  --sample-accession <SAMEA/SRSxxxx>
```

Antes de rodar em produção: troque o banco de taxonomia do mmseqs2 (usado
só pelo TaxVAMB) de Kalamari — pequeno, usado só para validar a mecânica —
pelo GTDB completo (comando em `config.yaml`).

## Ferramentas utilizadas

Este pipeline não reimplementa nenhuma dessas ferramentas — apenas as chama.
Cite os trabalhos originais, não este repositório, ao reportar resultados
científicos gerados com estas ferramentas.

| Etapa | Ferramenta | Referência |
|---|---|---|
| Montagem | [Flye / metaFlye](https://github.com/mikolmogorov/Flye) | Kolmogorov et al. 2020, *Nat. Biotechnol.* |
| Polimento | [Medaka](https://github.com/nanoporetech/medaka) | Oxford Nanopore Technologies |
| Mapeamento | [minimap2](https://github.com/lh3/minimap2) | Li 2018, *Bioinformatics* |
| Binning | [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) | Kang et al. 2019, *PeerJ* |
| Binning | [MaxBin2](https://sourceforge.net/projects/maxbin2/) | Wu et al. 2016, *Bioinformatics* |
| Binning | [CONCOCT](https://github.com/BinPro/CONCOCT) | Alneberg et al. 2014, *Nat. Methods* |
| Binning | [VAMB / TaxVAMB](https://github.com/RasmussenLab/vamb) | Nissen et al. 2021; Kutuzova et al. 2024 |
| Binning | [SemiBin2](https://github.com/BigDataBiology/SemiBin) | Pan et al. 2023, *Bioinformatics* |
| Binning | [MetaCoAG](https://github.com/metagentools/MetaCoAG) | Mallawaarachchi & Lin 2022 |
| Binning | [LRBinner](https://github.com/anuradhawick/LRBinner) | Wickramarachchi & Lin 2021/2022 |
| Consenso de bins | [DAS Tool](https://github.com/cmks/DAS_Tool) | Sieber et al. 2018, *Nat. Microbiol.* |
| Quimerismo | [GUNC](https://github.com/grp-bork/gunc) | Orakov et al. 2021, *Genome Biol.* |
| Qualidade | [CheckM2](https://github.com/chklovski/CheckM2) | Chklovski et al. 2023, *Nat. Methods* |
| Descontaminação | [Deepurify](https://github.com/zoubohao/Deepurify) | Zou et al. 2024, *Nat. Mach. Intell.* |
| Descontaminação | [MAGpurify2](https://github.com/apcamargo/magpurify2) | Camargo et al. 2023, *ISME J.* |
| rRNA | [barrnap](https://github.com/tseemann/barrnap) | Torsten Seemann |
| Taxonomia (TaxVAMB) | [mmseqs2](https://github.com/soedinglab/MMseqs2) | Steinegger & Söding 2017; Mirdita et al. 2021 |
| Submissão ENA | [genome_uploader](https://github.com/EBI-Metagenomics/genome_uploader) | EBI-Metagenomics |
| Orquestração | [Snakemake](https://github.com/snakemake/snakemake) | Mölder et al. 2021, *F1000Research* |

## Licença

O código de orquestração neste repositório (Snakefile e scripts em `scripts/`)
está sob licença MIT — veja [LICENSE](LICENSE). Cada ferramenta chamada pelo
pipeline mantém sua própria licença; consulte o repositório de origem de cada
uma antes de redistribuir.
