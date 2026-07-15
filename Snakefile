import yaml
from pathlib import Path

configfile: "config.yaml"

OUT = Path(config["outdir"])
ENVS = config["envs"]
DB = config["db"]
THR = config["thresholds"]
SAMPLES = list(config["samples"].keys())
THREADS = config["threads"]

def env(name):
    return ENVS[name]

def run_in(name, cmd):
    return f"conda run --no-capture-output -p {env(name)} {cmd}"

rule all:
    input:
        expand(str(OUT / "{sample}" / "checkm2_1" / "quality_report.tsv"), sample=SAMPLES),
        expand(str(OUT / "{sample}" / "rrna" / "rrna_presence.tsv"), sample=SAMPLES),
        expand(str(OUT / "{sample}" / "dastool" / "dastool_DASTool_summary.tsv"), sample=SAMPLES)

# ------------------------------------------------------------------
# 1. Assembly: Flye
# ------------------------------------------------------------------
rule flye_assembly:
    input:
        reads=lambda wc: config["samples"][wc.sample]
    output:
        asm=str(OUT / "{sample}" / "assembly" / "assembly.fasta"),
        graph=str(OUT / "{sample}" / "assembly" / "assembly_graph.gfa"),
        info=str(OUT / "{sample}" / "assembly" / "assembly_info.txt")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "assembly"),
        mode=config["flye"]["read_mode"]
    threads: THREADS
    shell:
        run_in("flye", "flye {params.mode} {input.reads} --out-dir {params.outdir} "
                        "--threads {threads} --meta")

# ------------------------------------------------------------------
# 2. Polish: Medaka
# ------------------------------------------------------------------
rule medaka_polish:
    input:
        asm=str(OUT / "{sample}" / "assembly" / "assembly.fasta"),
        reads=lambda wc: config["samples"][wc.sample]
    output:
        polished=str(OUT / "{sample}" / "polish" / "consensus.fasta")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "polish")
    threads: THREADS
    shell:
        run_in("flye", "medaka_consensus -i {input.reads} -d {input.asm} "
                        "-o {params.outdir} -t {threads}")

# ------------------------------------------------------------------
# 3. Cobertura compartilhada: minimap2 -> BAM ordenado (reusado por todos os binners)
# ------------------------------------------------------------------
rule map_reads:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta"),
        reads=lambda wc: config["samples"][wc.sample]
    output:
        bam=str(OUT / "{sample}" / "mapping" / "reads.sorted.bam"),
        bai=str(OUT / "{sample}" / "mapping" / "reads.sorted.bam.bai")
    threads: THREADS
    shell:
        run_in("flye",
            "bash -c 'minimap2 -ax map-ont -t {threads} {input.asm} {input.reads} "
            "| samtools sort -@ {threads} -o {output.bam} - "
            "&& samtools index {output.bam}'")

rule depth_jgi:
    input:
        bam=str(OUT / "{sample}" / "mapping" / "reads.sorted.bam")
    output:
        depth=str(OUT / "{sample}" / "mapping" / "depth.txt")
    shell:
        run_in("binning", "jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}")

# ------------------------------------------------------------------
# 4a. MetaBAT2
# ------------------------------------------------------------------
rule metabat2:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta"),
        depth=str(OUT / "{sample}" / "mapping" / "depth.txt")
    output:
        done=str(OUT / "{sample}" / "binning" / "metabat2" / "DONE")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "binning" / "metabat2")
    threads: THREADS
    shell:
        "mkdir -p {params.outdir} && (" + \
        run_in("binning", "metabat2 -i {input.asm} -a {input.depth} "
                           "-o {params.outdir}/bin -m 1500 -t {threads}") + \
        " || true) ; touch {output.done}"

# ------------------------------------------------------------------
# 4b. MaxBin2
# ------------------------------------------------------------------
rule maxbin2:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta"),
        depth=str(OUT / "{sample}" / "mapping" / "depth.txt")
    output:
        done=str(OUT / "{sample}" / "binning" / "maxbin2" / "DONE")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "binning" / "maxbin2"),
        abund=lambda wc: str(OUT / wc.sample / "binning" / "maxbin2" / "abund.tsv")
    threads: THREADS
    shell:
        "mkdir -p {params.outdir} && "
        "cut -f1,3 {input.depth} | tail -n +2 > {params.abund} ; (" + \
        run_in("binning", "run_MaxBin.pl -contig {input.asm} -abund {params.abund} "
                           "-out {params.outdir}/bin -thread {threads}") + \
        " || true) ; touch {output.done}"

# ------------------------------------------------------------------
# 4c. CONCOCT
# ------------------------------------------------------------------
rule concoct:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta"),
        bam=str(OUT / "{sample}" / "mapping" / "reads.sorted.bam")
    output:
        clusters=str(OUT / "{sample}" / "binning" / "concoct" / "clustering_merged.csv")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "binning" / "concoct")
    threads: THREADS
    shell:
        "mkdir -p {params.outdir} && (" + \
        run_in("binning",
            "bash -c '"
            "cut_up_fasta.py {input.asm} -c 10000 -o 0 --merge_last -b {params.outdir}/contigs_10K.bed > {params.outdir}/contigs_10K.fa "
            "&& concoct_coverage_table.py {params.outdir}/contigs_10K.bed {input.bam} > {params.outdir}/coverage_table.tsv "
            "&& concoct --composition_file {params.outdir}/contigs_10K.fa --coverage_file {params.outdir}/coverage_table.tsv -b {params.outdir}/ -t {threads} "
            "&& merge_cutup_clustering.py {params.outdir}/clustering_gt1000.csv > {output.clusters}'") + \
        " || true) ; touch {output.clusters}"

# ------------------------------------------------------------------
# 4d. VAMB (default) + TaxVAMB
# ------------------------------------------------------------------
rule vamb_default:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta"),
        bam=str(OUT / "{sample}" / "mapping" / "reads.sorted.bam")
    output:
        clusters=str(OUT / "{sample}" / "binning" / "vamb" / "vae_clusters_unsplit.tsv")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "binning" / "vamb"),
        bamdir=lambda wc: str(OUT / wc.sample / "mapping")
    threads: THREADS
    shell:
        "rm -rf {params.outdir} ; (" + \
        run_in("vamb", "vamb bin default --outdir {params.outdir} --fasta {input.asm} "
                       "--bamdir {params.bamdir} -m 1500 -p {threads} -o ''") + \
        " || true) ; mkdir -p {params.outdir} ; touch {output.clusters}"

rule mmseqs_taxonomy:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta")
    output:
        tax=str(OUT / "{sample}" / "taxonomy" / "mmseqs_lineage.tsv")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "taxonomy"),
        gtdb_db=config["db"]["mmseqs2_gtdb"]
    threads: THREADS
    shell:
        "mkdir -p {params.outdir} && (" + \
        run_in("mmseqs2",
            "mmseqs easy-taxonomy {input.asm} {params.gtdb_db} "
            "{params.outdir}/mmseqs {params.outdir}/tmp --threads {threads} "
            "--tax-lineage 1 --search-type 3") + \
        " && cp {params.outdir}/mmseqs_lca.tsv {output.tax}) || true ; touch {output.tax}"

rule convert_taxonomy_for_vamb:
    input:
        tax=str(OUT / "{sample}" / "taxonomy" / "mmseqs_lineage.tsv")
    output:
        vamb_tax=str(OUT / "{sample}" / "taxonomy" / "vamb_taxonomy.tsv")
    script:
        "scripts/convert_mmseqs_to_vamb_taxonomy.py"

rule taxvamb:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta"),
        tax=str(OUT / "{sample}" / "taxonomy" / "vamb_taxonomy.tsv")
    output:
        clusters=str(OUT / "{sample}" / "binning" / "taxvamb" / "vaevae_clusters_unsplit.tsv")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "binning" / "taxvamb"),
        bamdir=lambda wc: str(OUT / wc.sample / "mapping")
    threads: THREADS
    shell:
        "rm -rf {params.outdir} ; (" + \
        run_in("vamb", "vamb bin taxvamb --outdir {params.outdir} --fasta {input.asm} "
                       "--bamdir {params.bamdir} --taxonomy {input.tax} -m 1500 -p {threads} -o ''") + \
        " || true) ; mkdir -p {params.outdir} ; touch {output.clusters}"

# ------------------------------------------------------------------
# 4e. SemiBin2
# ------------------------------------------------------------------
rule semibin2:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta"),
        bam=str(OUT / "{sample}" / "mapping" / "reads.sorted.bam")
    output:
        done=str(OUT / "{sample}" / "binning" / "semibin2" / "DONE")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "binning" / "semibin2")
    threads: THREADS
    shell:
        "mkdir -p {params.outdir} ; (" + \
        run_in("semibin", "SemiBin2 single_easy_bin -i {input.asm} -b {input.bam} "
                          "-o {params.outdir} --sequencing-type long_read -t {threads}") + \
        " || true) ; touch {output.done}"

# ------------------------------------------------------------------
# 4f. MetaCoAG (usa grafo do Flye)
# ------------------------------------------------------------------
rule metacoag:
    input:
        asm=str(OUT / "{sample}" / "assembly" / "assembly.fasta"),
        graph=str(OUT / "{sample}" / "assembly" / "assembly_graph.gfa"),
        info=str(OUT / "{sample}" / "assembly" / "assembly_info.txt"),
        depth=str(OUT / "{sample}" / "mapping" / "depth.txt")
    output:
        done=str(OUT / "{sample}" / "binning" / "metacoag" / "DONE")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "binning" / "metacoag"),
        abund=lambda wc: str(OUT / wc.sample / "binning" / "metacoag" / "abund.tsv")
    threads: THREADS
    shell:
        "mkdir -p {params.outdir} && "
        "cut -f1,3 {input.depth} | tail -n +2 > {params.abund} ; (" + \
        run_in("metacoag",
            "metacoag --assembler flye --graph {input.graph} --contigs {input.asm} "
            "--paths {input.info} --abundance {params.abund} --output {params.outdir} "
            "--min_bin_size 50000") + \
        " || true) ; touch {output.done}"

# ------------------------------------------------------------------
# 4g. LRBinner (modo contigs)
# ------------------------------------------------------------------
rule lrbinner:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta"),
        reads=lambda wc: config["samples"][wc.sample]
    output:
        done=str(OUT / "{sample}" / "binning" / "lrbinner" / "DONE")
    params:
        outdir=lambda wc: str(OUT / wc.sample / "binning" / "lrbinner"),
        repo=config["lrbinner_repo"],
        reads_plain=lambda wc: str(OUT / wc.sample / "binning" / "lrbinner_reads.fastq")
    threads: THREADS
    shell:
        # LRBinner so detecta o tipo de arquivo pela extensao literal (nao entende .gz)
        "mkdir -p {params.outdir} && "
        "(zcat {input.reads} > {params.reads_plain} 2>/dev/null || cp {input.reads} {params.reads_plain}) && " + \
        run_in("lrbinner", "python {params.repo}/lrbinner.py contigs "
                           "-r {params.reads_plain} -c {input.asm} -o {params.outdir} "
                           "--separate -t {threads}") + \
        " ; touch {output.done}"

# ------------------------------------------------------------------
# 5. Padronizar todos os binners em contig2bin.tsv (formato DAS Tool)
# ------------------------------------------------------------------
rule contig2bin_metabat2:
    input: str(OUT / "{sample}" / "binning" / "metabat2" / "DONE")
    output: str(OUT / "{sample}" / "contig2bin" / "metabat2.tsv")
    params: bindir=lambda wc: str(OUT / wc.sample / "binning" / "metabat2")
    shell:
        "(" + run_in("dastool", "Fasta_to_Contig2Bin.sh -e fa -i {params.bindir}") + " > {output}) || true ; touch {output}"

rule contig2bin_maxbin2:
    input: str(OUT / "{sample}" / "binning" / "maxbin2" / "DONE")
    output: str(OUT / "{sample}" / "contig2bin" / "maxbin2.tsv")
    params: bindir=lambda wc: str(OUT / wc.sample / "binning" / "maxbin2")
    shell:
        "(" + run_in("dastool", "Fasta_to_Contig2Bin.sh -e fasta -i {params.bindir}") + " > {output}) || true ; touch {output}"

rule contig2bin_semibin2:
    input: str(OUT / "{sample}" / "binning" / "semibin2" / "DONE")
    output: str(OUT / "{sample}" / "contig2bin" / "semibin2.tsv")
    params: bindir=lambda wc: str(OUT / wc.sample / "binning" / "semibin2" / "output_bins")
    shell:
        "(" + run_in("dastool", "Fasta_to_Contig2Bin.sh -e gz -i {params.bindir}") + " > {output}) || true ; touch {output}"

rule contig2bin_metacoag:
    input: str(OUT / "{sample}" / "binning" / "metacoag" / "DONE")
    output: str(OUT / "{sample}" / "contig2bin" / "metacoag.tsv")
    params: bindir=lambda wc: str(OUT / wc.sample / "binning" / "metacoag" / "bins")
    shell:
        "(" + run_in("dastool", "Fasta_to_Contig2Bin.sh -e fasta -i {params.bindir}") + " > {output}) || true ; touch {output}"

rule contig2bin_lrbinner:
    input: str(OUT / "{sample}" / "binning" / "lrbinner" / "DONE")
    output: str(OUT / "{sample}" / "contig2bin" / "lrbinner.tsv")
    params: bindir=lambda wc: str(OUT / wc.sample / "binning" / "lrbinner" / "binned_contigs")
    shell:
        "(" + run_in("dastool", "Fasta_to_Contig2Bin.sh -e fasta -i {params.bindir}") + " > {output}) || true ; touch {output}"

rule contig2bin_concoct:
    input: str(OUT / "{sample}" / "binning" / "concoct" / "clustering_merged.csv")
    output: str(OUT / "{sample}" / "contig2bin" / "concoct.tsv")
    shell:
        "awk -F',' 'NR>1{{print $1\"\\t\"$2}}' {input} > {output}"

rule contig2bin_vamb:
    input: str(OUT / "{sample}" / "binning" / "vamb" / "vae_clusters_unsplit.tsv")
    output: str(OUT / "{sample}" / "contig2bin" / "vamb.tsv")
    shell:
        "awk -F'\\t' 'NR>1{{print $2\"\\t\"$1}}' {input} > {output}"

rule contig2bin_taxvamb:
    input: str(OUT / "{sample}" / "binning" / "taxvamb" / "vaevae_clusters_unsplit.tsv")
    output: str(OUT / "{sample}" / "contig2bin" / "taxvamb.tsv")
    shell:
        "awk -F'\\t' 'NR>1{{print $2\"\\t\"$1}}' {input} > {output}"

# ------------------------------------------------------------------
# 6. DAS Tool: consolida os binners (7 de 8 -- taxvamb excluido)
# ------------------------------------------------------------------
# taxvamb precisa de um banco mmseqs2/GTDB completo (config db.mmseqs2_gtdb).
# O download do GTDB completo travou/crashou em rodadas anteriores deste
# projeto (fora deste repo) -- por decisao explicita, essa etapa fica de fora
# por enquanto (GTDB-Tk roda depois, separadamente, num sistema online).
# Como BINNER_NAMES controla o expand() de inputs do DAS Tool, remover
# "taxvamb" daqui desliga toda a cadeia mmseqs_taxonomy -> convert_taxonomy_for_vamb
# -> taxvamb do DAG (nada mais depende dela), sem precisar guardar/comentar
# as rules individualmente.
BINNER_NAMES = ["metabat2", "maxbin2", "concoct", "vamb",
                "semibin2", "metacoag", "lrbinner"]

rule dastool:
    input:
        asm=str(OUT / "{sample}" / "polish" / "consensus.fasta"),
        c2b=expand(str(OUT / "{{sample}}" / "contig2bin" / "{binner}.tsv"),
                    binner=BINNER_NAMES)
    output:
        summary=str(OUT / "{sample}" / "dastool" / "dastool_DASTool_summary.tsv")
    params:
        outbase=lambda wc: str(OUT / wc.sample / "dastool" / "dastool"),
        binner_names=BINNER_NAMES
    threads: THREADS
    run:
        import subprocess
        non_empty = [(name, path) for name, path in zip(params.binner_names, input.c2b)
                     if Path(path).stat().st_size > 0]
        skipped = [name for name, path in zip(params.binner_names, input.c2b)
                   if Path(path).stat().st_size == 0]
        if skipped:
            print(f"AVISO: binners sem bins (pulados no DAS Tool): {', '.join(skipped)}")
        if not non_empty:
            raise RuntimeError("Nenhum binner produziu bins -- nao ha nada para o DAS Tool consolidar.")
        c2b_list = ",".join(path for _, path in non_empty)
        labels = ",".join(name for name, _ in non_empty)
        outdir = Path(params.outbase).parent
        outdir.mkdir(parents=True, exist_ok=True)
        score_threshold = config.get("dastool_score_threshold", 0.5)
        cmd = run_in("dastool",
            f"DAS_Tool -i {c2b_list} -l {labels} -c {input.asm} -o {params.outbase} "
            f"--search_engine diamond --write_bins --threads {threads} "
            f"--score_threshold {score_threshold}")
        subprocess.run(cmd, shell=True, check=True)

# ------------------------------------------------------------------
# 7. CheckM2 #1 (filtro inicial: so vale a pena limpar quem ja bate 50%/10%)
# ------------------------------------------------------------------
rule checkm2_first_pass:
    input:
        summary=str(OUT / "{sample}" / "dastool" / "dastool_DASTool_summary.tsv")
    output:
        report=str(OUT / "{sample}" / "checkm2_1" / "quality_report.tsv")
    params:
        bindir=lambda wc: str(OUT / wc.sample / "dastool" / "dastool_DASTool_bins"),
        outdir=lambda wc: str(OUT / wc.sample / "checkm2_1"),
        db=DB["checkm2"]
    threads: THREADS
    shell:
        "CHECKM2DB={params.db} " + \
        run_in("checkm2", "checkm2 predict -i {params.bindir} -x fa "
                          "-o {params.outdir} -t {threads} --force")

# ------------------------------------------------------------------
# 8. barrnap: presenca de rRNA (campo exigido pelo ENA/genome_uploader)
# ------------------------------------------------------------------
rule barrnap_check:
    input:
        report=str(OUT / "{sample}" / "checkm2_1" / "quality_report.tsv")
    output:
        tsv=str(OUT / "{sample}" / "rrna" / "rrna_presence.tsv")
    params:
        bindir=lambda wc: str(OUT / wc.sample / "dastool" / "dastool_DASTool_bins"),
        outdir=lambda wc: str(OUT / wc.sample / "rrna")
    shell:
        "mkdir -p {params.outdir} && "
        "echo -e 'bin_id\\thas_16S\\thas_23S\\thas_5S' > {output.tsv} && "
        "for f in {params.bindir}/*.fa; do "
        "  bid=$(basename $f .fa); "
        "  " + run_in("barrnap", "barrnap --quiet $f") + " > {params.outdir}/$bid.gff 2>/dev/null; "
        "  h16=$(grep -c '16S' {params.outdir}/$bid.gff || true); "
        "  h23=$(grep -c '23S' {params.outdir}/$bid.gff || true); "
        "  h5=$(grep -c '\\b5S' {params.outdir}/$bid.gff || true); "
        "  echo -e \"$bid\\t$h16\\t$h23\\t$h5\" >> {output.tsv}; "
        "done"
