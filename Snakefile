# If rules are split across multiple snakefiles, list them here
# include: "rules-A"
# include: "rules-B"
PHYLOFLASH_DBHOME="/ebio/abt2_projects/ag-swart-loxodes/db/phyloFlash/132"

rule all:
    input:
        expand("qc/phyloFlash/{lib}_rnaseq.phyloFlash.tar.gz", lib=config['rnaseq']),
        expand("assembly/trinity.{lib}.Trinity.fasta", lib=config['rnaseq'])

rule phyloflash_rnaseq:
    input:
        fwd=lambda wildcards: config['rnaseq'][wildcards.lib]['fwd'],
        rev=lambda wildcards: config['rnaseq'][wildcards.lib]['rev']
    output:
        "qc/phyloFlash/{lib}_rnaseq.phyloFlash.tar.gz"
    conda: "envs/phyloflash.yml"
    threads: 16
    log: "logs/phyloflash_rnaseq.{lib}.log"
    params:
        db=PHYLOFLASH_DBHOME
    shell:
        # Use -readlimit 1000000 for transcriptome libraries
        "phyloFlash.pl -lib {wildcards.lib}_rnaseq -readlength 150 -readlimit 1000000 -read1 {input.fwd} -read2 {input.rev} -CPUs {threads} -almosteverything -dbhome {params.db} 2> {log};"
        "mv {wildcards.lib}_rnaseq.phyloFlash* qc/phyloFlash/;"

rule trim_reads:
    input:
        fwd=lambda wildcards: config['rnaseq'][wildcards.lib]['fwd'],
        rev=lambda wildcards: config['rnaseq'][wildcards.lib]['rev']
    output:
        fwd="data/reads_trim/{lib}.trim.q24.R1.fastq.gz",
        rev="data/reads_trim/{lib}.trim.q24.R2.fastq.gz"
    threads: 8
    log: "logs/trim_reads.{lib}.log"
    conda: "envs/bbmap.yml"
    shell:
        "bbduk.sh -Xmx10g threads={threads} ref=resources/adapters.fa,resources/phix174_ill.ref.fa.gz in={input.fwd} in2={input.rev} ktrim=r qtrim=rl trimq=24 minlength=25 out={output.fwd} out2={output.rev} 2> {log}"


rule trinity:
    input:
        fwd="data/reads_trim/{lib}.trim.q24.R1.fastq.gz",
        rev="data/reads_trim/{lib}.trim.q24.R2.fastq.gz"
    output: "assembly/trinity.{lib}.Trinity.fasta"
    conda: "envs/trinity.yml"
    threads: 24
    params:
        prefix="assembly/trinity.{lib}"
    log: "logs/trinity.{lib}.log"
    shell:
        "Trinity --seqType fq --max_memory 64G --bflyHeapSpaceMax 40G --CPU {threads} --full_cleanup --left {input.fwd} --right {input.rev} --output {params.prefix} &> {log};"


