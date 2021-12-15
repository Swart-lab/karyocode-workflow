# If rules are split across multiple snakefiles, list them here
# include: "rules-A"
# include: "rules-B"
PHYLOFLASH_DBHOME="/ebio/abt2_projects/ag-swart-loxodes/db/phyloFlash/132"

rule all:
    input:
        expand("qc/phyloFlash/{lib}_rnaseq.phyloFlash.tar.gz", lib=config['rnaseq'])

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
