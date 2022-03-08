# If rules are split across multiple snakefiles, list them here
# include: "rules-A"
# include: "rules-B"
PHYLOFLASH_DBHOME="/ebio/abt2_projects/ag-swart-loxodes/db/phyloFlash/132"

rule all:
    input:
        expand(
            "qc/phyloFlash/{lib}_rnaseq_{readlim}.phyloFlash.tar.gz",
            lib=[i for i in config['rnaseq'] if i not in ['kn01_yy','cmag_mm','pard_mm','bjap_mm','fsal_mm']], readlim=[1000000]), # skip libraries where SPAdes step fails
        expand("qc/phyloFlash/{lib}_rnaseq_maponly_1000000.phyloFlash.tar.gz", lib=config['rnaseq']),
        expand("qc/phyloFlash/{lib}_rnaseq_se_maponly_1000000.phyloFlash.tar.gz", lib=config['rnaseq_se']),
        expand("qc/phyloFlash/{lib}_dnaseq.phyloFlash.tar.gz", lib=config['mda']),
        expand("assembly/trinity.{lib}.Trinity.fasta", lib=config['rnaseq']),
        expand("qc/trinity_assem/trinity.{lib}.v.{dbprefix}.blastx.out6.w_pct_hit_length", lib=config['rnaseq'], dbprefix=['uniprot_sprot','bsto_mac']),
        "assembly/trinity_gg.bsto.Trinity-GG.fasta", # genome-guided assembly
        "assembly/trinity_se.bsp_mk.Trinity.fasta",    # single-end read library
        "qc/trinity_assem/trinity_se.bsp_mk.v.bsto_mac.blastx.out6.w_pct_hit_length",    # single-end read library

rule phyloflash_rnaseq:
    input:
        fwd=lambda wildcards: config['rnaseq'][wildcards.lib]['fwd'],
        rev=lambda wildcards: config['rnaseq'][wildcards.lib]['rev']
    output:
        "qc/phyloFlash/{lib}_rnaseq_{readlim}.phyloFlash.tar.gz"
    conda: "envs/phyloflash.yml"
    threads: 16
    log: "logs/phyloflash_rnaseq.{lib}.{readlim}.log"
    wildcard_constraints: readlim="\d+"
    params:
        db=PHYLOFLASH_DBHOME
    shell:
        # Use -readlimit 1000000 for transcriptome libraries
        "phyloFlash.pl -lib {wildcards.lib}_rnaseq_{wildcards.readlim} -readlength 150 -readlimit {wildcards.readlim} -read1 {input.fwd} -read2 {input.rev} -CPUs {threads} -almosteverything -dbhome {params.db} 2> {log};"
        "mv {wildcards.lib}_rnaseq_{wildcards.readlim}.phyloFlash* qc/phyloFlash/;"

rule phyloflash_rnaseq_se:
    input:
        lambda wildcards: config['rnaseq_se'][wildcards.lib]
    output:
        "qc/phyloFlash/{lib}_rnaseq_se_{readlim}.phyloFlash.tar.gz"
    conda: "envs/phyloflash.yml"
    threads: 16
    wildcard_constraints: readlim="\d+"
    log: "logs/phyloflash_rnaseq_se.{lib}.{readlim}.log"
    params:
        db=PHYLOFLASH_DBHOME
    shell:
        # Do not run -poscov for SE libraries
        # Use -readlimit 1000000 for transcriptome libraries
        "phyloFlash.pl -lib {wildcards.lib}_rnaseq_se_{wildcards.readlim} -readlength 100 -readlimit {wildcards.readlim} -read1 {input} -CPUs {threads} -html -zip -treemap -log -dbhome {params.db} 2> {log};"
        "mv {wildcards.lib}_rnaseq_se_{wildcards.readlim}.phyloFlash* qc/phyloFlash/;"

rule phyloflash_rnaseq_maponly:
    input:
        fwd=lambda wildcards: config['rnaseq'][wildcards.lib]['fwd'],
        rev=lambda wildcards: config['rnaseq'][wildcards.lib]['rev']
    output:
        "qc/phyloFlash/{lib}_rnaseq_maponly_{readlim}.phyloFlash.tar.gz"
    conda: "envs/phyloflash.yml"
    threads: 16
    wildcard_constraints: readlim="\d+"
    log: "logs/phyloflash_rnaseq_maponly.{lib}.{readlim}.log"
    params:
        db=PHYLOFLASH_DBHOME
    shell:
        "phyloFlash.pl -lib {wildcards.lib}_rnaseq_maponly_{wildcards.readlim} -readlength 150 -readlimit {wildcards.readlim} -skip_spades -read1 {input.fwd} -read2 {input.rev} -CPUs {threads} -html -treemap -log -poscov -zip -dbhome {params.db} 2> {log};"
        "mv {wildcards.lib}_rnaseq_maponly_{wildcards.readlim}.phyloFlash* qc/phyloFlash/;"

rule phyloflash_rnaseq_se_maponly:
    input:
        lambda wildcards: config['rnaseq_se'][wildcards.lib]
    output:
        "qc/phyloFlash/{lib}_rnaseq_se_maponly_{readlim}.phyloFlash.tar.gz"
    conda: "envs/phyloflash.yml"
    threads: 16
    wildcard_constraints: readlim="\d+"
    log: "logs/phyloflash_rnaseq_se_maponly.{lib}.{readlim}.log"
    params:
        db=PHYLOFLASH_DBHOME
    shell:
        # Do not run -poscov for SE libraries
        # Use -readlimit 1000000 for transcriptome libraries
        "phyloFlash.pl -lib {wildcards.lib}_rnaseq_se_maponly_{wildcards.readlim} -readlength 100 -readlimit {wildcards.readlim} -skip_spades -read1 {input} -CPUs {threads} -html -treemap -log -zip -dbhome {params.db} 2> {log};"
        "mv {wildcards.lib}_rnaseq_se_maponly_{wildcards.readlim}.phyloFlash* qc/phyloFlash/;"

rule phyloflash_dnaseq:
    input:
        fwd=lambda wildcards: config['mda'][wildcards.lib]['fwd'],
        rev=lambda wildcards: config['mda'][wildcards.lib]['rev']
    output:
        "qc/phyloFlash/{lib}_dnaseq.phyloFlash.tar.gz"
    conda: "envs/phyloflash.yml"
    threads: 16
    log: "logs/phyloflash_dnaseq.{lib}.log"
    params:
        db=PHYLOFLASH_DBHOME
    shell:
        "phyloFlash.pl -lib {wildcards.lib}_dnaseq -readlength 150 -read1 {input.fwd} -read2 {input.rev} -CPUs {threads} -almosteverything -dbhome {params.db} 2> {log};"
        "mv {wildcards.lib}_dnaseq.phyloFlash* qc/phyloFlash/;"

rule trim_reads_rnaseq:
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

rule trim_reads_dnaseq:
    input:
        fwd=lambda wildcards: config['dnaseq'][wildcards.lib]['fwd'],
        rev=lambda wildcards: config['dnaseq'][wildcards.lib]['rev']
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


rule trinity_fulllength_qc:
    # based on https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts
    # TODO use ciliate proteoms as reference db
    input:
        assem="assembly/{prefix}.Trinity.fasta",
        db="db/{dbprefix}.fasta"
    output:
        blastx="qc/trinity_assem/{prefix}.v.{dbprefix}.blastx.out6",
        pchitlen="qc/trinity_assem/{prefix}.v.{dbprefix}.blastx.out6.w_pct_hit_length"
    threads: 8
    conda: "envs/trinity.yml"
    log: "logs/trinity_fulllength_qc.{dbprefix}.{prefix}.log"
    params:
        db_prefix="db/{dbprefix}"
    shell:
        "cat {input.assem} | parallel --gnu -j {threads} --recstart '>' -N 100 --pipe blastx -query - -db {params.db_prefix} -evalue 1e-20 -max_target_seqs 1 -outfmt 6 > {output.blastx};"
        "analyze_blastPlus_topHit_coverage.pl {output.blastx} {input.assem} {input.db} &> {log};"


rule trim_reads_rnaseq_se: # Single-end sequencing
    input:
        lambda wildcards: config['rnaseq_se'][wildcards.lib]
    output:
        "data/reads_trim/{lib}.trim.q24.U.fastq.gz"
    threads: 8
    log: "logs/trim_reads_se.{lib}.log"
    conda: "envs/bbmap.yml"
    shell:
        "bbduk.sh -Xmx10g threads={threads} ref=resources/adapters.fa,resources/phix174_ill.ref.fa.gz in={input} ktrim=r qtrim=rl trimq=24 minlength=25 out={output} 2> {log}"

rule trinity_se:
    input:
        "data/reads_trim/{lib}.trim.q24.U.fastq.gz"
    output: "assembly/trinity_se.{lib}.Trinity.fasta"
    conda: "envs/trinity.yml"
    threads: 24
    params:
        prefix="assembly/trinity_se.{lib}"
    log: "logs/trinity_se.{lib}.log"
    shell:
        "Trinity --seqType fq --max_memory 64G --bflyHeapSpaceMax 40G --CPU {threads} --full_cleanup --single {input} --output {params.prefix} &> {log};"

rule trinity_gg:
    # Genome guided assembly from mapping file
    input:
        lambda wildcards: config['rnaseq_bam'][wildcards.lib]
    output: "assembly/trinity_gg.{lib}.Trinity-GG.fasta"
    conda: "envs/trinity.yml"
    threads: 24
    params:
        prefix="assembly/trinity_gg.{lib}"
    log: "logs/trinity_gg.{lib}.log"
    shell:
        # jaccard clip option for gene-dense genomes - without this option, we have about 1/3 of contigs with poly-A on both ends
        "Trinity --genome_guided_bam {input} --genome_guided_max_intron 50 --max_memory 64G --CPU {threads} --bflyHeapSpaceMax 40G --output {params.prefix} --jaccard_clip --full_cleanup &> {log};"
        "mv {params.prefix}/Trinity-GG.fasta {output};" # manual cleanup because --full_cleanup doesn't work with genome guided mode
        "mv {params.prefix}/Trinity-GG.fasta.gene_trans_map {output}.gene_trans_map;"
        # "rm -r {prefix}/;"
