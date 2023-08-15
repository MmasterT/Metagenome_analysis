#!/usr/bin/env python

# set configurations
ASSAMBLER = ["spades", "megahit"]
SAMPLES =   config["samples"]
INDEX_SIZE = config.get("index_size", "12G")
MM_MEM = config.get("minimap_mem", "35gb")
MM_PPN = config.get("minimap_ppn", "10")
VAMB_MEM = config.get("vamb_mem", "20gb")
VAMB_PPN = config.get("vamb_ppn", "10")
SAMPLE_DATA = contig("samples2data.txt")
CONTIGS = contig("contigs.txt")
VAMB_PARAMS = config.get("vamb_params", "-o C -m 2000 --minfasta 500000")
VAMB_PRELOAD = config.get("vamb_preload", "")
SENSE = [1,2]
# parse if GPUs is needed #
VAMB_split = VAMB_PPN.split(":") 
VAMB_threads = VAMB_split[0]
VAMB = range(0,4)



# read in list of per-sample assemblies

## start of snakemake rules ##
# targets
rule all:
    input:
        "jgi_matrix/jgi.abundance.dat",
        "binning/vamb/clusters.tsv"


rule cat_contigs:
    input:
        expand("data/contigs_{assembler}_{{sample}}/contigs.fasta", assembler=ASSEMBLER)
    output:
        "contigs.flt.fna.gz"
    params:
        walltime="864000", nodes="1", ppn="1", mem="10gb"
    threads:
        int(1)
    log:
        "log/contigs/catcontigs.log"
    shell:"""
        concatenate.py {output} {input} -m 2000
        """

rule index:
    input:
        contigs = "contigs.flt.fna.gz"
    output:
        mmi = "contigs.flt.mmi"
    params:
        walltime="864000",
        nodes="1",
        ppn="1",
        mem="90gb"
    threads:
        int(1)
    log:
        "log/contigs/index.log"
    shell:"""
        minimap2 -I {INDEX_SIZE} -d {output.mmi} {input.contigs} 2> {log}
        """

rule dict:
    input:
        contigs = "contigs.flt.fna.gz"
    output:
        dict = "contigs.flt.dict"
    params:
        walltime="864000",
        nodes="1",
        ppn="1",
        mem="10gb"
    threads:
        int(config["sam_threads"])
    log:
        "log/contigs/dict.log"
    shell:"""
        samtools dict {input.contigs} | cut -f1-3 > {output.dict} 2> {log}
        """

rule minimap:
    input:
        fq = ,
        mmi = "contigs.flt.mmi",
        dict = "contigs.flt.dict"
    output:
        bam = temp("mapped/{sample}.bam")
    params:
        walltime="864000",
        nodes="1",
        ppn=MM_PPN,
        mem=MM_MEM
    threads:
        int(MM_PPN)
    log:
        "log/map/{sample}.minimap.log"
    shell:"""
            minimap2 -t {threads} -ax sr {input.mmi} {input.fq} | \
            grep -v "^@" | cat {input.dict} - | \
            samtools view -F 3584 -b - > {output.bam} 2>{log}
        """

rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        temp("mapped/{sample}.bam.sort")
    params:
        walltime="864000",
        nodes="1",
        ppn="2",
        mem="15gb",
        prefix="mapped/tmp.{sample}"
    threads:
        int(config["sam_threads"])
    log:
        "log/map/{sample}.sort.log"
    shell:"""
        samtools sort {input} -T {params.prefix} --threads {threads} \ 
        -m 3G -o {output} 2>{log}
        """

rule jgi:
    input:
        bam = "mapped/{sample}.bam.sort"
    output:
        jgi = temp("jgi/{sample}.raw.jgi")
    params:
        walltime="864000",
        nodes="1",
        ppn="1",
        mem="10gb"
    threads:
        int(1)
    log:
        "log/jgi/{sample}.jgi"
    shell:"""
        jgi_summarize_bam_contig_depths --noIntraDepthVariance \
        --outputDepth {output} {input} 2>{log}
        """

rule cut_column1to3: 
    input:
        "jgi/%s.raw.jgi" % SAMPLES[0] 
    output:
        "jgi/jgi.column1to3"
    params:
        walltime="86400",
        nodes="1",
        ppn="1",
        mem="1gb"
    log:
        "log/jgi/column1to3"
    shell:""" 
        cut -f1-3 {input} > {output} 2>{log}
        """

rule cut_column4to5:
    input:
        "jgi/{sample}.raw.jgi"
    output:
        "jgi/{sample}.cut.jgi"
    params:
        walltime="86400",
        nodes="1",
        ppn="1",
        mem="1gb"
    log:
        "log/jgi/{sample}.cut.log"
    shell:"""
        cut -f1-3 --complement {input} > {output} 2>{log}
    """

rule paste_abundances:
    input:
        column1to3="jgi/jgi.column1to3",
        data=expand("jgi/{sample}.cut.jgi", sample= set(SAMPLES))
    output:
        abundance = "jgi_matrix/jgi.abundance.dat" 
    params:
        walltime="86400",
        nodes="1",
        ppn="1",
        mem="1gb"
    log:
        "log/jgi/paste_abundances.log"
    shell:"""
        paste {input.column1to3} {input.data} > {output.abundance} 2>{log}
        """

rule vamb:
    input:
        jgi = "jgi_matrix/jgi.abundance.dat",
        contigs = "contigs.flt.fna.gz"
    output:
        ""
    params:
        walltime = "86400",
        nodes = "1",
        ppn = VAMB_PPN,
        mem = VAMB_MEM,
        outdir = "binning/vamb"
    log:
        "log/vamb/vamb.log"
    threads:
        int(config["vamb_threads"])
    shell:"""
        rm -rf vamb
        vamb --outdir {params.outdir} --fasta {input.contigs} --jgi {input.jgi} \ 
        {VAMB_PARAMS} 2>{log} && touch {output.tch}
        """
