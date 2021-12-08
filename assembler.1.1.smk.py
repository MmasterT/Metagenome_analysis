# -*- coding: utf-8 -*-

#runs in snakemake version 5.8.2

import os
#Set location of the configfile.

configfile: "/u02/Mofedi/Test/config.json"
#Set the snakefile working directory.

wordkir: config["workdir"]

#Create list with the wildcards to use.
SAMPLES = []
SENSE = ["1", "2"]
ASSAMBLER = ["spades", "megahit"]
BINNER = ["metabat", "maxbin", "vamb"]
EXTENSION = []
VAMB = range(0,4)

#Get the filename of te samples

file_names = [f for f in os.listdir(os.path.join(config["workdir"], "data/raw_reads"))
              if os.path.isfile(os.path.join(config["workdir"] + "/data/raw_reads", f))]

#From the filename get the sample name
for name in file_names:
    
    section_1= name.split(".")
    section_1.pop(-1)
    file = []
    for i in section_1:
        section_2 = i.split("_")
        for j in section_2:
            file.append(j)
    
    samp = file[0]
    SAMPLES.append(samp)
    
    with open(config['workdir'] + "/data/raw_reads/" + name, 'r') as f:
        first_line = f.readline().strip()
        sense = first_line.replace('1:N:0','')
        if first_line == sense:
            sense = 2
        else:
            sense = 1
        file_name = str(samp) + "_" + str(sense) + ".fastq"
        file_dir = os.path.join(config["workdir"], "data/raw_reads")
        os.rename(os.path.join(file_dir, name),os.path.join(file_dir, file_name))
        f.close()
SAMPLES = list(set(SAMPLES))

if config['assembler_method'] == 'concated':
    
    pair_1 = "data/procesed_reads/reads_concatfoward.fastq"
    pair_2 = "data/procesed_reads/reads_concatreverse.fastq"
    out_spades = "data/contigs_spades/contigs.fasta"
    outdir_spades = "data/contigs_spades/"
    out_megahit = "data/contigs_megahit/final.contigs.fa"
    outdir_megahit= "data/temp"
    finaldir_megahit = "data/contigs_megahit/"
    move_1 = "data/contigs_megahit/final.contigs.fa"
    move_2 = "data/contigs_megahit/contigs.fasta"
    all = "data/contigs_{assembler}/contigs.fasta"

else:
    pair_1 =  "data/procesed_reads/{sample}_1_final.fastq"
    pair_2 =  "data/procesed_reads/{sample}_2_final.fastq"
    out_spades = "data/contigs_spades_{sample}/contigs.fasta"
    outdir_spades = "data/contigs_spades_{sample}/"
    out_megahit = "data/contigs_megahit_{sample}/final.contigs.fa"
    outdir_megahit= "data/temp_{sample}"
    finaldir_megahit = "data/contigs_megahit_{sample}/"
    move_1 = "data/contigs_megahit_{sample}/final.contigs.fa"
    move_2 = "data/contigs_megahit_{sample}/contigs.fasta"
    all = "data/contigs_{assembler}_{sample}/contigs.fasta"


rule all:
    input:
        all = expand(all, assembler = ASSAMBLER)
        
rule move_contigs:
    input:
        move_1 = move_1
    output:
        move_2 = move_2
    shell:
        """
        mv {input.move_1} {output.move_2}
        """

#Runs megahit. Snakemake creates the output folder and megahit cont overwrite data so you have
#to run a command that creates a temporal output and then change it to the final name.
rule megahit:
    input:
        pair_1 = pair_1, pair_2 = pair_2
    output:
        out_megahit = out_megahit
    params:
        outdir_megahit = outdir_megahit, finaldir_megahit = finaldir_megahit, walltime="86400", node = "2",
        kmer = config["k-list"], mem = config["mega_mem"]
    threads:
        int(config["megahit_threads"])
    conda:
        "envs/megahit.yaml"
    shell:
        """
        rm -r -f {params.outdir_megahit}
        touch {output.out_megahit}
        megahit -1 {input.pair_1} -2 {input.pair_2} --k-list {params.kmer} -m {params.mem} -t {threads} -o {params.outdir_megahit}                                                                                                                                             
        rm -r -f {params.finaldir_megahit}                                                                                                                             
        mv {params.outdir_megahit} {params.finaldir_megahit}                                                                                                                                   
        """
#Runs metaSpades.      
rule spades:
    input:
        pair_1 = pair_1, pair_2 = pair_2
    output:
        out_spades = out_spades
    params:
        outdir_spades = outdir_spades, walltime="86400", node = "1",
        kmer = config["k-list"],mem = config["spades_mem"]
    threads:
        int(config["spades_threads"])
    conda:
        "envs/spades.yaml"
    shell:
        """
        rm -r -f {params.outdir_spades}
        spades.py --meta -1 {input.pair_1} -2 {input.pair_2} -t {threads} -m {params.mem} -k {params.kmer} -o {params.outdir_spades}
        """
rule concat_pairs:
    input:
        expand("data/procesed_reads/{{sample}}_{sense}_final.fastq", sense = SENSE)        
    output:
        temp("data/procesed_reads/{sample}_merged_reads.fastq")
    shell:
        "cat {input} > {output}"

rule concat_sample:
    input:
        read_1 = expand("data/procesed_reads/{sample}_1_final.fastq", sample=set(SAMPLES)),
        read_2 = expand("data/procesed_reads/{sample}_2_final.fastq", sample=set(SAMPLES))
    output:
        output_1 = "data/procesed_reads/reads_concatfoward.fastq",
        output_2 = "data/procesed_reads/reads_concatreverse.fastq"
    shell:
        "cat {input.read_1} > {output.output_1};"
        "cat {input.read_2} > {output.output_2}"

#Runs custom script. Filter the reads with the specifications of the config file. It checks
#that all reads has its pair.
rule filter_reads:
    input:
        fow_2 = "data/procesed_reads/{sample}_1trim.fastq",
        rev_2 = "data/procesed_reads/{sample}_2trim.fastq"
    output:
        filter_1 = "data/procesed_reads/{sample}_1_final.fastq",
        filter_2 = "data/procesed_reads/{sample}_2_final.fastq"
    script:
        "scripts/Q_Filter.4.py"

#Runs custom script. Cut the reads with the specifications of the config.json file.
rule trim:
    input:
        fow_1 = "data/raw_reads/{sample}_{sense}.fastq"
    output:
        out_1 = temp("data/procesed_reads/{sample}_{sense}trim.fastq")
    script:
        "scripts/triming.py"
