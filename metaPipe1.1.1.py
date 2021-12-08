## -*- coding: utf-8 -*-

import os
import pandas as pd
#Set location of the configfile.

configfile: "/u02/Mofedi/Test/config.json"

#Set the snakefile working directory. If not set the snakefile use the directory where it resides.

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
    name = name.rstrip()
    section = name.split("_", 1)    
    SAMPLES.append(section[0])
    
#getextension od the samples should be .fq or .fastq
for name in file_names:
    name = name.rstrip()
    section = name.split(".")
    EXTENSION.append(section[-1])
        
#no repeted elements to avoid repeted works in some rules

SAMPLES = list(set(SAMPLES))
#call the rules to make the contigs.


subworkflow vamb:
    configfile:
        "/u02/Mofedi/Test/config.json"
    workdir:
        config['workdir']
    snakefile:
        "vamb.snake.conda.py"

subworkflow contig:
    configfile:
        "/u02/Mofedi/Test/config.json"
    workdir:
        config['workdir']
    snakefile:
        "assembler.1.1.smk.py"

rule all:
    input:
        "drep/data_tables/Widb.csv"

rule drep:
    input:
        "done.txt"
    output:
        "drep/data_tables/Widb.csv"
    params:
        max_cont = config["max_contamination"], min_comp = config["min_completness"],
        sa = config ["sa"]
    conda:
        "envs/drep.yaml"
    threads:
        int(config["drep_thread"])
    shell:
        """
        rm -f done.txt
        dRep dereplicate -p {threads} -g genomes/* -comp {params.min_comp} -con {params.max_cont} -sa {params.sa} drep
        """

rule move_and_change_extension:
    input:
       vamb("vamb/clusters.tsv"),
        expand("metabat/bin_metabat_{assembler}.1.fa", assembler = ASSAMBLER),
        expand("maxbin/bin_maxbin_{assembler}.001.fasta", assembler = ASSAMBLER)
    output:
        "genomes/bin_vamb.0.fasta",
        expand("genomes/bin_maxbin_{assembler}.001.fasta", assembler = ASSAMBLER),
        expand("genomes/bin_metabat_{assembler}.1.fasta", assembler = ASSAMBLER),
        "done.txt"
    run:
        shell("mv maxbin/*.fasta genomes")
        
        directory = config['workdir'] #set directory
        final_dir = os.path.join(directory, "genomes")  #set the destination directory
        #get the path of the directory, the subdirectorys and all the file in it.
        
        for root, dirs, files in os.walk(directory):
            #iterate the subdirs
            for subdir in dirs:
                i=0
                for file in os.listdir(os.path.join(root,subdir)):
                    if file.endswith(".fa") and file.startswith("bin"):#if file is Metabat 2 bin
                        split  = file.split(".", 2)
                        split.pop()
                        name_no_extension = ".".join(split)
                        name_new_extension = name_no_extension + ".fasta" #change extension
                        destination = os.path.join(final_dir, name_new_extension)#set path
                        os.rename(os.path.join(root,subdir,file), destination)#move to destination
                    elif file.endswith(".fna") and file.startswith("S"):#if file VAMB bin
                        #do the same but number the not numbered VAMB files
                        name_new_extension = 'bin_vamb.' + str(i) + ".fasta"
                        destination = os.path.join(final_dir, name_new_extension)
                        os.rename(os.path.join(root,subdir,file), destination)
                        i += 1
                    else:
                        continue
        shell("touch done.txt")
        shell("rm -rf metabat")
        shell("rm -rf maxbin")
        
#runs metabat2
rule metabat2:
    input:
        abundance = "abundance/jgi_abundance_{assembler}.tsv"
    output:
        "metabat/bin_metabat_{assembler}.1.fa",
    params:
        contigs = "data/contigs_{assembler}/contigs.fasta", 
        walltime="86400"
    conda:
        "envs/metabat2.yaml"
    threads:
        int(config["meta_threads"])
    shell:
        """
        metabat -i {params.contigs} -t {threads} -a {input.abundance} -o metabat/bin_metabat_{wildcards.assembler} -m 2000
        """

#runs maxbin2       
rule maxbin2:
    input:
        abundance = "abundance/abund_{assembler}_maxbin_paths.txt"
    output:
        "maxbin/bin_maxbin_{assembler}.001.fasta"
    params:
        walltime="86400", mcl = config["min_contig_length"],
        contigs = "data/contigs_{assembler}/contigs.fasta"
    threads:
        int(config["maxb2_threads"])
    conda:
        "envs/maxbin.yaml"
    shell:
        """
        mkdir -p maxbin
        run_MaxBin.pl -contig {params.contigs} -out maxbin/bin_maxbin_{wildcards.assembler} -abund_list {input.abundance} -min_contig_length {params.mcl} -thread {threads}
        """

rule maxbin2_abundancepath:
    input: 
        path = expand("abundance/{num}_{{assembler}}_abund.tsv", num = list(range(len(SAMPLES))))
    output:
        "abundance/abund_{assembler}_maxbin_paths.txt"
    shell:
        """
        readlink -f {input.path} > {output}
        """


rule sample_abundance:
    input:
        "abundance/jgi_abundance_{assembler}.tsv", "abundance/temp_{assembler}.tsv"
    output:
        "abundance/{num}_{assembler}_abund.tsv"
    run:
        df = pd.read_csv(config['workdir'] + '/' + str(input[0]), sep = "\t")
        columns = df.columns.tolist()
        for i in columns:
            if i.endswith('var'):
                df = df.drop(i, axis =1)
            else:
                continue
        df.to_csv("abundance/jgi_abundance_nointra.tsv", index = False, sep = "\t")
        for i in range(len(SAMPLES)):
            shell("cut -f 1,"+ str((i+4)) + " abundance/jgi_abundance_nointra.tsv > abundance/" + str(i) + "_{wildcards.assembler}_abund.tsv")
        shell("for i in abundance/*{wildcards.assembler}_abund*; do sed -i '1d' $i;done")

rule general:
    input:
        "abundance/jgi_abundance_{assembler}.tsv"
    output:
        temp("abundance/temp_{assembler}.tsv")
    shell:
        "cut -f 1,2,3 {input} > {output}"
        
#Runs jgi and get abundance data
rule jgi:
    input:
        expand("data/mapped_reads/{sample}_{{assembler}}_sorted.bam", sample = set(SAMPLES)),
    output:
        "abundance/jgi_abundance_{assembler}.tsv",
    threads:
        int(1)
    params:
        walltime="86400"
    conda:
        "envs/metabat2.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input}
        """

#sort BAM file
rule samtools_sort:
    input:
        pair2 = "data/procesed_reads/{sample}_{assembler}.bam"
    output:
        "data/mapped_reads/{sample}_{assembler}_sorted.bam"
    params:
        walltime="86400"
    threads:
        int(config["sam_threads"])
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} -O bam -o {output} {input.pair2}
        """

#Creates BAM file from SAM
rule samtools_view:
    input:
        pairs = "data/procesed_reads/{sample}_{assembler}.sam"
    output:
        names = temp("data/procesed_reads/{sample}_{assembler}.bam")
    params: 
        walltime="86400"
    threads:
        int(config["sam_threads"])
    conda:
        "envs/samtools.yaml"
    shell: 
        """
        samtools view -@ {threads} -S -b -o {output.names} {input.pairs}
        """

#Maps reads to the contigs and creates a SAM file with the info    ation
rule bowtie:
    input:
        index = "{assembler}.bt_index",
        foward = contig("data/procesed_reads/{sample}_1_final.fastq"),
        rev = contig("data/procesed_reads/{sample}_2_final.fastq")
    output:
        pairs = temp("data/procesed_reads/{sample}_{assembler}.sam")    
    params:
        walltime="86400"
    threads:
        int(config["bwt_threads"])  
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        touch {output.pairs}
        bowtie2 -q -p {threads} -1 {input.foward} -2 {input.rev} -x index/{wildcards.assembler}.bt_index --no-unal -S {output.pairs}
        """
#Builds an index to the contig data

rule build_index:
    input:
        contig("data/contigs_{assembler}/contigs.fasta")
    output:
        temp("{assembler}.bt_index")
    threads:
        int(config["bwt_threads"])
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        touch {output}
        mkdir -p index
        bowtie2-build --threads {threads} {input} index/{output}
        """
