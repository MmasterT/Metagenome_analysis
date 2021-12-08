# Metagenome_analysis

SnakeMake Pipeline to retrive MAGs from Illumina pair-end metagenomic sequencing data. Download the repositorie and copy your raw reads in the raw_reads directory. 

## Pipeline

This pipeline works with snakemake version 5.8.2. The main snakefile is metaPipe1.1.1.py. You run the main snakefile and wait for the results. The number after the -j argument is the amount of threads you want to assign to the whole run. 

```
snakemake --use-conda -j 8 -s metaPipe1.1.1.py
```

For more information of SnakeMake use the oficial documentation: https://snakemake.readthedocs.io/en/v5.8.2/


## Configuration

The configuration of parameters is via the config.json file in the folder with thees parameters. It is ***very important*** to not change the left side of the : but the rigth side. Always the changes are between the quotes. The first thing it is to change the workdir to the directory path of the file where.

```{p Carga de Datos, echo = True}
{
        "workdir": "/u02/Mofedi/MAG_retrive_60genomas/", #path of the directory with all the files
        "t5" : "10", #amount of bases to cut from 5'
        "t3" : "5", #amount of bases to cut from 3'
        "Lmin" : "50", #minimum lenth of the reads
        "Lmax" : "600", #maximum length of the reads
        "Q" : "28", #minimum average quality of the reads 
        "N" : "False", #do not print amount of N bases
        "min_read_length" : "50",
        "spades_threads" : "42", #threads to use in metaSPAdes
        "spades_mem" : "350", #memory assigned to metaSPAdes in Gb
        "k-list" : "35,47,69,91", #list of kmeres to use
        "megahit_threads" : "42",
        "mega_mem" : "100",
        "min_contig_length" : "2000", #minimum length of the ontigs to use in binning
        "zip" : "--nozip", #concatenate contigs to use in vamb with no compression
        "bwt_threads" : "21",
        "sam_threads" : "21",
        "maxb2_threads" : "42",
        "meta_threads" : "42",
        "vamb_threads" : "42",
        "min_length_bin" : "500000", #filter contigs with this length or less
        "drep_thread" : "42",
        "max_contamination" : "15", #Filter bins with more than this contaminatiom
        "min_completness" : "75", #Filter bins with less than this completness
        "sa" : "0.98", #ANI similarity threshold for dRep
        "contigs" : "contigs.txt", #do not change
        "sample_data" : "samples2data.txt", #do not change
        "index_size" : "3G",# do not change
        "minimap_mem" : "80gb",
        "minimap_ppn" : "21",
        "vamb_mem" : "80gb",
        "vamb_ppn" : "42",
        "vamb_params" : "-o C -m 2000 --minfasta 500000",
        "vamv_preload" : "",
        "assembler_method" : "concated" #assembler suports concated reads method and single sample method. 
}
```
For the pipeline to work you have te change in all the .py files the path from the configfile variable to the path to your config.json file.
