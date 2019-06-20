import os
import collections
from Bio import SeqIO

configfile: "config.yaml"

##### Target rules #####

rule all:
    input:
        expand("pipeline_output/consensus_genomes/{barcode}.fasta",barcode=config["barcodes"])
        

##### Modules #####
include: "rules/gather.smk"
include: "rules/nanopolish_index.smk"
include: "rules/demultiplex.smk"
include: "rules/bin.smk"
include: "rules/mapping.smk"
include: "rules/sorting_calling_generate_cns.smk"
include: "rules/minion.smk"
include: "rules/generate_consensus.smk"

onstart:
    print("Setting up the artic package")
    shell("cd fieldbioinformatics && python setup.py install")
    shell("export PATH=$PATH:`pwd`/artic")
    shell("cd .. ")
