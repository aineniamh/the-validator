rule create_diversity:
    input:
        tree="references/norovirus/tree.nw"
    output:
        "diverse_references/reverse_genetics_gii6.fasta"
    shell:
        "seq-gen -l 7573 -s 1 -m HKY -z 12345 -t 0.2 -of -k 1 {input.tree} > {output}"
 

rule makeblastdb:
    input:
        "primer-schemes/noro2kb/V2/noro2kb.amplicons.fasta"
    output:
        "primer-schemes/noro2kb/V2/noro2kb.amplicons.fasta.nhr"
    shell:
        "makeblastdb -in {input} -dbtype nucl"

rule fastq_to_fasta:
    input:
        "pipeline_output/demultiplexed/{barcode}.fastq"
    output:
        "pipeline_output/demultiplexed/{barcode}.fasta"
    shell:
        "seqtk seq -A {input} > {output}"

rule blastn:
    input:
        db_hidden="primer-schemes/noro2kb/V2/noro2kb.amplicons.fasta.nhr",
        db="primer-schemes/noro2kb/V2/noro2kb.amplicons.fasta",
        reads="pipeline_output/demultiplexed/{barcode}.fasta" 
    output:
        "pipeline_output/blast_results/{barcode}.blast.csv"
    shell:
        "blastn -task blastn -db {input.db} "
        "-query {input.reads} -out {output} "
        "-num_threads 1 -outfmt 10"
            
rule bin:
    input:
        blast="pipeline_output/blast_results/{barcode}.blast.csv",
        reads="pipeline_output/demultiplexed/{barcode}.fastq",
        refs="primer-schemes/noro2kb/V2/noro2kb.amplicons.fasta",
        amps="primer-schemes/noro2kb/V2/noro2kb.amplicons.csv"
    params:
        outdir="pipeline_output/binned/{barcode}_bin/reads/",
        sample="{barcode}"
    output:
        summary="pipeline_output/binned/{barcode}_bin/binning_report.txt",
        ref=expand("pipeline_output/binned/{{barcode}}_bin/primer-schemes/minion/V_{amplicon}/minion.reference.fasta", amplicon=config["amplicons"]),
        bed=expand("pipeline_output/binned/{{barcode}}_bin/primer-schemes/minion/V_{amplicon}/minion.scheme.bed", amplicon=config["amplicons"]),
        reads=expand("pipeline_output/binned/{{barcode}}_bin/reads/{amplicon}.fastq", amplicon=config["amplicons"])
    shell:
        "python scripts/bin.py --blast_file {input.blast} "
        "--reference_file {input.refs} --reads {input.reads} "
        "--output_dir {params.outdir} --summary {output.summary} --bed_file {input.amps} "
        "--sample {params.sample}"
