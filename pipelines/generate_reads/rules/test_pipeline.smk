rule organise_simulations:
    input:
        config["outputDir"]+"/{reference}/simulated_sample/distance_{distance}.simulated.fastq"
    params:
        ref ="{reference}",
        distance = "{distance}"
    output:
        

rule run_pipeline:
    input:
        simulated_reads=config["outputDir"]+"/{reference}/simulated_sample/distance_{distance}.simulated.fastq",
        config = "../artic-noro/pipelines/master_consensus/config.yaml"
    params:
        barcode="{distance}"
        pipeline_output= config["outputDir"]
    output:
        config["outputDir"] + "/{reference}/consensus_sequences/{distance}.fasta"
    shell:
        "snakemake --nolock --snakefile ../artic-noro/pipelines/polish_consensus/Snakefile "
        "--configfile {input.config} "
        "--config outputPath={params.pipeline_output} "
        "barcode={params.barcode}"