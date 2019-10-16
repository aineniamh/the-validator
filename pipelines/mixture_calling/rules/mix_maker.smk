rule mix_maker:
    input:
        sabin = "test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        other = "test/{iteration}/diverse_reads/coverage_200_distance_{distance}.simulated.fastq"
    output:
        "test/{iteration}/mixture_prop_varied/sabin_{coverage}_vs_distance_{distance}_200.simulated.fastq"
    shell:
        "cat {input.sabin} {input.other} > {output}"
