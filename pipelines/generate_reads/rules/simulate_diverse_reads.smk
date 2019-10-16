rule simulate_diverse_reads:
    input:
        "test/{iteration}/diverse_refs/distance_{distance}.fasta"
    params:
        coverage=config["coverage"]
    output:
        "test/{iteration}/diverse_reads/coverage_200_distance_{distance}.simulated.fastq"
    shell:
        "badread simulate --reference {input} --quantity {params.coverage}X --length 1107,4 --identity 80,95,5 > {output}"

rule make_mixture:
    input:
        diverse_reads = "test/{iteration}/diverse_reads/coverage_200_distance_{distance}.simulated.fastq",
        coverage_varied = "test/{iteration}/simulated_reads/coverage_200.simulated.fastq"
    output:
        "test/{iteration}/mix_reads_distance_varied/sabin_and_distance_{distance}_coverage_200.simulated.fastq"
    shell:
        "cat {input.diverse_reads} {input.coverage_varied} > {output}"