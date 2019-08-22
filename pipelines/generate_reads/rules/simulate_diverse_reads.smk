rule simulate_d_reads:
    input:
        "{iteration}/diverse_refs/distance_{distance}.fasta"
    output:
        "{iteration}/simulated_diverse_reads/distance_{distance}.simulated.fastq"
    shell:
        "badread simulate --reference {input} --quantity 100X --length 1107,4 --identity 80,95,5 > {output}"

rule make_mixture:
    input:
        diverse_reads = "{iteration}/simulated_diverse_reads/distance_{distance}.simulated.fastq",
        coverage_varied = "{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq"
    output:
        "{iteration}/simulated_mixture/distance_{distance}_coverage_{coverage}.simulated.fastq"
    shell:
        "cat {input.diverse_reads} {input.coverage_varied} > {output}"