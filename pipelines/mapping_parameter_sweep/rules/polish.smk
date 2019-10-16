rule minimap2:
    input:
        reads="{iteration}/simulated_reads/coverage_200.simulated.fastq",
        ref="{iteration}/diverse_refs/distance_{distance}.fasta"
    params:
        
    output:
        temp("{iteration}/polishing/distance_{distance}_{coverage}X/mapped.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"
