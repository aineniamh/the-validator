
    
rule minimap2_racon0:
    input:
        reads="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        ref="{iteration}/diverse_refs/distance_{distance}.fasta"
    output:
        temp("{iteration}/polishing_medaka/distance_{distance}_{coverage}X/mapped.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

# rule racon1:
#     input:
#         reads="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
#         fasta="{iteration}/diverse_refs/distance_{distance}.fasta",
#         paf= "{iteration}/polishing/distance_{distance}_{coverage}X/mapped.paf"
#     output:
#         temp("{iteration}/polishing/distance_{distance}_{coverage}X/racon1.fasta")
#     shell:
#         "~/Documents/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"

# rule minimap2_racon1:
#     input:
#         reads="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
#         ref= "{iteration}/polishing/distance_{distance}_{coverage}X/racon1.fasta"
#     output:
#         temp("{iteration}/polishing/distance_{distance}_{coverage}X/mapped.racon1.paf")
#     shell:
#         "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

# rule racon2:
#     input:
#         reads="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
#         fasta= "{iteration}/polishing/distance_{distance}_{coverage}X/racon1.fasta",
#         paf= "{iteration}/polishing/distance_{distance}_{coverage}X/mapped.racon1.paf"
#     output:
#         temp("{iteration}/polishing/distance_{distance}_{coverage}X/racon2.fasta")
#     shell:
#         "~/Documents/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"


# rule minimap2_racon2:
#     input:
#         reads="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
#         ref= "{iteration}/polishing/distance_{distance}_{coverage}X/racon2.fasta"
#     output:
#         temp("{iteration}/polishing/distance_{distance}_{coverage}X/mapped.racon2.paf")
#     shell:
#         "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

# rule racon3:
#     input:
#         reads="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
#         fasta= "{iteration}/polishing/distance_{distance}_{coverage}X/racon2.fasta",
#         paf= "{iteration}/polishing/distance_{distance}_{coverage}X/mapped.racon2.paf"
#     output:
#         temp("{iteration}/polishing/distance_{distance}_{coverage}X/racon3.fasta")
#     shell:
#         "~/Documents/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"


# rule minimap2_racon3:
#     input:
#         reads="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
#         ref= "{iteration}/polishing/distance_{distance}_{coverage}X/racon3.fasta"
#     output:
#         temp("{iteration}/polishing/distance_{distance}_{coverage}X/mapped.racon3.paf")
#     shell:
#         "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

# rule racon4:
#     input:
#         reads="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
#         fasta= "{iteration}/polishing/distance_{distance}_{coverage}X/racon3.fasta",
#         paf= "{iteration}/polishing/distance_{distance}_{coverage}X/mapped.racon3.paf"
#     output:
#         temp("{iteration}/polishing/distance_{distance}_{coverage}X/racon4.fasta")
#     shell:
#         "~/Documents/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"

# rule minimap2_racon4:
#     input:
#         reads="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
#         ref= "{iteration}/polishing/distance_{distance}_{coverage}X/racon4.fasta"
#     output:
#         temp("{iteration}/polishing/distance_{distance}_{coverage}X/mapped.racon4.paf")
#     shell:
#         "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule medaka:
    input:
        # draft="{iteration}/polishing/distance_{distance}_{coverage}X/racon4.fasta",
        basecalls="{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        draft= "{iteration}/diverse_refs/distance_{distance}.fasta"
    params:
        outdir="{iteration}/consensus_medaka_sequences/distance_{distance}_coverage_{coverage}"
    output:
        consensus = "{iteration}/consensus_medaka_sequences/distance_{distance}_coverage_{coverage}/consensus.fasta"
    threads:
        2
    shell:
        "medaka_consensus -i {input.basecalls} -d {input.draft} -o {params.outdir} -t 2  || touch {output}"

