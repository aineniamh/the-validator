rule minimap2_racon0:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        ref="test/{iteration}/diverse_refs/distance_0.01.fasta"
    output:
        temp("test/{iteration}/polishing/distance_0.01_{coverage}X/mapped.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon1:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        fasta="test/{iteration}/diverse_refs/distance_0.01.fasta",
        paf= "test/{iteration}/polishing/distance_0.01_{coverage}X/mapped.paf"
    output:
        "test/{iteration}/polishing/distance_0.01_{coverage}X/racon1.fasta"
    shell:
        "~/Documents/test_racon/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"

rule minimap2_racon1:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        ref= "test/{iteration}/polishing/distance_0.01_{coverage}X/racon1.fasta"
    output:
        temp("test/{iteration}/polishing/distance_0.01_{coverage}X/mapped.racon1.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon2:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        fasta= "test/{iteration}/polishing/distance_0.01_{coverage}X/racon1.fasta",
        paf= "test/{iteration}/polishing/distance_0.01_{coverage}X/mapped.racon1.paf"
    output:
        "test/{iteration}/polishing/distance_0.01_{coverage}X/racon2.fasta"
    shell:
        "~/Documents/test_racon/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"


rule minimap2_racon2:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        ref= "test/{iteration}/polishing/distance_0.01_{coverage}X/racon2.fasta"
    output:
        temp("test/{iteration}/polishing/distance_0.01_{coverage}X/mapped.racon2.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon3:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        fasta= "test/{iteration}/polishing/distance_0.01_{coverage}X/racon2.fasta",
        paf= "test/{iteration}/polishing/distance_0.01_{coverage}X/mapped.racon2.paf"
    output:
        "test/{iteration}/polishing/distance_0.01_{coverage}X/racon3.fasta"
    shell:
        "~/Documents/test_racon/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"


rule minimap2_racon3:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        ref= "test/{iteration}/polishing/distance_0.01_{coverage}X/racon3.fasta"
    output:
        temp("test/{iteration}/polishing/distance_0.01_{coverage}X/mapped.racon3.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon4:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        fasta= "test/{iteration}/polishing/distance_0.01_{coverage}X/racon3.fasta",
        paf= "test/{iteration}/polishing/distance_0.01_{coverage}X/mapped.racon3.paf"
    output:
        "test/{iteration}/polishing/distance_0.01_{coverage}X/racon4.fasta"
    shell:
        "~/Documents/test_racon/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"

rule minimap2_racon4:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        ref= "test/{iteration}/polishing/distance_0.01_{coverage}X/racon4.fasta"
    output:
        "test/{iteration}/polishing/distance_0.01_{coverage}X/mapped.racon4.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule medaka:
    input:
        draft="test/{iteration}/polishing/distance_0.01_{coverage}X/racon4.fasta",
        basecalls="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq"
        # draft= "test/{iteration}/diverse_refs/distance_0.01.fasta"
    params:
        outdir="test/{iteration}/consensus/distance_0.01_coverage_{coverage}"
    output:
        consensus = "test/{iteration}/consensus/distance_0.01_coverage_{coverage}/consensus.fasta"
    threads:
        2
    shell:
        "medaka_consensus -i {input.basecalls} -d {input.draft} -o {params.outdir} -t 2  || touch {output}"

