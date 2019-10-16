rule minimap2_racon0:
    input:
        reads="real_reads/sabin_1.fastq",
        ref="test/{iteration}/diverse_refs/distance_{distance}.fasta"
    output:
        temp("test/{iteration}/polishing/distance_{distance}_real/mapped.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon1:
    input:
        reads="real_reads/sabin_1.fastq",
        fasta="test/{iteration}/diverse_refs/distance_{distance}.fasta",
        paf= "test/{iteration}/polishing/distance_{distance}_real/mapped.paf"
    output:
        temp("test/{iteration}/polishing/distance_{distance}_real/racon1.fasta")
    shell:
        "~/Documents/test_racon/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"

rule minimap2_racon1:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/polishing/distance_{distance}_real/racon1.fasta"
    output:
        temp("test/{iteration}/polishing/distance_{distance}_real/mapped.racon1.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon2:
    input:
        reads="real_reads/sabin_1.fastq",
        fasta= "test/{iteration}/polishing/distance_{distance}_real/racon1.fasta",
        paf= "test/{iteration}/polishing/distance_{distance}_real/mapped.racon1.paf"
    output:
        temp("test/{iteration}/polishing/distance_{distance}_real/racon2.fasta")
    shell:
        "~/Documents/test_racon/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"


rule minimap2_racon2:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/polishing/distance_{distance}_real/racon2.fasta"
    output:
        temp("test/{iteration}/polishing/distance_{distance}_real/mapped.racon2.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon3:
    input:
        reads="real_reads/sabin_1.fastq",
        fasta= "test/{iteration}/polishing/distance_{distance}_real/racon2.fasta",
        paf= "test/{iteration}/polishing/distance_{distance}_real/mapped.racon2.paf"
    output:
        temp("test/{iteration}/polishing/distance_{distance}_real/racon3.fasta")
    shell:
        "~/Documents/test_racon/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"


rule minimap2_racon3:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/polishing/distance_{distance}_real/racon3.fasta"
    output:
        temp("test/{iteration}/polishing/distance_{distance}_real/mapped.racon3.paf")
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon4:
    input:
        reads="real_reads/sabin_1.fastq",
        fasta= "test/{iteration}/polishing/distance_{distance}_real/racon3.fasta",
        paf= "test/{iteration}/polishing/distance_{distance}_real/mapped.racon3.paf"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/racon4.fasta"
    shell:
        "~/Documents/test_racon/racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"

rule minimap2_racon4:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/polishing/distance_{distance}_real/racon4.fasta"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/mapped.racon4.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule medaka:
    input:
        draft="test/{iteration}/polishing/distance_{distance}_real/racon4.fasta",
        basecalls="real_reads/sabin_1.fastq"
        # draft= "test/{iteration}/diverse_refs/distance_{distance}.fasta"
    params:
        outdir="test/{iteration}/consensus/distance_{distance}_coverage_200"
    output:
        consensus = "test/{iteration}/consensus/distance_{distance}_coverage_200/consensus.fasta"
    threads:
        2
    shell:
        "medaka_consensus -i {input.basecalls} -d {input.draft} -o {params.outdir} -t 2  || touch {output}"

