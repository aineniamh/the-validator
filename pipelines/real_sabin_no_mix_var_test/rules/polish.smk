rule minimap2_racon0:
    input:
        reads="real_reads/sabin_1.fastq",
        ref="test/{iteration}/diverse_refs/distance_{distance}.fasta"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/mapped.paf"
    shell:
        "minimap2 -x map-ont --secondary=no {input.ref} {input.reads} > {output}"

rule racon1:
    input:
        reads="real_reads/sabin_1.fastq",
        fasta="test/{iteration}/diverse_refs/distance_{distance}.fasta",
        paf= "test/{iteration}/polishing/distance_{distance}_real/mapped.paf"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/racon1.fasta"
    shell:
        "racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"

rule minimap2_racon1:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/polishing/distance_{distance}_real/racon1.fasta"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/mapped.racon1.paf"
    shell:
        "minimap2 -x map-ont --secondary=no {input.ref} {input.reads} > {output}"

rule racon2:
    input:
        reads="real_reads/sabin_1.fastq",
        fasta= "test/{iteration}/polishing/distance_{distance}_real/racon1.fasta",
        paf= "test/{iteration}/polishing/distance_{distance}_real/mapped.racon1.paf"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/racon2.fasta"
    shell:
        "racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"


rule minimap2_racon2:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/polishing/distance_{distance}_real/racon2.fasta"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/mapped.racon2.paf"
    shell:
        "minimap2 -x map-ont --secondary=no {input.ref} {input.reads} > {output}"

rule racon3:
    input:
        reads="real_reads/sabin_1.fastq",
        fasta= "test/{iteration}/polishing/distance_{distance}_real/racon2.fasta",
        paf= "test/{iteration}/polishing/distance_{distance}_real/mapped.racon2.paf"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/racon3.fasta"
    shell:
        "racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"


rule minimap2_racon3:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/polishing/distance_{distance}_real/racon3.fasta"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/mapped.racon3.paf"
    shell:
        "minimap2 -x map-ont --secondary=no {input.ref} {input.reads} > {output}"

rule racon4:
    input:
        reads="real_reads/sabin_1.fastq",
        fasta= "test/{iteration}/polishing/distance_{distance}_real/racon3.fasta",
        paf= "test/{iteration}/polishing/distance_{distance}_real/mapped.racon3.paf"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/racon4.fasta"
    shell:
        "racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output} || touch {output}"

rule minimap2_racon4:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/polishing/distance_{distance}_real/racon4.fasta"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/mapped.racon4.paf"
    shell:
        "minimap2 -x map-ont --secondary=no {input.ref} {input.reads} > {output}"

rule medaka:
    input:
        draft="test/{iteration}/polishing/distance_{distance}_real/racon4.fasta",
        basecalls="real_reads/sabin_1.fastq"
        # draft= "test/{iteration}/diverse_refs/distance_{distance}.fasta"
    params:
        outdir="test/{iteration}/consensus/distance_{distance}_real"
    output:
        consensus = "test/{iteration}/consensus/distance_{distance}_real/consensus.fasta"
    threads:
        2
    shell:
        "medaka_consensus -i {input.basecalls} -d {input.draft} -o {params.outdir} -t 2  || touch {output}"


rule minimap2_medaka_paf:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/consensus/distance_{distance}_real/consensus.fasta"
    output:
        "test/{iteration}/polishing/distance_{distance}_real/mapped.medaka.paf"
    shell:
        "minimap2 -x map-ont --secondary=no {input.ref} {input.reads} > {output}"


rule medaka_from_scratch:
    input:
        draft="test/{iteration}/diverse_refs/distance_{distance}.fasta",
        basecalls="real_reads/sabin_1.fastq"
    params:
        outdir="test/{iteration}/consensus/distance_{distance}_real_scratch"
    output:
        consensus = "test/{iteration}/consensus/distance_{distance}_real_scratch/consensus.fasta"
    threads:
        2
    shell:
        "medaka_consensus -i {input.basecalls} -d {input.draft} -o {params.outdir} -t 2  || touch {output}"


rule minimap2_medaka_paf_scratch:
    input:
        reads="real_reads/sabin_1.fastq",
        ref= "test/{iteration}/consensus/distance_{distance}_real_scratch/consensus.fasta"
    output:
        "test/{iteration}/polishing/distance_{distance}_real_scratch/mapped.medaka.paf"
    shell:
        "minimap2 -x map-ont --secondary=no {input.ref} {input.reads} > {output}"