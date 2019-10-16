rule map_to_ref:
    input:
        reads="test/{iteration}/mix_reads_distance_varied/sabin_and_distance_{distance}_coverage_200.simulated.fastq",
        reference = "test/{iteration}/diverse_refs/distance_{distance}.fasta"
    output:
        "test/{iteration}/mix/map_ref/distance_{distance}_coverage_200X.bam"
    shell:
        "minimap2 -ax map-ont {input.reference} {input.reads} | samtools sort -O bam -o {output}"

rule index_map_to_ref:
    input:
        "test/{iteration}/mix/map_ref/distance_{distance}_coverage_200X.bam"
    output:
        "test/{iteration}/mix/map_ref/distance_{distance}_coverage_200X.bam.bai"
    shell:
        "samtools index {input}"

rule map_to_cns:
    input:
        reads="test/{iteration}/mix_reads_distance_varied/sabin_and_distance_{distance}_coverage_200.simulated.fastq",
        reference = "test/{iteration}/mix/consensus/distance_{distance}_coverage_200/consensus.fasta"
    output:
        "test/{iteration}/mix/map_cns/distance_{distance}_coverage_200X.bam"
    shell:
        "minimap2 -ax map-ont {input.reference} {input.reads} | samtools sort -O bam -o {output}"

rule index_map_to_cns:
    input:
        "test/{iteration}/mix/map_cns/distance_{distance}_coverage_200X.bam"
    output:
        "test/{iteration}/mix/map_cns/distance_{distance}_coverage_200X.bam.bai"
    shell:
        "samtools index {input}"

rule polish_call:
    input:
        consensus = "test/{iteration}/mix/consensus/distance_{distance}_coverage_200/consensus.fasta",
        map = "test/{iteration}/mix/map_cns/distance_{distance}_coverage_200X.bam",
        index = "test/{iteration}/mix/map_cns/distance_{distance}_coverage_200X.bam.bai"
    params:
        dir = "test/{iteration}/mix/calls_cns/distance_{distance}_coverage_200"
    output:
        "test/{iteration}/mix/calls_cns/distance_{distance}_coverage_200/round_2_final_unphased.vcf"
    shell:
        "medaka_variant -f {input.consensus} -i {input.map} -o {params.dir} -m r941_trans || touch {output}"

rule call:
    input:
        ref = "test/{iteration}/diverse_refs/distance_{distance}.fasta",
        map = "test/{iteration}/mix/map_ref/distance_{distance}_coverage_200X.bam",
        index = "test/{iteration}/mix/map_ref/distance_{distance}_coverage_200X.bam.bai"
    params:
        dir = "test/{iteration}/mix/calls_ref/distance_{distance}_coverage_200"
    output:
        "test/{iteration}/mix/calls_ref/distance_{distance}_coverage_200/round_2_final_unphased.vcf"
    shell:
        "medaka_variant -f {input.ref} -i {input.map} -o {params.dir} -m r941_trans || touch {output}"

rule make_fasta_sabin_vs_ref_vs_cns:
    input:
        ref = "test/{iteration}/diverse_refs/distance_{distance}.fasta",
        sabin = "references/poliovirus/Sabin_1_amplicon.fasta",
        cns = "test/{iteration}/mix/consensus/distance_{distance}_coverage_200/consensus.fasta"
    output:
        three = "test/{iteration}/mix/alignments/distance_{distance}_coverage_200X.fasta",
        ref_cns = "test/{iteration}/mix/alignments/distance_{distance}_ref_cns.fasta",
        sabin_cns = "test/{iteration}/mix/alignments/distance_{distance}_sabin_cns.fasta"
    run:
        fw=open(str(output.three),"w")
        fw2=open(str(output.ref_cns),"w")
        fw3=open(str(output.sabin_cns),"w")
        for record in SeqIO.parse(str(input.ref),"fasta"):
            fw.write(">{}\n{}\n".format("diverged_ref",record.seq))
            fw2.write(">{}\n{}\n".format("diverged_ref",record.seq))
        for record in SeqIO.parse(str(input.sabin),"fasta"):
            fw.write(">{}\n{}\n".format("sabin",record.seq))
            fw3.write(">{}\n{}\n".format("sabin",record.seq))
        for record in SeqIO.parse(str(input.cns),"fasta"):
            fw.write(">{}\n{}\n".format("consensus",record.seq))
            fw2.write(">{}\n{}\n".format("consensus",record.seq))
            fw3.write(">{}\n{}\n".format("consensus",record.seq))
        fw.close()
        fw2.close()
        fw3.close()

rule align_sabin_ref_cns:
    input:
        "test/{iteration}/mix/alignments/distance_{distance}_coverage_200X.fasta"
    output:
        "test/{iteration}/mix/alignments/distance_{distance}_coverage_200X.aln.fasta"
    shell:
        "mafft {input} > {output}"


rule align_ref_cns:
    input:
        "test/{iteration}/mix/alignments/distance_{distance}_ref_cns.fasta"
    output:
        "test/{iteration}/mix/alignments/distance_{distance}_ref_cns.aln.fasta"
    shell:
        "mafft {input} > {output}"

rule align_sabin_cns:
    input:
        "test/{iteration}/mix/alignments/distance_{distance}_sabin_cns.fasta"
    output:
        "test/{iteration}/mix/alignments/distance_{distance}_sabin_cns.aln.fasta"
    shell:
        "mafft {input} > {output}"