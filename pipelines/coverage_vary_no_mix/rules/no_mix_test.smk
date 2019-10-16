rule map_to_ref:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        reference = "test/{iteration}/diverse_refs/distance_0.01.fasta"
    output:
        "test/{iteration}/no_mix/map_ref/distance_0.01_coverage_{coverage}X.bam"
    shell:
        "minimap2 -ax map-ont {input.reference} {input.reads} | samtools sort -O bam -o {output}"

rule index_map_to_ref:
    input:
        "test/{iteration}/no_mix/map_ref/distance_0.01_coverage_{coverage}X.bam"
    output:
        "test/{iteration}/no_mix/map_ref/distance_0.01_coverage_{coverage}X.bam.bai"
    shell:
        "samtools index {input}"

rule map_to_cns:
    input:
        reads="test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq",
        reference = "test/{iteration}/consensus/distance_0.01_coverage_{coverage}/consensus.fasta"
    output:
        "test/{iteration}/no_mix/map_cns/distance_0.01_coverage_{coverage}X.bam"
    shell:
        "minimap2 -ax map-ont {input.reference} {input.reads} | samtools sort -O bam -o {output}"

rule index_map_to_cns:
    input:
        "test/{iteration}/no_mix/map_cns/distance_0.01_coverage_{coverage}X.bam"
    output:
        "test/{iteration}/no_mix/map_cns/distance_0.01_coverage_{coverage}X.bam.bai"
    shell:
        "samtools index {input}"

rule polish_call:
    input:
        consensus = "test/{iteration}/consensus/distance_0.01_coverage_{coverage}/consensus.fasta",
        map = "test/{iteration}/no_mix/map_cns/distance_0.01_coverage_{coverage}X.bam",
        index = "test/{iteration}/no_mix/map_cns/distance_0.01_coverage_{coverage}X.bam.bai"
    params:
        dir = "test/{iteration}/no_mix/calls_cns/distance_0.01_coverage_{coverage}"
    output:
        "test/{iteration}/no_mix/calls_cns/distance_0.01_coverage_{coverage}/round_2_final_phased.vcf"
    shell:
        "medaka_variant -f {input.consensus} -i {input.map} -o {params.dir} -m r941_trans || touch {output}"

rule call:
    input:
        ref = "test/{iteration}/diverse_refs/distance_0.01.fasta",
        map = "test/{iteration}/no_mix/map_ref/distance_0.01_coverage_{coverage}X.bam",
        index = "test/{iteration}/no_mix/map_ref/distance_0.01_coverage_{coverage}X.bam.bai"
    params:
        dir = "test/{iteration}/no_mix/calls_ref/distance_0.01_coverage_{coverage}"
    output:
        "test/{iteration}/no_mix/calls_ref/distance_0.01_coverage_{coverage}/round_2_final_phased.vcf"
    shell:
        "medaka_variant -f {input.ref} -i {input.map} -o {params.dir} -m r941_trans || touch {output}"

rule make_fasta_sabin_vs_ref_vs_cns:
    input:
        ref = "test/{iteration}/diverse_refs/distance_0.01.fasta",
        sabin = "references/poliovirus/Sabin_1_amplicon.fasta",
        cns = "test/{iteration}/consensus/distance_0.01_coverage_{coverage}/consensus.fasta"
    output:
        three = "test/{iteration}/alignments/distance_0.01_coverage_{coverage}X.fasta",
        ref_cns = "test/{iteration}/alignments/distance_0.01_coverage_{coverage}_ref_cns.fasta",
        sabin_cns = "test/{iteration}/alignments/distance_0.01_coverage_{coverage}_sabin_cns.fasta"
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
        "test/{iteration}/alignments/distance_0.01_coverage_{coverage}X.fasta"
    output:
        "test/{iteration}/alignments/distance_0.01_coverage_{coverage}X.aln.fasta"
    shell:
        "mafft {input} > {output}"


rule align_ref_cns:
    input:
        "test/{iteration}/alignments/distance_0.01_coverage_{coverage}_ref_cns.fasta"
    output:
        "test/{iteration}/alignments/distance_0.01_coverage_{coverage}_ref_cns.aln.fasta"
    shell:
        "mafft {input} > {output}"

rule align_sabin_cns:
    input:
        "test/{iteration}/alignments/distance_0.01_coverage_{coverage}_sabin_cns.fasta"
    output:
        "test/{iteration}/alignments/distance_0.01_coverage_{coverage}_sabin_cns.aln.fasta"
    shell:
        "mafft {input} > {output}"

