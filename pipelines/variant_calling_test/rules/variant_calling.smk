rule:
    input:
        mixture="{iteration}/simulated_mixture/distance_{distance}_coverage_{coverage}.simulated.fastq",
        reference = "references/poliovirus/Sabin_1_amplicon.fasta"
    output:
        