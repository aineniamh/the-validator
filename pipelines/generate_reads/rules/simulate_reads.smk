
# rule split_to_amplicons:
#     input:
#         fasta = config["outputDir"]+"/{reference}/diversity/distance_{distance}.fasta",
#         amplicons = config["ampliconFilePath"]
#     params:
#         amps=config["amplicons"]
#     output:
#         expand(config["outputDir"]+"/{{reference}}/amplicons/distance_{{distance}}_{amplicon}.fasta", amplicon=config["amplicons"])
#     run:
#         ref,div,distance= str(input.fasta).rstrip('.fasta').lstrip('/').split('/')[-3:]
#         for record in SeqIO.parse(str(input.fasta), "fasta"):
#             for key in params.amps:
#                 header = record.description
#                 header += " amplicon={} coords={}".format(key,params.amps[key])
#                 outfile = config["outputDir"]+"/{}/amplicons/{}_{}.fasta".format(ref,distance,key)
#                 fw = open(outfile,"w")
#                 start,end = params.amps[key]
#                 amp_seq = record.seq[start:end]
#                 fw.write(">{}\n{}\n".format(header, amp_seq))
#                 fw.close()

rule simulate_reads:
    input:
        "references/poliovirus/Sabin_1_amplicon.fasta"
    params:
        coverage="{coverage}"
    output:
        "test/{iteration}/simulated_reads/coverage_{coverage}.simulated.fastq"
    shell:
        "badread simulate --reference {input} --quantity {params.coverage}X --length 1107,4 --identity 80,95,5 > {output}"


# rule make_sample_file:
#     input:
#         expand(config["outputDir"]+"/{{reference}}/simulated_amplicons/distance_{{distance}}_{amplicon}.simulated.fastq", amplicon=config["amplicons"])
#     output:
#         config["outputDir"]+"/{reference}/simulated_sample/distance_{distance}.simulated.fastq"
#     shell:
#         "cat {input} > {output}"


