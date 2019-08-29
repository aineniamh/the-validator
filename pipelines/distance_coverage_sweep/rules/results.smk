
rule get_cns:
    input:
    params:
        input_dir = "{iteration}"
    output:
        "{iteration}/consensus_sequences.fasta"
    run:
        c = 0
        with open(str(output),"w") as fw:
            for r, d, f in os.walk(params.input_dir):
                for filename in f:
                    if filename=="consensus.fasta":
                        d = r.split("/")[-1]
                        try:
                            for record in SeqIO.parse(r + "/" + filename,"fasta"):
                                fw.write(">{}\n{}\n".format(d, record.seq))

                                c+=1
                        except:
                            continue
        print(c)

rule make_blast_db:
    input:
        "references/poliovirus/Sabin_1_amplicon.fasta"
    output:
        "references/poliovirus/Sabin_1_amplicon.fasta.nhr"
    shell:
        "makeblastdb -in {input} -dbtype nucl -out {input}"

rule blast:
    input:
        db= "references/poliovirus/Sabin_1_amplicon.fasta",
        db_hidden = "references/poliovirus/Sabin_1_amplicon.fasta.nhr",
        seqs = "{iteration}/consensus_sequences.fasta"
    output:
        "{iteration}/blast_against_sabin.csv"
    shell:
        "blastn -db {input.db} -query {input.seqs} -outfmt 10 -out {output}"

rule get_accuracy:
    input:
        "{iteration}/blast_against_sabin.csv"
    params:
        coverage = config["coverages"],
        distance = config["distances"]
    output:
        "{iteration}/distance_coverage_sweep.csv"
    run:
        combo_dict = {}
        for i in params.distance:
            for j in params.coverage:
                combo_dict["distance_{}_{}X".format(i,j)]=(i,j,0,0,0)
        with open(str(input),"r") as f:
            for l in f:
                l = l.rstrip('\n')
                tokens=l.split(',')
                d = tokens[0]
                distance = d.split('_')[1]
                coverage = d.split('_')[3]
                identity=tokens[2]
                gaps = tokens[5]
                mismatch = tokens[4]
                key = "distance_{}_{}X".format(distance, coverage)
                combo_dict[key]=(distance, coverage, identity, gaps, mismatch)
        fw = open(str(output),"w")
        fw.write("distance,coverage,identity,gaps,mismatch\n")
        for key in combo_dict:
            print(key, combo_dict[key])
            fw.write("{},{},{},{},{}\n".format(combo_dict[key][0],combo_dict[key][1],combo_dict[key][2],combo_dict[key][3],combo_dict[key][4]))
        fw.close()

rule make_figs:
    input:
        expand("{iteration}/distance_coverage_sweep.csv", iteration=config["iterations"])
    output:
        # distance ="distance_vs_accuracy.png",
        # coverage = "coverage_vs_accuracy.png",
        both = "coverage_distance_accuracy.png"
    run:

          # noqa: F401 unused import
        fig = plt.figure(figsize=(15,10))
        plt.rcParams.update({'font.size': 16})

        ax = fig.add_subplot(111, projection='3d')
        plt.xlabel("Distance")
        plt.ylabel("Coverage")
        ax.set_zlabel('%Identity')
        print(str(input))
        for filename in str(input).split(' '):
            data = pd.read_csv(filename)
            ax.scatter(data["distance"],data["coverage"],data["identity"],color="#294763" )
        fig.savefig(str(output.both))


