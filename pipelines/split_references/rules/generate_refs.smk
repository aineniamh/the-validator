rule split_reference_panel:
    input:
        ref= config["referencePanelPath"],
        config = config["configPath"]
    params:
        outputDir = config["outputDir"]
    output:
        config["outputDir"] + "/validator_config.yaml"
    run:
        if not os.path.exists(config["outputDir"]):
            os.mkdir(config["outputDir"])
        
        records = []
        for record in SeqIO.parse(str(input.ref), "fasta"):
            if not os.path.exists(config["outputDir"]+'/'+record.id):
                os.mkdir(config["outputDir"]+'/'+record.id)
            fw = open(str(params.outputDir) +"/"+record.id + '/diversity/distance_0.fasta',"w")
            fw.write(">{}\n{}\n".format(record.id, record.seq))
            fw.close()
            records.append(record.id)
        with open(str(output),"w") as fw:
            with open(str(input.config),"r") as f:
                for l in f:
                    l = l.rstrip('\n')
                    fw.write(l + '\n')
            fw.write("references:\n")
            for i in records:
                fw.write("  - {}\n".format(i))

