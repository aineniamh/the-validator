import os
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='organise the simulated files.')
parser.add_argument("--dir", action="store", type=str, dest="dir")
# parser.add_argument("--distance", action="store", type=float, dest="distance")
parser.add_argument("--config", action="store", type=str, dest="config")
# parser.add_argument("--outfile", action="store", type=str, dest="outfile")

args = parser.parse_args()
barcodes = []
for r,d,f in os.walk(str(args.dir)):
    path = r
    if path.endswith("simulated_sample"):
        dirs = path.split('/')
        ref = dirs[-2]
        for filename in f:
            if filename.endswith("fastq"):
                distance = filename.rstrip('.simulated.fastq')
                barcode = "{}-{}".format(ref,distance)
                barcodes.append(barcode)
                # records = []
                # for record in SeqIO.parse(r+ '/' + filename,"fastq"):
                #     header = str(record.description)
                #     header += " barcode={} reference_hit={} ".format(barcode, ref)
                #     record.description = header
                #     records.append(record)
                # new_file ="/Users/aineniamh/Documents/artic-noro/simulated_test/annotated_reads/barcode_{}.fastq".format(barcode)
                # with open (new_file,"w") as fw:
                #     SeqIO.write(records, fw,"fastq")

with open("/Users/aineniamh/Documents/artic-noro/simulated_test/config.yaml","w") as fw:
    fw.write("barcodes:\n")
    for i in barcodes:
        fw.write("  - {}\n".format(i))
    with open(str(args.config),"r") as f:
        for l in f:
            l = l.rstrip('\n')
            fw.write(l + '\n')

            
                        



