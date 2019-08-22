import os
import argparse
from Bio import SeqIO
barcodes = []

for r,d,f in os.walk("/Users/aineniamh/Documents/artic-noro/simulated_test/broken_reads"):
    for filename in f:
        if filename.endswith("fastq"):
            records = []
            barcode = filename.rstrip(".fastq").lstrip("barcode_")
            fixed = barcode.replace("_","-")
            for record in SeqIO.parse(r + '/' + filename, "fastq"):
                header = record.description
                # print(barcode, fixed)
                new_header = header.replace(barcode,fixed)
                # print(new_header)
                record.description=new_header
                records.append(record)
            
            new_file = "/Users/aineniamh/Documents/artic-noro/simulated_test/annotated_reads/barcode_" + fixed +".fastq"
            with open(new_file, "w") as fw:
                SeqIO.write(records,fw,"fastq")
            

