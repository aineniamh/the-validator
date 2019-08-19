import random
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Mutate a fasta sequence to a set distance.')
parser.add_argument("--fasta", action="store", type=str, dest="fasta")
parser.add_argument("--distance", action="store", type=float, dest="distance")
parser.add_argument("--outfile", action="store", type=str, dest="outfile")


args = parser.parse_args()

def mutate(seq, distance):
    mutations = ['A', 'T', 'C', 'G']
    i = random.randint(0, len(seq)-1)
    val = random.choice(mutations)
    indexes_to_mutate = random.sample(range(len(seq)), int(distance*len(seq)))

    sorted_indexes = sorted(indexes_to_mutate)

    return ''.join(val if c in sorted_indexes else a for c, a in enumerate(seq))

for record in SeqIO.parse(str(args.fasta),"fasta"):
    header = record.description
    header += "distance={}".format(args.distance)
    
    new_seq = mutate(record.seq,float(args.distance))
    with open(str(args.outfile),"w") as f:
        f.write(">{}\n{}\n".format(record.description, new_seq))

