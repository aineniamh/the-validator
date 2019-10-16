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
    
    num_muts = int(distance*len(seq))

    indexes_to_mutate = random.sample(range(len(seq)), num_muts)
    sorted_indexes = sorted(indexes_to_mutate)
    
    new_seq = ''
    for c,a in enumerate(seq):
        if c in sorted_indexes:
            new_mut = [i for i in mutations if i!=a]
            new_seq += random.choice(new_mut)
        else:
            new_seq += a
    return new_seq, num_muts

for record in SeqIO.parse(str(args.fasta),"fasta"):
    header = "Mutant_of_"+record.description
    header += " distance={}".format(args.distance)
    
    new_seq,num_muts = mutate(record.seq,float(args.distance))
    header += " num_mutants={}".format(num_muts)

    record.description = header
    with open(str(args.outfile),"w") as f:
        f.write(">{}\n{}\n".format(record.description, new_seq))

