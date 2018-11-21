import sys
from Bio import SeqIO
import argparse
import random

## this function is needed to restrict the two float arguments to the required range
def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

## add the positional arguments
parser = argparse.ArgumentParser(description='Given a fastq file and a proportion, return a file with with that proportion of all reads randomly selected.')
parser.add_argument("seq_file", help="file containing the nanopore reads, in fastq format",type=str)
parser.add_argument("fraction",metavar="[0.0-1.0]",help="Proportion of resads to randomly select, from 0 to 1. ",type=restricted_float)
parser.add_argument("out_file", help="output file containing the randomly selected nanopore reads, in fastq format",type=str)
args = parser.parse_args()

## read in the fastq file
records = []
for record in SeqIO.parse(args.seq_file,'fastq'):
    records.append(record)
num_records = len(records)

## determine the number of reads to keep
reads_to_keep = int(args.fraction*num_records)

## randomly add the number of reads we need
downsample = []
seen_indexes = []
for i in range(reads_to_keep):
    index = random.randint(0,num_records - 1)
    while(index in seen_indexes):
        index = random.randint(0,num_records - 1)
    downsample.append(records[index])
    seen_indexes.append(index)

SeqIO.write(downsample,args.out_file,'fastq')
