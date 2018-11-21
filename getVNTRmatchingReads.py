from nano_vntr import nanoVNTR
import sys
from Bio import SeqIO
import argparse

## add the positional arguments
parser = argparse.ArgumentParser(description='Get the reads that match the flanking regions and the VNTR region from a file with nanopore reads in fastq format.')
parser.add_argument("upstream_seq", help="upstream sequence flanking the VNTR",type=str)
parser.add_argument("downstream_seq", help="downstream sequence flanking the VNTR",type=str)
parser.add_argument("repeat_seq", help="sequence of a single repeat unit in the VNTR",type=str)
parser.add_argument("seq_file", help="file containing the nanopore reads, in fastq format",type=str)
parser.add_argument("out_file", help="output file containing the randomly selected nanopore reads, in fastq format",type=str)

## this function is needed to restrict the two float arguments to the required range
def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

## add the optional arguments
parser.add_argument("-cpu","--processors",default=1,choices = range(1,65),metavar="[1-64]",help="number of CPUs to use. Default = 1",type=int)
parser.add_argument("-fc","--flank_cutoff",default=0.75,help="threshold for accepting fuzzy matches of up and downstream sequences flanking the VNTR, from 0 to 1. Default = 0.75",type=restricted_float)
parser.add_argument("-rc","--repeat_cutoff",default=0.9,help="threshold for accepting a VNTR based on the portion comprised by the unit repeat sequence, from 0 to 1. Default = 0.9",type=restricted_float)
parser.add_argument("-pad","--padding",default=1,help="number of extra repeats to add to the repeat count in the VNTR, like a bias term. Default = 1",type=int)
args = parser.parse_args()

## find the number of VNTR repeats in each read
nv = nanoVNTR(args.upstream_seq,args.downstream_seq,args.repeat_seq,args.flank_cutoff,args.repeat_cutoff,args.padding)
nv.get_read_records_from_file(args.seq_file,'fastq')
nv.multiprocess_read_repeats(args.processors)

## get the records for which there are countable VNTRs, and print those as output
records = nv.get_all_countable_VNTR_reads()
SeqIO.write(records,args.out_file,'fastq')
