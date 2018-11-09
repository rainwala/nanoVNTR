import sys
from Bio import SeqIO
import difflib
import re
from collections import defaultdict
from fuzzywuzzy import fuzz
from sklearn.cluster import KMeans
import numpy as np

class nanoVNTR:

    def __init__(self,LF,RF,REP,CUTOFF,REP_CUTOFF,PADDING):
        self.LF = LF # sequence flanking the repeat region on the left side
        self.RF = RF # sequence flanking the repeat region on the right side
        self.REP = REP # single repeat sequence
        self.CUTOFF = CUTOFF # threshold for accepting a fuzzy match
        self.REP_CUTOFF = REP_CUTOFF # threshold for what percentage of the rep region needs to conform to rep multiples
        self.PADDING = PADDING # number of extra repeats to add to the repeat count in the VNTR, like a bias term
        self.REP_DICT = {} # stores read VNTR information

    ## gets the positions of fuzzy matches at least as good as the threshold of a query string within a larger string 
    ## provided that the repeat makes up at least cutoff proportion of the string
    def get_matches(self,seq, query, threshold):
        matches = {}
        s = difflib.SequenceMatcher(None, seq, query)
        for i, j, n in s.get_matching_blocks():
            match = ''.join(seq[i:i+n])
            if len(match) / float(len(query)) >= threshold:
                matches[i] = match
        return matches

    ## gets the number of occurrences of a repeat sequence in a string 
    def get_region_rep_count(self,rep_region,rep,cutoff):
        if rep not in rep_region:
            return None
        occurences = [m.start() for m in re.finditer(rep, rep_region)]
        if len(occurences)/(len(rep_region)/len(rep)) < cutoff:
            return None
        trim = rep_region[occurences[0]:occurences[-1]+len(rep)]
        return (len(trim)/3) + self.PADDING

    ## given a Bio::SeqRecord of a single read, gets the number of repeats 
    ## in the VNTR region, as well as the left and right flanking sequence and VNTR seqeunce, 
    ## or returns none if the flanking regions are not found to occur only once
    ## in a single orientation, with the intervening sequence harbouring a high 
    ## enough proportion of repeats
    def get_single_read_repeats(self,record):
        ## try to find a fuzzy match for both the left and right flank in either orientation
        read_seqs = [record.seq,record.seq.reverse_complement()]
        for seq in read_seqs:
            l_matches = self.get_matches(str(seq),self.LF,self.CUTOFF)
            if len(l_matches) == 1:
                r_matches = self.get_matches(str(seq),self.RF,self.CUTOFF)
                if len(r_matches) == 1:
                    l_index = list(l_matches.keys())[0] + len(l_matches[list(l_matches.keys())[0]])
                    r_index = list(r_matches.keys())[0]
                    mid = str(seq[l_index:r_index])
                    # find the first AND last instance of the REP in the mid sequence, and trim accordingly
                    num_reps = self.get_region_rep_count(mid,self.REP,self.REP_CUTOFF)
                    return {'LF': l_matches, 'RF': r_matches, 'VNTR': mid, 'REP_COUNT': num_reps}
        return None

    ## given a fasta or fastq file of an alignment, returns a dict of read ID -> output from get_single_read_repeats
    def get_all_read_repeats(self,read_file,file_format):
        self.REP_DICT = {} # to ensure idempotence
        for record in SeqIO.parse(read_file,file_format):
            repeat = self.get_single_read_repeats(record)
            if repeat is not None:
                self.REP_DICT[record.id] = repeat
    
    ## prints all read summaries for reads in which VNTRs were found
    def print_read_repeat_summaries(self):
        for read_id in sorted(self.REP_DICT):
            outline = ""
            outline += read_id + "\t"
            for key in ['LF','RF','VNTR','REP_COUNT']:
                outline += key + ': ' + str(self.REP_DICT[read_id][key]) + '\t'
            print(outline.rstrip('\t'))

    ## converts the read_repeat dict to a list of repeat counts for allele finding
    def convert_read_repeat_dict_to_repeat_count_list(self):
        return [self.REP_DICT[read_id]['REP_COUNT'] for read_id in self.REP_DICT if self.REP_DICT[read_id]['REP_COUNT'] is not None]
        
    ## estimates the number of alleles in a list of repeat counts, based on variance

    ## given a list of repeat counts, an expected number of alleles (1 <= k <= 2) and a measure of average, 
    ## returns k VNTR counts
    def get_VNTR_alleles(self,read_count_list,num_alleles,average_measure):
        km = KMeans(n_clusters=num_alleles)
        km.fit(np.array(read_count_list).reshape(-1,1))
        if average_measure == 'centroid':
            return sorted([round(cc[0]) for cc in km.cluster_centers_])
        else:
            clusters = defaultdict(list)
            for i in range(len(km.labels_)):
                clusters[km.labels_[i]].append(read_count_list[i])
            if average_measure == 'mean':
                return sorted([round(np.mean(clusters[l])) for l in clusters])
            elif average_measure == 'median':
                return sorted([round(np.median(clusters[l])) for l in clusters])
        return None

