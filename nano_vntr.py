import sys
from Bio import SeqIO
import difflib
import re
from collections import defaultdict
from sklearn.cluster import KMeans
import numpy as np
from multiprocessing import Process,Manager

class nanoVNTR:

    def __init__(self,LF,RF,REP,CUTOFF,REP_CUTOFF,PADDING):
        self.LF = LF # sequence flanking the repeat region on the left side
        self.RF = RF # sequence flanking the repeat region on the right side
        self.REP = REP # single repeat sequence
        self.CUTOFF = CUTOFF # threshold for accepting a fuzzy match
        self.REP_CUTOFF = REP_CUTOFF # threshold for what percentage of the rep region needs to conform to rep multiples
        self.PADDING = PADDING # number of extra repeats to add to the repeat count in the VNTR, like a bias term
        self.READ_RECS = [] # list of read records, in Bio::SeqRecord format
        self.REP_DICT = {} # stores read VNTR information
        self.NUM_ALLELES_BIAS = 0.3 # used when estimating whether there are 1 or 2 alleles

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

    ## given a fasta or fastq file of reads, stores a list of read records
    def get_read_records_from_file(self,read_file,file_format):
        self.READ_RECS = [] # to ensure idempotence
        for record in SeqIO.parse(read_file,file_format):
            self.READ_RECS.append(record)

    ## given a list of reads in Bio::SeqRecord format and a global idct, updates the dict with read ID -> output from get_single_read_repeats
    def get_read_repeats_for_seqrecords(self,record_list,return_dict):
        results = {}
        for record in record_list:
            repeat = self.get_single_read_repeats(record)
            if repeat is not None:
                results[record.id] = repeat
        return_dict.update(results)

    ## splits all the read records into multiple lists with one per CPU, given the number of CPUs to use, and updates
    ## the self.REP_DICT with all of the results
    def multiprocess_read_repeats(self,cpu_count):
        splits = np.array_split(self.READ_RECS,cpu_count)
        processes = []
        m = Manager()
        return_dict = m.dict() 
        for i in range(len(splits)):
            p = Process(target=self.get_read_repeats_for_seqrecords, args=(splits[i],return_dict))
            processes.append(p)
            p.start()
        for one_process in processes:
            one_process.join()
        self.REP_DICT = dict(return_dict)

    ## returns a list of records in which countable VNTRs were found
    def get_all_countable_VNTR_reads(self):
        reads = []
        for i in range(len(self.READ_RECS)):
            if (self.READ_RECS[i].id in self.REP_DICT) and (self.REP_DICT[self.READ_RECS[i].id]['REP_COUNT'] is not None):
                reads.append(self.READ_RECS[i])
        return reads

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
        
    ## estimates the number of alleles (1 or 2) in a list of repeat counts, based on analysis of variance
    def estimate_num_alleles(self,read_count_list):
        sil_scores = {}
        X = np.array(read_count_list).reshape(-1,1)
        km = KMeans(n_clusters=2)
        km.fit(X)
        clusters = defaultdict(list) 
        for i in range(len(km.labels_)):
            clusters[km.labels_[i]].append(read_count_list[i])
        # if the total variance is higher when we treat the list as two clusters than when we don't
        # we will assume there is only 1 allele. Otherwise we assume 2.
        two_allele_var = sum([np.var(clusters[i]) for i in clusters])
        if two_allele_var > (np.var(read_count_list) - self.NUM_ALLELES_BIAS):
            return 1
        return 2
        
    ## given a list of repeat counts, an expected number of alleles (1 <= k <= 2) and a measure of average, 
    ## returns k VNTR counts
    def get_VNTR_alleles(self,read_count_list,average_measure):
        num_alleles = self.estimate_num_alleles(read_count_list)
        # handle the special case where there is only one allele
        if num_alleles == 1:
            allele = round(np.mean(read_count_list))
            return [allele,allele]
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

