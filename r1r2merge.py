#!/usr/bin/env python3
import math
import tre
import sys
import sines_io
from sines_io import fastq_zst_records
from Bio.SeqRecord import SeqRecord
import math
from math import ceil

import os

# http://biopython.org/
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#current file generated with 15, 4
common_req = 20
init_err = 3
tot_err = 3
padding = 20
step = 10000

def merged_paired_ends(records1, records2):
    tot_good = 0
    tot_great = 0
    tot = 0
    print('in merged_paired_ends',records1,records2)
    for (rec1,rec2) in zip(records1, records2):
        tot += 1 
        str1 = str(rec1.seq)
        str2 = str(rec2.seq.reverse_complement())
        end1 = str1[-common_req:]
        re = tre.compile(end1, tre.EXTENDED)
        # we expect small errors here
        res_seq = None
        match = re.search(str2, tre.Fuzzyness(maxerr = init_err))
        if match:
            tot_good += 1
            match_loc = match.groups()[0][0]
            to_search_len = match_loc + common_req
            fuzzyness = max(tot_err, ceil(0.1*to_search_len))
            re = tre.compile(str1[-to_search_len :], tre.EXTENDED)
            match_tot = re.search(str2, tre.Fuzzyness(maxerr = fuzzyness))

            if match_tot:
            #    if (tot_good % 100 == 0):
            #        print('fuzzyness = ', fuzzyness)
                tot_great += 1
                # An arbitrary decision: take the common string from r2
                res_str = end1[ : -to_search_len] + str2 
                res_seq = SeqRecord(Seq(res_str), id=rec1.id, name=rec1.name, description = rec1.description, letter_annotations =
                                    {"phred_quality":[30 for i in range(len(res_str))]})
                if (tot_great % step == 0):
                    print('nicely matched ',str1,'\n',str2, to_search_len, match_tot.group(0), match.group(0), match_tot.cost, match.cost)
                yield res_seq

        res_str = str1 + ''.join(['N' for i in range(padding)]) + str2
        res_seq = SeqRecord(Seq(res_str), id=rec1.id, name=rec1.name,
                   description = rec1.description, letter_annotations = {"phred_quality":[30 for i in range(len(res_str))]}) 
        if  (tot % step == 0):
            print (tot, tot_good, tot_great)
           # print('matched ',str1,'\n',str2, len(str1), len(str2))
        yield res_seq
                

                
            
### MAIN ###
[base1, base2, part, out_base] = sys.argv[1:]
def write_merged():
    out_fname = f'{out_base}.part{part}e8.merged.fastq'
    in_fname1 = f'{base1}.part{part}e8.fastq.zst'
    in_fname2 = f'{base2}.part{part}e8.fastq.zst'
    # We assume for now all records match. Generally, need to verify.
    print('About to merge ',in_fname1, in_fname2)
    records1 = fastq_zst_records(in_fname1)
    records2 = fastq_zst_records(in_fname2)
    merged = merged_paired_ends(records1, records2)
    print('about to write ', out_fname + '.tmp')
    Bio.SeqIO.write(merged, out_fname + '.tmp', 'fastq')
    os.rename(out_fname + '.tmp', out_fname)    
        
write_merged()
