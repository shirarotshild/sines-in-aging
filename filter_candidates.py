#!/usr/bin/env python3.6

import sys
import tre

import Bio
from Bio import SeqIO
from Bio.Seq import Seq

import gene_lib
from gene_lib import log
from gene_lib import get_sine_forward

def filter_potential_sines(in_fname, sine_fname_list, sine_header=67, maxerr=19):
    """
    Finds candidate SINEs with a certain distance from a prefix length.
    To be used for preliminary screening (input for later steps).
    """
    regexp_parts = []
    for sine_fname in sine_fname_list:
        cur_sine = gene_lib.get_sine_forward(sine_fname)[:sine_header]
        cur_sine_rc = str(Seq(cur_sine).reverse_complement())
        regexp_parts += [cur_sine, cur_sine_rc]
    regexp = '|'.join(regexp_parts)
    
    re = tre.compile(regexp, tre.EXTENDED)
    fuzziness = tre.Fuzzyness(maxerr=maxerr)
    
    with gene_lib.open_any(in_fname, 'rt') as in_file_handle:
        records = SeqIO.parse(in_file_handle, format="fastq")

        for rec in records:
            match = re.search(str(rec.seq), fuzziness)
            if match:
                # log(rec.seq)
                #sine_location = match.groups() #returns tuple of tuples (in this case: ((2,78), ) for example
                SeqIO.write(rec, sys.stdout, 'fasta')

# Writes to stdout, uncompresed .fasta
[*sine_fname_list, merged_input_fname] = sys.argv[1:]
log(f'About to screen {merged_input_fname} for any of {sine_fname_list}, direct or reverse-complement.')
filter_potential_sines(merged_input_fname, sine_fname_list)
