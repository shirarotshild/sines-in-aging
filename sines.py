
# coding: utf-8

# In[1]:

from pprint import pprint

import collections
import difflib
from fractions import Fraction
import gzip
import io
import itertools
import os
import random
import re
import sys
from time import time

import sines_io
from sines_io import fastq_gz_strings
from sines_io import gz_strings
from sines_io import fastq_zst_strings
from sines_io import fastq_zst_records
from sines_io import zst_strings

import tre


# http://biopython.org/
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# pip3 install python-Levenshtein==0.12.0
import Levenshtein
from Levenshtein import distance

def test1(lim):
    nsegments = 0
    print('''starting test1''')
    with gzip.open("Young-lung/wt-lung_R1_001.fastq.gz", "rt") as handle:
        for r in SeqIO.parse(handle, "fastq"):
            print(r.seq,'''\n==============''')
            nsegments += 1
            if nsegments == lim:
                return



# DEPRECATED
def closest_sub(long_str, short_str):
    d = distance(long_str, short_str)
    return d - (len(long_str) - len(short_str))
    
def closest_sub_exact(long_str, short_str):
    min_d = 100
    if (len(long_str) < len(short_str)):
        print('''length error''')
        # return a very large value
        return 100
    nsteps = len(long_str) - len(short_str) + 1
    for i in range(nsteps):
        cur_substr = long_str[i : i+len(short_str)]
        d = distance(cur_substr, short_str)
        if (d < min_d):
            min_d = d
    return min_d        

def cutoff(long_seq, x, max_dist):
    lx = len(x)
    ll = len(long_seq)
    l_sub = int(lx / (max_dist + 1))
    for i in range(max_dist+1):
        # TODO: last one does not need to be of length l_sub
        cur_sub = x[i*l_sub: (i+1)*l_sub]
        if re.search(cur_sub, long_seq):
            return 1
    return 0    


start_time = None


# TODO: translate 'N' to r'[AGCT]' etc. for more accurate distances.
# https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
# https://biopython.org/DIST/docs/api/Bio.Data.CodonTable.AmbiguousCodonTable-class.html
#def iupac_notation_to_regexp(iupac_string):


def filter_potential_sines(records, sine_pattern, sine_header = 67, maxerr = 14):
    re = tre.compile(sine_pattern[:sine_header], tre.EXTENDED)
    fuzziness = tre.Fuzzyness(maxerr=maxerr)
    
    for rec in records:
        match = re.search(str(rec.seq), fuzziness)
        if match:
            yield rec


def write_filtered(base='Young-lung/wt-lung_R2_001', parts=range(30)):
    [B1] = get_sines_forward("B1.fasta")
    for i in parts:
        # Before running, copy/symlink input files into current dir:
        # ln -s /media/anat/2Tincognita/anti-aging/*-lung/*part* .
        in_fname = f'{base}.part{i}e8.fastq.zst'
        out_fname = f'{base}.part{i}e8.potential_B1_head67err14.fastq'
        records = fastq_zst_records(in_fname)
        filtered = filter_potential_sines(records, B1)
        print(f"Writing {out_fname}...")
        # Rename when done to easily distinguish unfinished files.
        SeqIO.write(filtered, out_fname + '.tmp', 'fastq')
        os.rename(out_fname + '.tmp', out_fname)


# TODO: Need to rewrite to get the sine + a barcode before (direct) / after (reversed-complemented)
# For now, we want bar-code
def search_sines(sines, r1_f, override = 0, upper_mut_dist = 20, step_print = 1000000, nlines = 100000000, sine_l = 70):

    print ('override =',override)
    sine_set = []
    stats = collections.Counter()

    global bar_codes
    bar_codes = {}
    
    global detailed_stats
    detailed_stats = collections.Counter()
    
    global distances_from_combined_regexp
    distances_from_combined_regexp = {}

    complete_regexp = '''|'''.join([sine[:sine_l] for sine in sines])
    p = tre.compile(complete_regexp, tre.EXTENDED)

    if override == 1:
        bases = ['A','C','G','T']
        ind_list = [random.randrange(4) for i in range(sine_l)]
        r_sine = ''.join( [bases[ipnd_list[i]] for i in range(sine_l)] )
        r_sine_rc = ''.join( [bases[3-ind_list[i]] for i in range(sine_l)] )
        sine_set = [r_sine, r_sine_rc]
        complete_regexp = '|'.join(sine_set)
        p = tre.compile(complete_regexp, tre.EXTENDED)

    # Also specifies the shift  range   
    if override > 1:
        if override > 2:
            d = override - 1 #random.randrange(2, override)
            print('skipping ',d)
            for (i,cur_seq) in enumerate(r1_f):
                if i == d:
                    break
                
        sine_set = []
        for (i,s) in enumerate(r1_f):
            cur_seq = Seq(s[:sine_l], IUPAC.IUPACAmbiguousDNA())
            cur_seq_rc = cur_seq.reverse_complement()
            sine_set.append(str(cur_seq))
            sine_set.append(str(cur_seq_rc))
            if i == 2:
                break
            
        complete_regexp = '|'.join(sine_set)
        p = tre.compile(complete_regexp, tre.EXTENDED)     

        
    total = 0
    cnt = 0
    start_time = time()
    print('''sequences = ''')

   # bar_code_len = 60                         
    for cur_seq in r1_f:
        total += 1
        m = p.search(cur_seq, tre.Fuzzyness(maxerr = sine_l - 10))
        if m:
            res = m.group(0)
            d = m.cost
            # Filter out strings that were cut out. Approximate by max-length matches
            # 10 is arbitrary, not very small
            # barcodes are not in place here

            
            stats[d] += 1
            bar_code_min_len = 23
          #  if (m.groups()[0][1] < len(cur_seq) - 5) and (d <= upper_mut_dist):
          
            if (m.groups()[0][0] >= bar_code_min_len) and (d <= upper_mut_dist):
                cnt += 1      
                detailed_stats[res] += 1
                bar_code = cur_seq[m.groups()[0][0] - bar_code_min_len : m.groups()[0][0]]
                
                bar_codes.setdefault(bar_code, 0)
                bar_codes[bar_code] +=  1

            
        #    distances_from_combined_regexp[res] = d 

        if (total % step_print == 0) or (total == nlines):
            print('''stats for first''', total, '''segments \n''')
            print('''========================''')
            print('''time elapsed''', (time() - start_time)/60.0, '''minutes''')
            
            for k in sorted(stats):
                print('edit distance =', k, 'matches =', stats[k], '''/''',cnt)
         #   pprint.pprint(collections.Counter(detailed_stats.values()))
        
        if (total == nlines):
            return bar_codes



def verify(frac, res, f1 = 1, pref_bound = 20):
    pref_len = min(res[0], res[1])
    return (frac <= f1) and (res[1] - pref_bound > 20)
  

# Looks for barcodes that have a neighborhood of size d=4 (arbitrary) which is empty.
# For each barcode we get the multiplicity
def search_sines2(sine, r1_f, frac_bound, pref_bound, start_line = 0, step_print = 1000000, nlines = 200000000, thresh = 9, pref = 60):

    global stats
    stats = {}
     
    print('step ',step_print, nlines)
    sine = sine[:pref]
    matcher = difflib.SequenceMatcher(isjunk=None, a=sine,b='',autojunk = False) 
    
    total = 0
    cnt = 0
    start_time = time()
    print('''condidates for sine = ''')

    if start_line > 0:
        for (i,cur_seq) in enumerate(r1_f):
            if i == start_line - 1:
                break

                
    for cur_seq in r1_f:
        
        if (total % step_print == 0 or total == nlines):
            print('''distances for first''', total, '''segments \n''')
            print('''========================''')
            print('''time elapsed''', (time() - start_time)/60.0, '''minutes''')
            for k in sorted(stats):
                n = sum([i for i in stats[k][1].values()])
                print('longest common =', k, 'num matches =', n, stats[k][0], '''/''',cnt)
                if (total >= nlines) and (k >= thresh):
                    for (i,frac) in enumerate(sorted(stats[k][1])):
                        print (k,'Fraction = ',frac)
                        if i == 20:
                            break
                        
        if (total == nlines):
            break
        
        total += 1
        matcher.set_seq2(cur_seq)
        res = matcher.find_longest_match(0, len(sine), 0, len(cur_seq))
        com = res[2]

        complete_regexp = sine[: res[0]] + '$'
        p = tre.compile(complete_regexp, tre.EXTENDED)
        max_fuzz = res[0] # int(frac_bound*res[0]) is better perhaps, but want to trivialize it for now
        m = p.search(cur_seq[:res[1]], tre.Fuzzyness(maxcost = max_fuzz,
                                                     delcost = int(1/4.0*max_fuzz)+1,
                                                     inscost = int(1/4.0*max_fuzz)+1))
        if m == None:
            continue
        
        start_p = m.groups()[0][0]
        d = m.cost
            
        # This is the fraction of edit distance out of all.
        # In most cases, this is the right edit distance for the overall prefix

        if (res[0] + com) == 0:
            print('How peculier!','com =',com,'res[0] = ',res[0], m.cost)
            continue
        
        frac = Fraction(d, res[0] + com)

        stats.setdefault(com, [0, collections.Counter()])
        stats[com][0] += 1

        try:
            if (start_p >= pref_bound) and Fraction(d, res[0]) <= frac_bound:
                stats[com][1][frac] += 1
                cnt += 1
        except (ZeroDivisionError):
            pass



'''|'''.join(['''foo''', '''bar''','''ooki'''])


def get_sines_forward(sine_fname):
    """Only in direction given in file."""
    for (i,sine_record) in enumerate(SeqIO.parse(sine_fname, "fasta")):
        cur_seq = Seq(str(sine_record.seq), IUPAC.IUPACAmbiguousDNA())
        yield str(cur_seq)

def get_sines(sine_fname):
    """As given in file + reverse complements."""
    for (i,sine_record) in enumerate(SeqIO.parse(sine_fname, "fasta")):
        cur_seq = Seq(str(sine_record.seq), IUPAC.IUPACAmbiguousDNA())
        yield str(cur_seq)
        cur_seq_rc = cur_seq.reverse_complement()
        yield str(cur_seq_rc)
        print(cur_seq, cur_seq_rc, '''\n ======================''')

print('''Here come the SINES!''')



good_lines = [
    '''TGATTATCAGGTGAGAAATCACGATGGGAATTAAAAGCATTCTGAAGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAAGCAGAGGCAGACGGATTTCTGAATTCGAGGCCAGCCTGGTCTACAGAGTGAGTTCCAGGAC''',
    '''GAATCCTTGTTTTACAGCTGGATACGATGTAGGCTTACAGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGTGGATTTCTGAGTTCGAGGCCAGCCTGGTCTACAAAGTGAGTTCCAGGACAGCCAGG''',
    '''TTTTGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGCGGATTTCAGAGTTTGAGGCCAGCCTGGTCTACAAAGTGAGTTCCAGGACAGCTGGGCTACAGAGAAATCCTGACTTAAAAAAACAAAAACA''',
    '''TTTTGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGCGGATTTCAGAGTTTGAGGCCAGCCTGGTCTACAAAGTGAGTTCCAGGACAGCTGGGCTACAGAGAAAGCCTGACTTAAAAAAACAAAAACA''',
    '''GCAGGTAAGAACCATCAAAGCGACCCTATTAGGTAAATCCTGATAATATTCCATTTTAAAAATGGTGAAAGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGCGGATTTCTGAGTTCGAGGCCAGCCT''',
    '''CATAAGAAAGAGCTGTGCGGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGTGGATTTCTGAGTTCGAGGCCAGCCTGGTCTACAAAGTGAGTTTCAGGACAGCCAGGGCTATACAGAGAAACCCTGT''',
]
#search_sines("mouse SINEs.fasta",good_lines)

 
def get_min_stats2(cur_stats, ok_inter = 23, good_inter = 11, good_d = 4, ok_d = 6):
    # TODO Improve over the complexity here: Idea in mind works if there are only a few close ones for each string. 
    start_time = time()
    global stats
    stats = collections.Counter()
    global good
    good = {}

    matcher = difflib.SequenceMatcher(isjunk=None, a='',b='',autojunk = False)

    #note we actually do count to vector itself as being close to itself
    tot = len(cur_stats)
    for (i,g) in enumerate(cur_stats):      
        if (cur_stats[g] > 1):
            stats[len(g)] += 1
        else:
            max_intersection = (-1,0)
            matcher.set_seq1(g)
            for (j,h) in enumerate(cur_stats):
                if i != j:
                    matcher.set_seq2(h)
                    res = matcher.find_longest_match(0, len(g), 0, len(h))
                    com = res[2]
                    if com > max_intersection[0]:
                        max_intersection = (com, h)
                        
            if (max_intersection[0] <= ok_inter):
                good.setdefault(max_intersection[0], [])
                d = distance(g,max_intersection[1])
                good[max_intersection[0]].append(distance(g,max_intersection[1]))
       
            stats[max_intersection[0]] += 1
            
        if (i % 1000 == 0) or (i >= tot - 1):
            print(i)
            cnt = i
            print('counter stats are ::::::::::::::::::::::: ')
            for j in sorted(stats):
                print('intersection =', j, 'matches =', stats[j], '''/''',i)
            
            print('''time elapsed''', (time() - start_time)/60.0, '''minutes''')

    

def get_min_stats(cur_stats, thresh = 6):
    # TODO Improve over the complexity here: Idea in mind works if there are only a few close ones for each string. 
    start_time = time()
    counters = collections.Counter()
    good = {}

    #note we actually do count to vector itself as being close to itself                               
    for (i,g) in enumerate(cur_stats):
        if (i % 1000 == 0):
            print(i)
        g_stats = collections.Counter()
        for (j,h) in enumerate(cur_stats):
            d = distance(g,h)
            g_stats[d] += cur_stats[h]

        v = sum([g_stats[k] for k in sorted(g_stats) if k <= thresh])
        if (v>0):
           counters[v] += cur_stats[g]
                               

    cnt = sum(counters[i] for i in counters)
    print('counter stats are ::::::::::::::::::::::: ')
    for i in sorted(counters):
        print('frequency =', i, 'matches =', counters[i], '''/''',cnt)
    
    print('''time elapsed''', (time() - start_time)/60.0, '''minutes''')    
        

def upper_level(frac = 0.25, pref_bound = 28, start_line = 0):                
    print('==============================================================================')

    #<<<<<<< Updated upstream
    #for sine in get_sines("B1.fasta"):
    #    search_sines2(sine,fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), frac, pref_bound, start_line)

    global barcode_set
    barcode_set = []
    for sine in get_sines("B1.fasta"):
        bar_codes = search_sines([sine], fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''))
        barcode_set.append(bar_codes)
 #       get_min_stats(detailed_stats)
                                   
#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), step_print=500_000, nlines=10_000_000)
#get_min_stats(bar_codes)
#>>>>>>> Stashed changes

#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 1)
#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 2)
#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 3)
#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 4)
#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 26)
 
#search_sines2("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''))

### MAIN ###
[base, part] = sys.argv[1:]
write_filtered(base, [int(part)])
