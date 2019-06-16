
# coding: utf-8

# In[1]:

from pprint import pprint

import gzip
import tre
import random
import difflib
import collections
import fractions
from fractions import Fraction

# http://biopython.org/
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import time
from time import time

import re

# pip3 install python-Levenshtein==0.12.0
import Levenshtein
from Levenshtein import distance

def test1(lim):
    nsegments = 0
    print('''starting test1''')
    with gzip.open("wt-lung_R1_001.fastq.gz", "rt") as handle:
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


def search_sines(sine_f, r1_f, override = 0, upper_mut_dist = 20, step_print = 10000, nlines = 500000, sine_l = 80):

    print ('override =',override)
    sine_set = []
    stats = collections.Counter()

    global bar_codes
    bar_codes = {}
    
    global detailed_stats
    detailed_stats = collections.Counter()
    
    global distances_from_combined_regexp
    distances_from_combined_regexp = {}

    matcher = difflib.SequenceMatcher()
    
    for sine_record in SeqIO.parse(sine_f, "fasta"):
        cur_seq = Seq(str(sine_record.seq)[:sine_l], IUPAC.IUPACAmbiguousDNA())
        cur_seq_rc = cur_seq.reverse_complement()
        sine_set.append(str(cur_seq))
        sine_set.append(str(cur_seq_rc))
        print(cur_seq, cur_seq_rc, '''\n ======================''')

    complete_regexp = '''|'''.join(sine_set)
    p = tre.compile(complete_regexp, tre.EXTENDED)

    if override == 1:
        bases = ['A','C','G','T']
        ind_list = [random.randrange(4) for i in range(sine_l)]
        r_sine = ''.join( [bases[ind_list[i]] for i in range(sine_l)] )
        r_sine_rc = ''.join( [bases[3-ind_list[i]] for i in range(sine_l)] )
        sine_set = [r_sine, r_sine_rc]
        complete_regexp = '''|'''.join(sine_set)
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
            
        complete_regexp = '''|'''.join(sine_set)
        p = tre.compile(complete_regexp, tre.EXTENDED)     

        
    total = 0
    cnt = 0
    start_time = time()
    print('''sequences = ''')

    bar_code_len = 60                         
    for cur_seq in r1_f:
        total += 1
        m = p.search(cur_seq, tre.Fuzzyness(maxerr = upper_mut_dist))
        if m:
            res = m.group(0)
            d = m.cost
            # Filter out strings that were cut out. Approximate by max-length matches
            # 10 is arbitrary, not very small
            
            if (m.groups()[0][1] < len(cur_seq) - 10) and (m.groups()[0][0] > 40):
                # print(m.groups(), len(cur_seq))
                cnt += 1      
                stats[d] += 1

                bar_code = cur_seq[m.groups()[0][0] - 40 : m.groups()[0][0]]

                if bar_code in bar_codes:
                   bar_codes[bar_code] +=  1
                else:
                    bar_codes[bar_code] = 1

            detailed_stats[res] += 1
            distances_from_combined_regexp[res] = d 

        if (total % step_print == 0 or total == nlines):
            print('''distances for first''', total, '''segments \n''')
            print('''========================''')
            print('''time elapsed''', (time() - start_time)/60.0, '''minutes''')
            for k in sorted(stats):
                print('edit distance =', k, 'matches =', stats[k], '''/''',cnt)
        
        if (total == nlines):
            break

##    print('''returning with nlines =''', nlines)
##    print('''detailed stats are''')    
##      for k in sorted(detailed_stats):
##       if detailed_stats[k][1] <= upper_mut_dist:
##            print(k, detailed_stats[k])


def verify(frac, res):
    pref_len = min(res[0], res[1])
    return (frac >= 0.8) and (res[1] - pref_len > 20)
    

def search_sines2(sine, r1_f, thresh = 8, pref = 60, step_print = 100000, nlines = 200000):

    global stats
    stats = {}
     
    
    sine = sine[:pref]
    matcher = difflib.SequenceMatcher(isjunk=None, a=sine,b='',autojunk = False) 
    
    total = 0
    start_time = time()
    print('''sequences for sine = ''')
                
    for cur_seq in r1_f:
        total += 1
        matcher.set_seq2(cur_seq)
        res = matcher.find_longest_match(0, len(sine), 0, len(cur_seq))
        pref_len = min(res[0], res[1])
        d = res[2]
        comp_1 = sine[res[0] - pref_len: res[0] + d]
        comp_2 = cur_seq[res[1] - pref_len: res[1] + d]          
        d1 = distance(comp_1, comp_2)

        stats.setdefault(d, [0, collections.Counter()])
        stats[d][0] += 1
        stats[d][1][Fraction(d1, pref_len + d)] += 1 
            
                        
        if (total % step_print == 0 or total == nlines):
            print('''distances for first''', total, '''segments \n''')
            print('''========================''')
            print('''time elapsed''', (time() - start_time)/60.0, '''minutes''')
            for k in sorted(stats):
                print('longest common =', k, 'num matches =', len(stats[k][1]), '''/''',total)
                if (total >= nlines) and (k >= thresh):
                    n = 0
                    for i in stats[k][1]:
                        if verify(stats[k][1][i], res):
                            n += 1
                     #   else:
                     #       if (stats[k][1][i] < 0.8):
                     #           print ('Fraction = ',stats[k][1][i])
                    print('number of >0.8 instances is',n)
                    #    print ('fraction ',i,' appears ',stats[k][1][i],' times')
            
        if (total == nlines):
            break



'''|'''.join(['''foo''', '''bar''','''ooki'''])


# In[5]:


def fastq_gz_strings(filename):
    with gzip.open(filename, "rt") as handle:    
        for r in SeqIO.parse(handle, "fastq"):
            yield str(r.seq)
             
def gz_strings(filename):
    with gzip.open(filename, "rt") as handle:    
        for line in handle:
            if line[0] not in '''@+#''':  # skip fastq headers/quality
                yield line

def get_sines(sine_f):
    for (i,sine_record) in enumerate(SeqIO.parse(sine_f, "fasta")):
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


def get_min_stats(bar_codes):
    start_time = time()
    distances = {}
    for (i,g) in enumerate(bar_codes):
        min = 40
        for (j,h) in enumerate(bar_codes):
            d = distance(g,h)
            if (i != j) and (d < min):
                min = d
        if min in distances:
            distances[min] += 1
        else:
            distances[min] = 1

    cnt = sum(distances[i] for i in distances)
    print('distance stats are ::::::::::::::::::::::: ')
    for k in sorted(distances):
        print('edit distance =', k, 'matches =', distances[k], '''/''',cnt)
    print('''time elapsed''', (time() - start_time)/60.0, '''minutes''')    
        
                
print('==============================================================================')

for sine in get_sines("B1.fasta"):
    search_sines2(sine,fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''))

#get_min_stats(bar_codes)

#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 1)
#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 2)
#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 3)
#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 4)
#search_sines("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''), 26)
 
#search_sines2("mouse SINEs.fasta",fastq_gz_strings('''wt-lung_R1_001.fastq.gz'''))

# In[26]:

# TODO: DO NOT USE?  20% faster but different results
#search_sines("mouse SINEs.fasta", gz_strings('''old_lung_R2_001.fastq.gz'''), 10000)
 

# # End of code, Saved results below

# In[ ]:


# Test results for 40, r1_old, 1000000
'''0 3 / 5994000
1 41 / 5994000
2 344 / 5994000
3 880 / 5994000
4 1214 / 5994000
5 1438 / 5994000
6 1549 / 5994000
7 1575 / 5994000
8 1794 / 5994000
9 2312 / 5994000
10 3458 / 5994000
11 4386 / 5994000
12 5159 / 5994000
13 6969 / 5994000
14 11904 / 5994000
15 32205 / 5994000
16 106350 / 5994000
17 313644 / 5994000
18 728476 / 5994000
19 1280632 / 5994000
20 1449496 / 5994000
21 1159602 / 5994000
22 623602 / 5994000
23 199311 / 5994000
24 43522 / 5994000
25 9264 / 5994000
26 2935 / 5994000
27 1247 / 5994000
28 480 / 5994000
29 146 / 5994000
30 31 / 5994000
31 16 / 5994000
32 4 / 5994000
33 6 / 5994000
34 5 / 5994000'''

# Test results for 40, r2_old, 1000000
'''0 9 / 6000000
1 41 / 6000000
2 231 / 6000000
3 724 / 6000000
4 1044 / 6000000
5 1318 / 6000000
6 1469 / 6000000
7 1527 / 6000000
8 1651 / 6000000
9 2176 / 6000000
10 3229 / 6000000
11 4060 / 6000000
12 5038 / 6000000
13 6622 / 6000000
14 11872 / 6000000
15 32432 / 6000000
16 107532 / 6000000
17 317798 / 6000000
18 738101 / 6000000
19 1281388 / 6000000
20 1452101 / 6000000
21 1156272 / 6000000
22 615449 / 6000000
23 199331 / 6000000
24 43301 / 6000000
25 9528 / 6000000
26 3194 / 6000000
27 1554 / 6000000
28 501 / 6000000
29 208 / 6000000
30 81 / 6000000
31 38 / 6000000
32 48 / 6000000
33 78 / 6000000
34 54 / 6000000'''

'''time elapsed  89.94379512866338 wt-r1, 1000000
0 9 / 6000000
1 79 / 6000000
2 428 / 6000000
3 1104 / 6000000
4 1554 / 6000000
5 1845 / 6000000
6 1885 / 6000000
7 1744 / 6000000
8 1984 / 6000000
9 2717 / 6000000
10 3987 / 6000000
11 4970 / 6000000
12 6004 / 6000000
13 7561 / 6000000
14 13355 / 6000000
15 35661 / 6000000
16 117950 / 6000000
17 343914 / 6000000
18 779759 / 6000000
19 1310694 / 6000000
20 1441609 / 6000000
21 1101603 / 6000000
22 571210 / 6000000
23 182370 / 6000000
24 43364 / 6000000
25 11148 / 6000000
26 5102 / 6000000
27 3058 / 6000000
28 1677 / 6000000
29 985 / 6000000
30 421 / 6000000
31 141 / 6000000
32 55 / 6000000
33 38 / 6000000
34 14 / 6000000
35 1 / 6000000'''

'''time elapsed  15.340233178933461 old-r2, 30000000
0 27 / 24000000
1 257 / 24000000
2 1431 / 24000000
3 4034 / 24000000
4 5505 / 24000000
5 6258 / 24000000
6 6354 / 24000000
7 6010 / 24000000
8 5526 / 24000000
9 5572 / 24000000
10 5902 / 24000000
11 5906 / 24000000
12 5881 / 24000000
13 5863 / 24000000
14 7052 / 24000000
15 11016 / 24000000
16 21563 / 24000000
17 41693 / 24000000
18 61513 / 24000000
19 64138 / 24000000
20 45925 / 24000000
21 22032 / 24000000
22 6930 / 24000000
23 1498 / 24000000
24 204 / 24000000
25 39 / 24000000
26 12 / 24000000
28 1 / 24000000'''

# ==================
'''time elapsed  19.913855942090354 1000000 young r2
0 8 / 6000000
1 190 / 6000000
2 623 / 6000000
3 1446 / 6000000
4 1821 / 6000000
returning with nlines =  1000000'''

#=============
'''time elapsed  20.19312702814738 1000000 old r2
0 9 / 6000000
1 110 / 6000000
2 419 / 6000000
3 973 / 6000000
4 1259 / 6000000
returning with nlines =  1000000'''

'''time elapsed 16.6222222050031 minutes
0 9 / 6000000
1 110 / 6000000
2 419 / 6000000
3 971 / 6000000
4 1258 / 6000000
returning with nlines = 1000000'''

