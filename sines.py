
# coding: utf-8

# In[1]:


import gzip
import tre

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
    print('starting test1')
    with gzip.open("old_lung_R1_001.fastq.gz", "rt") as handle:
        for r in SeqIO.parse(handle, "fastq"):
            print(r.seq,'\n==============')
            nsegments += 1
            if nsegments == lim:
                return

re.search(r'xzg', 'wxz') 


# In[2]:



p = tre.compile(r'foo|goo', tre.EXTENDED)
print(p)
m = p.search('fststoststosst', tre.Fuzzyness(maxerr=1))
print(m)
m = p.search('fststoststosst', tre.Fuzzyness(maxerr=2))
print(m.group(0))
print(m.groups())
print()
m = p.search('fxxfox', tre.Fuzzyness(maxerr=2))
print(m.group(0))
print(m.cost,m.fuzzyness)


# In[3]:


# DEPRECATED
def closest_sub(long_str, short_str):
    d = distance(long_str, short_str)
    return d - (len(long_str) - len(short_str))
    
def closest_sub_exact(long_str, short_str):
    min_d = 100
    if (len(long_str) < len(short_str)):
        print('length error')
        # return a very large value
        return 100
    nsteps = len(long_str) - len(short_str) + 1
    for i in range(nsteps):
        cur_substr = long_str[i : i+len(short_str)]
        d = distance(cur_substr, short_str)
        if (d < min_d):
            min_d = d
    return min_d        


# In[3]:


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
            
def search_sines(sine_f, r1_f,  upper_mut_dist = 4, step_print = 10000, nlines = 1000000, sine_l = 40,):
    sine_set = []
    stats = {}
    for sine_record in SeqIO.parse(sine_f, "fasta"):
        cur_seq = Seq(str(sine_record.seq)[:sine_l], IUPAC.IUPACAmbiguousDNA())
        cur_seq_rc = cur_seq.reverse_complement()
        sine_set.append(str(cur_seq))
        sine_set.append(str(cur_seq_rc))
        print(cur_seq, cur_seq_rc, '\n ======================')
        
    complete_regexp = '|'.join(sine_set)
    p = tre.compile(complete_regexp, tre.EXTENDED)
    
    cnt = 0
    start_time = time()
    print('sequences = ')
    for cur_seq in r1_f:
        
        cnt += 1
        m = p.search(cur_seq, tre.Fuzzyness(maxerr=upper_mut_dist))
        if m:
            res = m.group(0)
            d = m.cost
            if d in stats:
                stats[d] += 1
            else:
                stats[d] = 1 
                    
        #    if (cutoff(cur_seq, x, 4) > 0):
        #        d = closest_sub_exact(cur_seq, x)
        #        if d in stats:
        #            stats[d] += 1
        #        else:
        #            stats[d] = 1
        if (cnt % step_print == 0 or cnt == nlines):
            print('distances for first', cnt, 'segments \n')
            print('========================')
            print('time elapsed', (time() - start_time)/60.0, 'minutes')
            for k in sorted(stats):
                print (k, stats[k], '/',6*cnt)  
        if (cnt == nlines):
            print('returning with nlines =', nlines)
            return
    


# In[4]:


'|'.join(['foo', 'bar','ooki'])


# In[5]:


def fastq_gz_strings(filename):
    with gzip.open(filename, "rt") as handle:    
        for r in SeqIO.parse(handle, "fastq"):
            yield str(r.seq)
             
def gz_strings(filename):
    with gzip.open(filename, "rt") as handle:    
        for line in handle:
            if line[0] not in '@+#':  # skip fastq headers/quality
                yield line


# In[6]:


print('Hello World!')
test1(10)
print('Here come the SINES!')


#    search_sines("mouse SINEs.fasta",)
#    

good_lines = [
    'TGATTATCAGGTGAGAAATCACGATGGGAATTAAAAGCATTCTGAAGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAAGCAGAGGCAGACGGATTTCTGAATTCGAGGCCAGCCTGGTCTACAGAGTGAGTTCCAGGAC',
    'GAATCCTTGTTTTACAGCTGGATACGATGTAGGCTTACAGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGTGGATTTCTGAGTTCGAGGCCAGCCTGGTCTACAAAGTGAGTTCCAGGACAGCCAGG',
    'TTTTGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGCGGATTTCAGAGTTTGAGGCCAGCCTGGTCTACAAAGTGAGTTCCAGGACAGCTGGGCTACAGAGAAATCCTGACTTAAAAAAACAAAAACA',
    'TTTTGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGCGGATTTCAGAGTTTGAGGCCAGCCTGGTCTACAAAGTGAGTTCCAGGACAGCTGGGCTACAGAGAAAGCCTGACTTAAAAAAACAAAAACA',
    'GCAGGTAAGAACCATCAAAGCGACCCTATTAGGTAAATCCTGATAATATTCCATTTTAAAAATGGTGAAAGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGCGGATTTCTGAGTTCGAGGCCAGCCT',
    'CATAAGAAAGAGCTGTGCGGCCGGGCATGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGTGGATTTCTGAGTTCGAGGCCAGCCTGGTCTACAAAGTGAGTTTCAGGACAGCCAGGGCTATACAGAGAAACCCTGT',
]
search_sines("mouse SINEs.fasta",good_lines)


# In[ ]:


search_sines("mouse SINEs.fasta",fastq_gz_strings('old_lung_R2_001.fastq.gz'), 6, 100000)


# In[26]:


# TODO: DO NOT USE?  20% faster but different results
search_sines("mouse SINEs.fasta", gz_strings('old_lung_R2_001.fastq.gz'), 10000)


# # End of code, Saved results below

# In[ ]:


# Test results for 40, r1_old, 1000000
'0 3 / 5994000
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
34 5 / 5994000'

# Test results for 40, r2_old, 1000000
'0 9 / 6000000
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
34 54 / 6000000'

'time elapsed  89.94379512866338 wt-r1, 1000000
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
35 1 / 6000000'

'time elapsed  15.340233178933461 old-r2, 30000000
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
28 1 / 24000000'

# ==================
'time elapsed  19.913855942090354 1000000 young r2
0 8 / 6000000
1 190 / 6000000
2 623 / 6000000
3 1446 / 6000000
4 1821 / 6000000
returning with nlines =  1000000'

#=============
'time elapsed  20.19312702814738 1000000 old r2
0 9 / 6000000
1 110 / 6000000
2 419 / 6000000
3 973 / 6000000
4 1259 / 6000000
returning with nlines =  1000000'

'time elapsed 16.6222222050031 minutes
0 9 / 6000000
1 110 / 6000000
2 419 / 6000000
3 971 / 6000000
4 1258 / 6000000
returning with nlines = 1000000'

