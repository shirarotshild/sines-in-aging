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

# http://biopython.org/
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from gene_lib import log

def fastq_gz_strings(filename):
    with gzip.open(filename, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            yield str(record.seq)

def gz_strings(filename):
    with gzip.open(filename, "rt") as text:    
        for i, line in enumerate(text):
            # records of 4 lines: @header, dna, +, quality
            if i % 4 == 1:
                yield line.rstrip()

def fastq_zst_strings(filename):
    for rec in fastq_zst_records(filename):
        yield str(rec.seq)

def fastq_zst_records(filename):
    # https://github.com/indygreg/python-zstandard - pip3 install zstandard
    import zstandard as zstd
    log(f"Reading {filename}...") 
    with open(filename, 'rb') as fastq_zst_handle:
        fastq_handle = zstd.ZstdDecompressor().stream_reader(fastq_zst_handle)
        # wrapper adds support for .readline(), for line in ...
        fastq_text = io.TextIOWrapper(fastq_handle, encoding='ascii')
        for record in SeqIO.parse(fastq_text, "fastq"):
            yield record

def zst_strings(filename):
    # https://github.com/indygreg/python-zstandard - pip3 install zstandard
    import zstandard as zstd
    with open(filename, 'rb') as fastq_zst_handle:
        fastq_handle = zstd.ZstdDecompressor().stream_reader(fastq_zst_handle)
        # wrapper adds support for .readline(), for line in ...
        fastq_text = io.TextIOWrapper(fastq_handle, encoding='ascii')
        for i, line in enumerate(fastq_text):
            # records of 4 lines: @header, dna, +, quality
            if i % 4 == 1:
                yield line.rstrip()
