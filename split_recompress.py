#!/usr/bin/env python3.6
import gzip
import sys
#import zstandard as zstd
import itertools
import time
import os

# https://stackoverflow.com/a/24527424/1413499
def chunks(iterable, size):
    iterator = iter(iterable)
    for first in iterator:
        yield itertools.chain([first], itertools.islice(iterator, size - 1))

def split_recompress(basename, skip=[]):
    """
    skip - indexes of parts to skip (still need to decompress but not compress)
    """
    in_fname = f'{basename}.fastq.gz'
    with gzip.open(in_fname, 'rb') as in_fastq:
        for i, chunk in enumerate(chunks(in_fastq, 100_000_000)):
            out_fname = f'{basename}.part{i}e8.fastq.gz'
            if i in skip:
                print('skipping', out_fname)
                max(chunk)  # consume iterator
                continue
            print('writing', out_fname)
            t0 = time.time()
            with gzip.open(out_fname, 'wb', compresslevel=6) as out_fastq_gz:
            #with open(out_fname, 'wb') as out_fastq_gz:
                # a line-by-line python loop is subotimal, but ~60% CPU is spent in
                # ZSTD_compressBlock_lazy2_extDict so good enough for a one-off.
                for line in chunk:
                    out_fastq_gz.write(line)
            t1 = time.time()
            print(t1 - t0, 'sec.', 'Last line:', line)
    os.remove(in_fname)

[in_filename_base] = sys.argv[1:]
split_recompress(in_filename_base)
#split_recompress('Old-lung/old_lung_R1_001')
#split_recompress('Old-lung/old_lung_R2_001')
#split_recompress('Young-lung/wt-lung_R1_001')
#split_recompress('Young-lung/wt-lung_R2_001')
