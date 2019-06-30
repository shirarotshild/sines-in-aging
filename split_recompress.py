#!/usr/bin/env python3
import gzip
import zstandard as zstd
import itertools
import time
import os

# https://stackoverflow.com/a/24527424/1413499
def chunks(iterable, size):
    iterator = iter(iterable)
    for first in iterator:
        yield itertools.chain([first], itertools.islice(iterator, size - 1))

def split_recompress(basename, skip=0):
    in_fname = f'{basename}.fastq.gz'
    with gzip.open(in_fname, 'rb') as in_fastq:
        for i, chunk in enumerate(chunks(in_fastq, 100_000_000)):
            out_fname = f'{basename}.part{i}e8.fastq.zst'
            if i < skip:
                print('skipping', out_fname)
                max(chunk)  # consume iterator
                continue
            print('writing', out_fname)
            t0 = time.time()
            with open(out_fname, 'wb') as out_fastq_zstd:
                cctx = zstd.ZstdCompressor(level=12, threads=4)
                with cctx.stream_writer(out_fastq_zstd) as out_fastq:
                    # a line-by-line python loop is subotimal, but ~60% CPU is spent in
                    # ZSTD_compressBlock_lazy2_extDict so good enough for a one-off.
                    for line in chunk:
                        out_fastq.write(line)
            t1 = time.time()
            print(t1 - t0, 'sec.', 'Last line:', line)
    #os.remove(in_fname)

split_recompress('old_lung_R1_001', skip=28)
# split_recompress('old_lung_R2_001', skip=6)
# split_recompress('wt-lung_R1_001.fastq.gz')
# split_recompress('wt-lung_R2_001.fastq.gz')
