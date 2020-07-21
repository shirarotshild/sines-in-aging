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

def split_recompress(in_fname, out_basename, skip=[]):
    """
    skip - indexes of parts to skip (still need to decompress but not compress)
    """
    with gzip.open(in_fname, 'rb') as in_fastq:
        for i, chunk in enumerate(chunks(in_fastq, 100_000_000)):
            out_fname = f'{out_basename}.part{i}e8.fastq.gz'
            if i in skip:
                log('Skipping', out_fname)
                max(chunk)  # consume iterator
                continue

            log('Writing', out_fname)
            t0 = time.time()
            with gzip.open(out_fname + '.tmp', 'wb', compresslevel=2) as out_fastq_gz:
            #with open(out_fname, 'wb') as out_fastq_gz:
                out_fastq_gz.writelines(chunk)
            t1 = time.time()
            os.rename(out_fname + '.tmp', out_fname)
            log(t1 - t0, 'sec.')

#split_recompress('Old-lung/old_lung_R1_001.fastq.gz', 'Old-lung/old_lung_R1_001')

if __name__ == '__main__':
    [in_fname, out_basename] = sys.argv[1:]
    split_recompress(in_fname, out_basename)
    log('REMOVING', in_fname)
    os.remove(in_fname)  # may fail for /dev/stdin, /dev/fd/... etc. but that's OK
