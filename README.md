TODO document all stages

# Install tools

On ubuntu linux, try running `./INSTALL/ALL.sh`.

Note it installs python libraries for a _single_ python version.  All our *.py scripts should use that version.

See INSTALL/README.md for details on install script.

# Dowload input from AWS S3 storage

Choose an organ whose samples you want to process:

```
aws s3 ls s3://endogene/mice_wgs/Aging/
aws s3 ls s3://endogene/mice_wgs/Aging/Old-liver/
mkdir Old-liver
aws s3 cp s3://endogene/mice_wgs/Aging/Old-liver/old_liver_R1_001.fastq.gz Old-liver/
```

Now repeat for `_R2_` and for both Old and Young samples — but if you're short on disk space, you'll want to next steps for each one so you can delete inputs...

## streaming

The command supports `-` as either source (meaning stdin) or destination (meaning stdout):
```
aws s3 cp s3://endogene/mice_wgs/Aging/Old-liver/old_liver_R1_001.fastq.gz - | gunzip --stdout | head --lines=20
```

# breaking into chunks

for example in p53-liver
```
~/sines-in-aging/split_recompress.py p53-liver_R1_001.fastq.gz p53-liver_R1_001
~/sines-in-aging/split_recompress.py p52-liver_R2_001.fastq.gz p53-liver_R2_001
...
```
This cuts input into chunks of 100 000 000 (1e8) lines, writing files named:
```
p53-liver_R1_001.part0e8.fastq.gz
p53-liver_R1_001.part1e8.fastq.gz
...
p53-liver_R1_001.part27e8.fastq.gz
...
```
(last chunk will be smaller.)

This allows following steps to be run from the middle and/or parallelized.

# r1r2merge.py

```
parallel -j+0 --eta --joblog=old-liver-merge.joblog ~/sines-in-aging/r1r2merge.py 'old_liver_R1_001.part{}e8.fastq.gz' 'old_liver_R2_001.part{}e8.fastq.gz' 'old_liver_merged.part{}e8.fastq.gz' 
```

## Testing

```sh
./r1r2merge.py merge-test-R1.fastq merge-test-R2.fastq merge-test-merged.fastq.gz
gunzip merge-test-merged.fastq.gz
diff -s merge-test-merged.fastq merge-test-expected.fastq
```

## Process

- Currently, run_part_1.py is the high level script. This is messy and should be refactored in the future
- To generate potential sines run mode = 1.
- To generate barcodes run mode = 3
- cat $(ls -v old_liver_merged*/old_liver_merged.part*e8_sineBarcode.fastq.gz) > old_liver_merged_sineBarcode.fastq.gz

