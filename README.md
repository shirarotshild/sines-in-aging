TODO document all stages

# breaking into parts

for example in p53-liver
```
~/sines-in-aging/split_recompress.py p53-liver_R1_001
~/sines-in-aging/split_recompress.py p52-liver_R2_001
parallel -j+0 --eta --joblog=old-liver-merge.joblog ~/sines-in-aging/r1r2merge.py 'old_liver_R1_001.part{}e8.fastq.gz' 'old_liver_R2_001.part{}e8.fastq.gz' 'old_liver_merged.part{}e8.fastq.gz' 
...
```

# r1r2merge.py
 
## Testing

`sh
./r1r2merge.py merge-test-R1.fastq merge-test-R2.fastq merge-test-merged.fastq.gz
gunzip merge-test-merged.fastq.gz
diff -s merge-test-merged.fastq merge-test-expected.fastq
```
## Process

- Currently, run_part_1.py is the high level script. This is messy and should be refactored in the future
- To generate potential sines run mode = 1.
- To generate barcodes run mode = 3
- Rg
- cat $(ls -v old_liver_merged*/old_liver_merged.part*e8_sineBarcode.fastq.gz) > old_liver_merged_sineBarcode.fastq.gz

