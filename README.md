TODO document all stages

# r1r2merge.py

## Testing

```sh
./r1r2merge.py merge-test-R1.fastq merge-test-R2.fastq merge-test-merged.fastq.gz
gunzip merge-test-merged.fastq.gz
diff -s merge-test-merged.fastq merge-test-expected.fastq
```
