# illumina-align

A modular workflow for aligning Illumina sequencing data and summary statistics prior to downstream analysis.

## Usage

```
nextflow run j23414/illumina-align \
  --reads [path/*_{R1,R2}_001.fastq.gz \
  --reference [path/reference.fasta] \
  --outdir "align-results" \
  -profile stjude
```

## Optional post-processing merging

```
# (option A) Just get table of coverage
cat align-results/samtools_coverage/*tophit* > all-top-coverage.tsv
python illumina-align/bin/coverage_to_wide.py \
  --top-coverage all-top-coverage.tsv \
  --merged results-coverage.tsv

# (option B) Combine with table of initial read counts
cat align-results/countreads/*_readcount.tsv > all-readcount.tsv
python illumina-align/bin/coverage_to_wide.py \
  --top-coverage all-top-coverage.tsv \
  --readcount all-readcount.tsv \
  --merged results-coverage.tsv
```
