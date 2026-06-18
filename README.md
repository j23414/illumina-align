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

