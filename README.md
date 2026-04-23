# unnamed

rename this later :P

## Usage

```
nextflow run j23414/unnamed -r main \
  --reads [path/*_{R1,R2}_001.fastq.gz \
  [--reference [path/reference.fasta] \
  --outdir "results" \
  -resume \
  -profile stjude
```

* If only the reads are provided, run QC (fastqc, multiqc). 
* If reference is provided, run BWAMEM2...
