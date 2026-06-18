include { BWAMEM2_INDEX } from './modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'

process COUNT_READS {
  label 'process_low'

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("*.tsv"), emit: tsv, optional: true

  script:
  """
  export R1_count=\$(( \$(gzip -dc ${reads[0]} | wc -l) / 4 ))
  export R2_count=\$(( \$(gzip -dc ${reads[1]} | wc -l) / 4 ))

  printf "%s\\t%s\\t%s\\n" \
    "${meta.id}" \
    "\$R1_count" \
    "\$R2_count" \
    > ${meta.id}_readcount.tsv
  """
}

process SAMTOOLS_COVERAGE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

  input: tuple val(meta), path(bam)
  output: tuple path("${bam.baseName}_coverage.tsv"), path("${bam.baseName}_tophit.tsv")
  script:
  """
  samtools coverage ${bam} > ${bam.baseName}_coverage.tsv
  cat ${bam.baseName}_coverage.tsv | sort -k6,6nr > ${bam.baseName}_sorted.tsv

  # 1=reference; 4=numreads; 6=coverage; 7=meandepth
  for segment in HA NA PB2 PB1 PA NP MP NS; do
    grep "|\$segment" ${bam.baseName}_sorted.tsv \
        | awk -F'\t' -v seg="\$segment" 'OFS="\t" {{print "${bam.baseName}", seg, \$1, \$4, \$6, \$7}}' \
        | head -n1 \
        >> ${bam.baseName}_tophit.tsv
  done
  """
}

workflow {
    main:
    // Load illumina reads
    if (params.samplesheet) {
      reads_ch = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            tuple([id: row.sample], [file(row.R1), file(row.R2)])
          }
    } else if (params.reads) {
        reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
          | map { name, reads -> tuple([id: name], reads) }
    } else {
        error "Please specify either --samplesheet samplesheet.csv or --reads 'data/*_{R1,R2}.fastq.gz'"
    }
    reads_ch | COUNT_READS

    // Load reference
    if (params.reference){
      reference_ch = channel.fromPath(params.reference, checkIfExists:true)
      | map { n -> tuple(n.baseName, n)
        }
    } else {
      error "Please specify --reference reference.fasta to align against"
    }

    BWAMEM2_INDEX(reference_ch)
    bwamem2mem_input = reads_ch
    | combine(BWAMEM2_INDEX.out.index)
    | combine(reference_ch)
    | map {
       n ->
       return tuple([n.get(0), n.get(1)], [n.get(2), n.get(3)], [n.get(4), n.get(5)], true )
     }

    BWAMEM2_MEM(
      bwamem2mem_input | map{ n -> n.get(0)},
      bwamem2mem_input | map{ n -> n.get(1)},
      bwamem2mem_input | map{ n -> n.get(2)},
      bwamem2mem_input | map{ n -> n.get(3)}
    )

    BWAMEM2_MEM.out.bam
    | SAMTOOLS_COVERAGE

}
