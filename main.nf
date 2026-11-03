include { BWAMEM2_INDEX } from './modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH } from './modules/nf-core/samtools/depth/main'

process COUNT_READS {
  tag "${meta.id}"
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
    tag "${meta.id}"
    label 'process_low'

    conda "${workflow.projectDir}/modules/nf-core/samtools/index/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

  input: tuple val(meta), path(bam)
  output:
  tuple val(meta), path("${bam.baseName}_coverage.tsv"), emit: coverage
  tuple val(meta), path("${bam.baseName}_tophit.tsv"), emit: tophit

  script:
  def segments = params.segments
  """
  samtools coverage ${bam} > ${bam.baseName}_coverage.tsv
  # Sort by coverage, then by mapped depth
  cat ${bam.baseName}_coverage.tsv | sort -k6,6nr -k7,7nr > ${bam.baseName}_sorted.tsv

  # 1=reference; 4=numreads; 6=coverage; 7=meandepth
  for segment in ${segments}; do
    grep "|\$segment" ${bam.baseName}_sorted.tsv \
        | awk -F'\t' -v seg="\$segment" 'OFS="\t" {{print "${bam.baseName}", seg, \$1, \$4, \$6, \$7}}' \
        | head -n1 \
        >> ${bam.baseName}_tophit.tsv
  done
  """
}

process COVERAGE_TO_WIDE {
  conda 'conda-forge::python=3.12 conda-forge:pandas=2.2'

  container "${workflow.containerEngine in ['singularity', 'apptainer']
    ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/8829bf123fc6b4b654fc8ab8b619d7f0e895993227a81150fa39aa590da91b1c/data'
    : 'community.wave.seqera.io/library/pandas:3.1.0--f6c858a02e385828'}"

  input:
    tuple path(top_coverage), path(readcount)

  output:
    path("results-coverage.tsv")

  script:
  """
  python ${workflow.projectDir}/bin/coverage_to_wide.py \\
    --top-coverage ${top_coverage} \\
    --readcount ${readcount} \\
    --merged results-coverage.tsv
  """
}

process PLOT_DEPTH {
  tag "${meta.id}"
  label 'process_low'

  conda 'conda-forge::r-base=4.4 conda-forge::r-ggplot2'

  input: tuple val(meta), path(depth), path(tophit)
  output: tuple val(meta), path("${depth.baseName}_depth.png")
  script:
  """
  Rscript ${workflow.projectDir}/bin/plot_depth.R ${depth} ${tophit} ${depth.baseName}_depth.png
  """
}

process FILTER_TOPHIT_DEPTH {

    input:
    tuple val(meta), path(depth), path(tophit)

    output:
    tuple val(meta), path("${meta.id}.top-hit-depth.tsv"), emit: tsv

    script:
    """
    awk -v sample="${meta.id}" '
        BEGIN {
            FS = OFS = "\\t"
        }

        # Read top-hit references (column 3)
        NR == FNR {
            refs[\$3] = 1
            next
        }

        # Retain only matching reference sequences
        \$1 in refs {
            print sample, \$0
        }
    ' $tophit $depth > ${meta.id}.top-hit-depth.tsv
    """
}

// process ALL_SAMPLE_DEPTH {
//     input:
//     path(depth_tsvs)

//     output:
//     path("all-combined-depth.tsv")

//     script:
//     """
//     printf 'Sample\tCHROM\tPOS\tDepth\n' > all-combined-depth.tsv

//     for depth_file in ${depth_tsvs}; do
//         sample=\$(basename "\$depth_file" .tsv)
//         awk -v sample="\$sample" 'BEGIN {OFS="\\t"} { print sample, \$0 }' "\$depth_file" >> all-combined-depth.tsv
//     done
//     """
// }

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
       return tuple([n.get(0), n.get(1)], [n.get(2), n.get(3)], [n.get(4), n.get(5)], params.keep_multimapped)
     }

    BWAMEM2_MEM(
      bwamem2mem_input | map{ n -> n.get(0)},
      bwamem2mem_input | map{ n -> n.get(1)},
      bwamem2mem_input | map{ n -> n.get(2)},
      bwamem2mem_input | map{ n -> n.get(3)}
    )

    // // Generate processed samplesheet
    // BWAMEM2_MEM.out.bam
    // .map { meta, bam -> "${meta.id},bwamem2_mem/${bam.fileName}" }
    // .collectFile(
    //   name: "${params.outdir}/processed_samplesheet.csv",
    //   seed: "sample,bam",
    //   newLine: true,
    // )

    // Calculate a single coverage value
    BWAMEM2_MEM.out.bam
    | SAMTOOLS_COVERAGE

    // Generate processed samplesheet with top-hit references
    bam_ch = BWAMEM2_MEM.out.bam
    .map { meta, bam -> tuple(meta.id, meta, bam) }

    tophit_ch = SAMTOOLS_COVERAGE.out.tophit
    .map { meta, tophit -> tuple(meta.id, tophit) }

    bam_tophit_ch = bam_ch
    .join(tophit_ch)
    .map {id, meta, bam, tophit ->
        // Extract reference IDs (column 3)
        def references = tophit.text.readLines()
          .collect { it.split('\t')[2] }
          .join(';')

        "${meta.id},bwamem2_mem/${bam.fileName},${references}"
    }

    bam_tophit_ch
    .collectFile(
        name: "${params.outdir}/processed_samplesheet.csv",
        seed: "sample,bam,references",
        newLine: true
    )

    // Get read depth across all positions of reference
    BWAMEM2_MEM.out.bam
    | SAMTOOLS_INDEX

    // Wide coverage table
    top_coverage_ch = SAMTOOLS_COVERAGE.out.tophit
    .map { _meta, tsv -> tsv }
    .collectFile(
      name: "${params.outdir}/merge/all-top-coverage.tsv"
    )
    readcount_ch = COUNT_READS.out.tsv
    .map { _meta, tsv -> tsv }
    .collectFile(
      name: "${params.outdir}/merge/all-readcount.tsv"
    )
    merge_inputs_ch = top_coverage_ch.combine(readcount_ch)
    COVERAGE_TO_WIDE(
      merge_inputs_ch
    )

    // Get depth across entire segment
    SAMTOOLS_DEPTH(
      BWAMEM2_MEM.out.bam | join(SAMTOOLS_INDEX.out.index),
      [[],[]] // All positions, so pass an empty positional value, this might get long if there are many ~200 references
    )

    left_depth_ch = SAMTOOLS_DEPTH.out.tsv
    .map { meta, tsv -> tuple(meta.id, meta, tsv)}
    right_tophit_ch = SAMTOOLS_COVERAGE.out.tophit
    .map { meta, tsv -> tuple(meta.id, tsv)}

    depth_tophit_ch = left_depth_ch
    .join(right_tophit_ch)
    .map { id, meta, depth, tophit -> tuple(meta, depth, tophit)}

    FILTER_TOPHIT_DEPTH(
      depth_tophit_ch
    )

    FILTER_TOPHIT_DEPTH.out.tsv
    .map { _meta, tsv -> tsv }
    .collectFile(
      name: "${params.outdir}/results-depth.tsv",
      seed: "Sample\tCHROM\tPOS\tDepth",
      newLine: true,
    )

    //ALL_SAMPLE_DEPTH(depth_tsv_ch)

    if (params.plot_depth){
      SAMTOOLS_DEPTH.out.tsv
      | join(SAMTOOLS_COVERAGE.out.tophit)
      | PLOT_DEPTH
    }
}
