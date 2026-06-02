include { BWAMEM2_INDEX } from './modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_VIEW } from './modules/nf-core/samtools/view/main'

process GET_REFERENCE_NAMES {
  input: path(reference)
  output: path("${reference.baseName}_names.txt")
  script:
  """
  grep ">" ${reference} | sed 's/>//g' | sed 's/ //g' > ${reference.baseName}_names.txt
  """
}

process MY_SAMTOOLS_VIEW {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input: tuple val(reference_name), path(bam)
    output: tuple val("${reference_name.tokenize('|')[0]}"), path("${bam.baseName}_*.bam")

    script:
    def refname = reference_name.tokenize('|')[0]

    """
    samtools index ${bam}
    samtools view -h ${bam} | grep -e "${refname}" -e "^@PG" | samtools view -bS > ${bam.baseName}_${refname}.bam
    """
}

process VIRAL_CONSENSUS {
    maxForks 5
    memory 513.GB
    label 'process_high'

    input: tuple val(refname), path(bam), path(ref)
    output: path("${bam.baseName}.fasta")
    script:
    """
    grep -A1 ${refname} ${ref} > ${refname}.fasta

    viral_consensus \
    -i ${bam} \
    -r ${refname}.fasta \
    -o ${bam.baseName}.fasta
    """
}

workflow {
    main:
    // Illumina Reads
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
//      | take(3)
      | map { n ->
            def meta = [id: n.get(0) ]
            return tuple(meta, n.get(1))
        }
//      | view

    // Load reference
    if(params.reference){
      reference_ch = channel.fromPath(params.reference, checkIfExists:true)
      | map {
        n ->
        return tuple(n.baseName, n)
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

      // Split reference by segment name
      reference_names_ch = reference_ch
      | map { n -> n.get(1) }
      | GET_REFERENCE_NAMES
      | map { n -> n.readLines()}
      | flatten

      reference_names_ch
      | combine(BWAMEM2_MEM.output.bam | map{ n -> n.get(1)})
      | MY_SAMTOOLS_VIEW
      | combine(reference_ch | map{ n -> n.get(1)})
      | VIRAL_CONSENSUS
      | view
    }
}
