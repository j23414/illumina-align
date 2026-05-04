include { FASTQC } from './modules/nf-core/fastqc/main' 
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { BWAMEM2_INDEX } from './modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'

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

    // Run QC
    FASTQC(reads_ch)
    FASTQC.out.html

    multiqc_input = FASTQC.out.html.map{meta, files -> files}
    | mix(FASTQC.out.zip.map{meta, files -> files})
    | flatten
    | collect
    | map {
      n ->
      def meta = [id: 'all']
      return tuple(meta, n, [], [], [], [])
    }

    MULTIQC(multiqc_input)
    MULTIQC.out.report.view()

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
    }

}
