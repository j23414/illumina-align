include { FASTQC } from './modules/nf-core/fastqc/main' 
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { BWAMEM2_MEM } from '../modules/nf-core/bwamem2/mem/main'

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
}