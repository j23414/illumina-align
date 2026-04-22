include { FASTQC } from './modules/nf-core/fastqc/main' 
include { MULTIQC } from './modules/nf-core/multiqc/main'     

workflow {
    main:
    // Illumina Reads
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
      | take(3)
      | map { n ->
            def meta = [id: n.get(0) ]
            return tuple(meta, n.get(1))
        }
      | view

    // Run QC
    FASTQC(reads_ch)
    FASTQC.out.html.view()
}