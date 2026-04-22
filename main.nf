include { FASTQC } from './modules/nf-core/fastqc/main' 
include { MULTIQC } from './modules/nf-core/multiqc/main'     

workflow {
    main:
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
      | view {files -> "Read files : $files "}

    // reads_ch
    // | FASTQC
    // | MULTIQC
    // | view
}