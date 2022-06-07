/*
QC pipeline
    * FastQC
    * Kraken2 filter
    * Fastp quality + length trimming
    * MultiQC report
*/

// Import pipeline functions
include { fastqc } from '../nf-modules/fastqc/0.11.8/fastqc'
include { kraken2 } from '../nf-modules/kraken2/2.1.2/kraken2'
include { fastp_paired } from '../nf-modules/fastp/0.23.2/fastp_paired'
include { multiqc } from '../nf-modules/multiqc/1.12/multiqc'

// Sub-workflow
workflow QC {
    main:

    // Get sequence data
    Channel
        .fromFilePairs(
            [params.seqdir.path, params.seqdir.pattern].join('/'),
            size: params.seqdir.nfiles,
        )
        .ifEmpty { exit 1, "Can't find read files." }
        .set { ch_reads }

    // Handle optional arguments
    def bq_phred =  params.containsKey('bq_phred') ?
        params.bq_phred :
        15
    
    def n_base_limit =  params.containsKey('n_base_limit') ?
        params.n_base_limit :
        5
    
    def average_qual =  params.containsKey('average_qual') ?
        params.average_qual :
        0
    
    def length_required =  params.containsKey('length_required') ?
        params.length_required :
        15

    // FastQC
    fastqc(
        ch_reads,
        params.out_prefix,
        params.outdir
    )

    // Kraken filtering
    kraken2(
        ch_reads,
        params.out_prefix,
        params.krakendb,
        params.outdir
    )

    // Fastp
    fastp_paired(
        kraken2.out.unclassified,
        params.out_prefix,
        params.platform,
        bq_phred,
        n_base_limit,
        average_qual,
        length_required,
        params.outdir
    )

    // Collect all the files and pass to multiqc
    fastqc.out.zip.collect().set { ch_fqc }
    kraken2.out.report.collect().set { ch_kr2 }
    fastp_paired.out.json.collect().set { ch_fp }

    // Generate MultiQC report
    multiqc(
        ch_fqc,
        ch_kr2,
        ch_fp,
        params.out_prefix,
        params.outdir
    )
}