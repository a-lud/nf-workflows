/*
QC pipeline
    * FastQC
    * Kraken2 filter
    * Fastp quality + length trimming
    * MultiQC report
*/

// Import pipeline functions
include { fastqc } from '../nf-modules/fastqc/0.11.8/fastqc'
include { bbduk } from '../nf-modules/bbmap/39.01/bbduk'
include { kraken2 } from '../nf-modules/kraken2/2.1.2/kraken2'
include { fastp_paired } from '../nf-modules/fastp/0.23.2/fastp_paired'
include { multiqc } from '../nf-modules/multiqc/1.14/multiqc'

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

    def mink = params.containsKey('mink') ?
        params.mink : 11
    
    def rcomp = params.containsKey('rcomp') ?
        params.rcomp : 't'

    def ktrim = params.containsKey('ktrim') ?
        params.ktrim : 'r'

    def hdist = params.containsKey('hdist') ?
        params.hdist : 0
    
    def minlength = params.containsKey('minlength') ?
        params.minlength : 25

    // Handle optional arguments
    def bq_phred = params.containsKey('bq_phred') ?
        params.bq_phred :
        15
    
    def n_base_limit = params.containsKey('n_base_limit') ?
        params.n_base_limit :
        5
    
    def average_qual = params.containsKey('average_qual') ?
        params.average_qual :
        0
    
    def length_required = params.containsKey('length_required') ?
        params.length_required :
        15

    // out_prefix is used as a parent directory. Adjust the outdir variable
    def outdir = [params.outdir, params.out_prefix].join('/')

    // Introductory text for MQC report
    def intro = 'intro_text: Nextflow QC pipeline. Aggregation of FastQC, BBduk, Kraken2 and Fastp results'

    // FastQC
    fastqc(
        ch_reads,
        outdir + '/fastqc'
    )

    // BBduk - PCR primer in some libraries that needs removing
    bbduk(
        ch_reads,
        params.ref,
        params.k,
        mink,
        rcomp,
        ktrim,
        hdist,
        minlength,
        outdir + '/bbduk'
    )

    // Kraken filtering
    kraken2(
        bbduk.out.fq,
        params.krakendb,
        outdir + '/kraken2'
    )

    // Fastp
    fastp_paired(
        kraken2.out.unclassified,
        params.platform,
        bq_phred,
        n_base_limit,
        average_qual,
        length_required,
        outdir + '/fastp'
    )

    // Collect all the files and pass to multiqc
    fastqc.out.zip.collect().set { ch_fqc }
    bbduk.out.log.collect().set { ch_bbdk }
    kraken2.out.report.collect().set { ch_kr2 }
    fastp_paired.out.json.collect().set { ch_fp }

    // Concat all the channels so a single input can be passed to multiqc
    ch_fqc.concat(ch_kr2, ch_fp, ch_bbdk).collect().set { ch_multiqc }

    // Generate MultiQC report
    multiqc(
        ch_multiqc,
        projectDir + '/conf/multiqc-configs/mqc-qc.yml',
        intro,
        outdir
    )
}