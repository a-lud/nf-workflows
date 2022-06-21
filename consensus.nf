/*
Consensus pipeline
This pipeline has been written to have as few temporary files as possible.
Currently the pileup and variant calling are collapsed into a single
process, while all filtering processes are conducted in a separate single
process. With small data, the intermediate stages could be separated out
for more granular control, but with large files it'll become problematic with
large VCF files being present at each stage of the pipeline.
    * Genotypes called with mpileup before variants are called with call
    * A range of filtering options are available for the user to provide.
*/

// Import pipeline functions
include { bcftools_mpileup_call } from '../nf-modules/bcftools/1.15.1/bcftools_mpileup_call'
include { bcftools_consensus } from '../nf-modules/bcftools/1.15.1/bcftools_consensus'

// Sub-workflow
workflow CONSENSUS {
    main:

    // Get BAM file: [ bam.basename, [ bam ] ]
    Channel
        .fromFilePairs(
            [params.bamdir.path, params.bamdir.pattern].join('/'),
            size: params.bamdir.nfiles,
        )
        .ifEmpty { exit 1, "Can't find BAM files." }
        .set { ch_bam }

    // Parse CSV [ bam.basename, [reference] ]
    Channel
        .fromPath(
            params.sheet
        )
        .splitCsv()
        .map { row -> 
            tuple(row[0], file(row[1]))
        }
        .set { ch_csv }
    
    // Join BAMs with their reference file [ bam.id, [bam], [reference] ]
    ch_bam.join(ch_csv).set { ch_input }

    /*
    Optional arguments
        * Attempted to set some sensible defaults for consensus pipelines.
          Typically you don't want SNPs near INDELs, or the INDELs themselves.
          Variant normalisation is also needed to prevent MNPs from causing strife
          at the consensus stage. Finally, first ALT allele is used as the new
          sequence in the consensus.
    */
    def mpq =  params.containsKey('mapq') ?
        params.mapq :
        10
    
    def bsq = params.containsKey('baseq') ?
        params.baseq :
        10

    def ploidy = params.containsKey('ploidy') ?
        params.ploidy :
        2

    def mpileupOpt = params.containsKey('mpileup_opt') ?
        params.mpileup_opt :
        ''
    
    def callOpt = params.containsKey('call_opt') ?
        params.call_opt :
        ''
    
    def filterOpt = params.containsKey('filter_opt') ?
        params.filter_opt :
        '--SnpGap 10'

    def viewOpt = params.containsKey('view_opt') ?
        params.view_opt :
        '--exclude-types indels'

    def normOpt = params.containsKey('norm_opt') ?
        params.norm_opt :
        '-m +any'
    
    def sortOpt = params.containsKey('sort_opt') ?
        params.sort_opt :
        ''

    def consensusOpt = params.containsKey('consensus_opt') ?
        params.consensus_opt :
        '-H 1'

    // Outprefix is used as a parent directory here, so adjust the outpath
    def outdir = [params.outdir, params.out_prefix].join('/')

    // Genotype and Call
    bcftools_mpileup_call(
        ch_input,
        params.mapq,
        params.baseq,
        ploidy,
        mpileupOpt,
        callOpt,
        outdir
    )

    // Make consensus sequence using VCFs
    bcftools_consensus(
        bcftools_mpileup_call.out,
        filterOpt,
        viewOpt,
        normOpt,
        sortOpt,
        consensusOpt,
        outdir
    )
}