/*
HyPhy development pipeline pipeline

This workflow is designed to run 'HyPhy analyses' pipelines, which are custom written pipelines
that extend the default HyPhy methods.
    * 1. Get codon MSAs
    * 2. Match MSAs with tree
    * 3. Run BUSTED-PH for each gene
*/


// Import pipeline modules
include { busted_ph } from '../nf-modules/hyphy_analyses/16bc859/busted_ph'

// Sub-workflow
workflow HYPHY_ANALYSES {
    main:
        // Get MSA fasta files
    Channel
        .fromFilePairs(
            [params.msa.path, params.msa.pattern].join('/'),
            size: params.msa.nfiles,
        )
        .ifEmpty { exit 1, "Can't find MSA files." }
        .set { ch_msa }
    
    // Get tree file
    Channel
        .fromPath(
            params.tree
        )
        .ifEmpty { exit 1, "Can't import tree file" }
        .set { ch_tree }

    // Define some variables
    def libpath = params.exedir + '/res'

    def outdir = params.outdir + '/' + params.out_prefix
    def outhyphy = outdir + '/hyphy'
    def outlog = outdir + '/logs-hyphy'

    // Combine MSA files with tree file
    ch_msa.combine(ch_tree).set { ch_msa_tree }

    // Submit each MSA as separate job
    busted_ph(
        ch_msa_tree,
        params.exedir,
        libpath,
        params.batchFile,
        params.testLabel,
        outhyphy,
        outlog
    )
}