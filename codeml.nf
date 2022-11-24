/*
CodeML pipeline

This pipeline runs PAML selection tests via ETE3 evol. I've implemented
this version to also run drop-out analyses if Branch-Site models are run.
See here for more information: https://doi.org/10.1101/2021.10.26.465984 

Pipeline overview:
    * Run all user models
    * Build drop-out tree file if Branch-Site models are specified
    * Run drop-out branch site models
    * Build summary tables
    * For BS models -> compare to drop out?
*/

// Import pipeline functions
include { remove_foreground } from '../nf-modules/general/remove_foreground'
include { codeml } from '../nf-modules/ete3/3.1.2/codeml'
include { codeml as codeml_dropout } from '../nf-modules/ete3/3.1.2/codeml'
include { etetools } from '../nf-modules/general/etetools'
include { etetools as etetools_dropout } from '../nf-modules/general/etetools'
include { compareLRT } from '../nf-modules/general/compareLRT'

// Sub-workflow
workflow CODEML {
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
    def outdir = params.outdir + '/' + params.out_prefix
    def dropout = params.containsKey('dropout') ?: false
    
    // Combine MSA files with tree file
    ch_msa.combine(ch_tree).set { ch_msa_tree }

    // Outdir for general run
    outdir_codeml = outdir + '/codeml'

    // Run codeml
    codeml(
        ch_msa_tree,
        params.models,
        outdir_codeml
    )

    // Build summary tables
    codeml.out.cml.collect().set { ch_codeml }
    etetools(ch_codeml, 'standard', outdir_codeml)

    // Run dropout analysis
    if (dropout) {
        if (params.models.tokenize(" ").containsAll(["bsA", "bsA1"])) {
            // Remove foreground branch and connected tips from tree
            remove_foreground(ch_tree, outdir)

            // Combine MSA files with foreground removed tree file
            ch_msa.combine(remove_foreground.out.no_fg).set { ch_msa_tree_no_fg }

            // Run ONLY site models
            outdir_dropout = outdir + '/codeml_dropout'
            
            // Dropout codeml run
            codeml_dropout(
                ch_msa_tree_no_fg,
                'M1 M2',
                outdir_dropout
            )

            codeml_dropout.out.cml.collect().set { ch_dropout }

            // Build summary tables
            etetools_dropout(ch_dropout, 'dropout', outdir_dropout)

            // compare drop out site models to branch-site models
            etetools.out.summary.set { ch_etetool_summary }
            etetools_dropout.out.summary.set { ch_dropout_summary }

            compareLRT(ch_etetool_summary, ch_dropout_summary, outdir_dropout)
        }
    }

}