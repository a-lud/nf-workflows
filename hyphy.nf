/*
HyPhy pipeline

Implementation of HyPhy standard analyses.
*/

// Import pipeline modules
// include { fel } from '../nf-modules/hyphy/2.5.25/fel'
// include { slac } from '../nf-modules/hyphy/2.5.25/slac'
// include { fubar } from '../nf-modules/hyphy/2.5.25/fubar'
// include { meme } from '../nf-modules/hyphy/2.5.25/meme'
// include { absrel } from '../nf-modules/hyphy/2.5.25/absrel'
// include { busted } from '../nf-modules/hyphy/2.5.25/busted'
include { relax } from '../nf-modules/hyphy/2.5.42/relax'

// Sub-workflow
workflow HYPHY {
    main:
    // Data channel - Fasta files
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
    def outhyphy = outdir + '/hyphy'
    def outlog = outdir + '/logs-hyphy'

    // Combine MSA files with tree file
    ch_msa.combine(ch_tree).set { ch_msa_tree }

    // Submit each MSA as separate job
    switch(params.analysis) {
        case "RELAX":
            relax(
                ch_msa_tree,
                params.testLabel,
                outhyphy,
                outlog
            )
            break;
        default:
            println("Shiet...")
    }
}