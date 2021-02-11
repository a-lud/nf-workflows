/*
HyPhy pipeline
    * Conduct selection analyses using HyPhy
*/

// Import utility functions
include {checkHyphyArgs;printHyphyArgs} from '../lib/utils'

// Check data
checked = checkHyphyArgs(params)
printHyphyArgs(checked, params.pipeline)

// Import pipeline modules
include { fel } from '../nf-modules/hyphy/2.5.25/fel'
include { fel } from '../nf-modules/hyphy/2.5.25/slac'
include { fel } from '../nf-modules/hyphy/2.5.25/fubar'
include { fel } from '../nf-modules/hyphy/2.5.25/meme'
include { fel } from '../nf-modules/hyphy/2.5.25/absrel'
include { fel } from '../nf-modules/hyphy/2.5.25/busted'
include { fel } from '../nf-modules/hyphy/2.5.25/relax'

// Sub-workflow
workflow HYPHY {
    main:
        // Data channel - Fasta files
        files_path = params.files_dir + '/' + params.files_ext
        Channel
            .fromPath(files_path)
            .ifEmpty { exit 1, "Can't import files at ${files_path}"}
            .collect()
            .toList()
            .set { ch_aln }
        
        // Data channel - Tree file
        Channel
            .fromPath(checked.tree)
            .ifEmpty { exit 1, "Can't import tree file ${params.tree}"}
            .set { ch_tree }
        
        // Combine alignment files + trees
        ch_aln
            .combine(ch_tree)
            .set { ch_inputs }

        if(checked.method.any { it == 'fel' }) {
            fel(ch_inputs, params.outdir, params.fel_optional)
        }

        if(checked.method.any { it == 'slac' }) {
            slac(ch_inputs, params.outdir, params.fel_optional)
        }

        if(checked.method.any { it == 'fubar' }) {
            fubar(ch_inputs, params.outdir, params.fel_optional)
        }

        if(checked.method.any { it == 'meme' }) {
            meme(ch_inputs, params.outdir, params.fel_optional)
        }

        if(checked.method.any { it == 'absrel' }) {
            absrel(ch_inputs, params.outdir, params.fel_optional)
        }

        if(checked.method.any { it == 'busted' }) {
            busted(ch_inputs, params.outdir, params.fel_optional)
        }

        if(checked.method.any { it == 'relax' }) {
            relax(ch_inputs, params.outdir, params.fel_optional)
        }
}