/*
HyPhy pipeline
    * Conduct selection analyses using HyPhy
*/

// Import utility functions
include {checkHyphyArgs;printHyphyArgs} from '../lib/utils'

// Import pipeline modules
include { fel } from '../nf-modules/hyphy/2.5.25/fel'

// Check data
checked = checkHyphyArgs(params)
printHyphyArgs(checked, params.pipeline)

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

        // Data channel - ID
        // ch_files
        //     .map { id, file ->
        //         return file
        //     }
        //     .collect()
        //     .set { ch_aln }
        
        // Combine alignment files + trees
        ch_aln
            .combine(ch_tree)
            .set { ch_inputs }

        if(checked.method.any { it == 'fel' }) {
            fel(ch_inputs, params.outdir, params.fel_optional)
        } else if(checked.method.any { it == 'absrel' } ) {
            println('aBSREL')
        } else if(checked.method.any { it == 'meme' } ) {
            println('MEME')
        }

        // 

        // Obtain MSA files


        // Obtain tree files


        // HyPhy analyses
        
}