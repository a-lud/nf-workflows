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

        ch_aln.combine(ch_tree).view()

        // fel(ch_aln, ch_tree, params.outdir, params.fel_optional)

        // Obtain MSA files


        // Obtain tree files


        // HyPhy analyses
        
}