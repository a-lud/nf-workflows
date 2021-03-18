/*
CodeML pipeline
    * Conduct selection analyses using CodeML from the PAML package
*/

// Import utility functions
include {checkCodemlArgs;printCodemlArgs} from '../lib/utils'

// Import utility functions
include { codeml } from '../nf-modules/ete3/3.1.2/codeml'

// Check data
checked = checkCodemlArgs(params)
printCodemlArgs(checked, params.pipeline)

// Sub-workflow
workflow CODEML {
    main:

    files_path = params.files_dir + '/' + params.files_ext
    Channel
        .fromFilePairs(files_path, size: 1)
        .ifEmpty { exit 1, "Can't import files at ${files_path}" }
        .set { ch_aln }
    
    Channel
        .fromPath(params.trees.tokenize(','))
        .ifEmpty { exit 1, "Can't import tree files" }
        .set { ch_tree }
    
    ch_aln
        .combine( ch_tree )
        .set { ch_seq_tree }

    // Run codeml
    codeml(ch_seq_tree,
           params.models,
           params.tests,
           params.codeml_optional)
}