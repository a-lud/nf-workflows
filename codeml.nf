/*
CodeML pipeline
    * Selection testing using ETE3 evol's CodeML implementation
*/

// Import pipeline functions
include { codeml } from '../nf-modules/ete3/3.1.2/codeml'

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
    
    // Combine MSA files with tree file
    ch_msa.combine(ch_tree).set { ch_msa_tree }

    // Run codeml
    codeml(ch_msa_tree,
           params.models,
           params.outdir)
}