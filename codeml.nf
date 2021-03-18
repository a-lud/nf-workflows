/*
CodeML pipeline
    * Conduct selection analyses using CodeML from the PAML package
*/

// Import utility functions
include { codeml } from '../nf-modules/ete3/3.1.2/codeml'

// Check data
checked = checkCodemlArgs(params)
printCodemlArgs(checked, params.pipeline)

// Sub-workflow
workflow CODEML {
    // main:

    // Configure ete evol to work
    // File path = new File("${FASTDIR}/nf-conda_envs/ete")
    // Channel.value(path.isDirectory()).set { check }
    // setup_ete(workflow, check)
    // setup_ete.out.ifEmpty('exists').set { setup }

    // Run ete-evol
    // run_codeml(setup,
    //            seqs_tree_mark,
    //            params.outdir,
    //            params.models,
    //            params.tests,
    //            params.leaves,
    //            params.internals,
    //            params.codeml_param,
    //            workflow)
}