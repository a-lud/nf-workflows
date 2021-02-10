/*
HyPhy pipeline
    * Conduct selection analyses using HyPhy
*/

// Import utility functions
include {checkHyphyArgs;printHyphyArgs} from '../lib/utils'
checked_args = checkHyphyArgs(params)
printHyphyArgs(checked_args)

// Import pipeline modules


// Sub-workflow
workflow HYPHY {
    main:

    println('HELLO')

        // Obtain MSA files


        // Obtain tree files


        // HyPhy analyses
        
}