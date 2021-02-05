/*
MSA pipeline
    * Align multi-fasta files using MAFFT,TCOFFE,MUSCLE,CLUSTAL
    * Convert peptide alignments to nucleotide
    * Clean alignments using GBlocks
*/

// Import utility functions
include {checkMsaArgs} from '../lib/utils'
checkMsaArgs(params)

// Import pipeline modules
include {runMSA} from '../nf-modules/general/runMSA'

workflow MSA {    
    main:

        // Data channel - Fasta files
        files_path = params.files_dir + '/' + params.files_ext
        Channel
            .fromFilePairs(files_path, size: 1)
            .ifEmpty { exit 1, "Can't import files at ${files_path}"}
            .set { ch_files }
        
        // Data channel - ID and files
        ch_files
            .map { id, file ->
                return id
            }
            .collect()
            .set { ch_ids }
        
        ch_files
            .map { id, file ->
                return file
            }
            .collect()
            .set { ch_files}

        // Align sequences
        runMSA(ch_files, ch_ids, params.outdir,
               params.aligner, params.aligner_args)

}