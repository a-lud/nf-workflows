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

        // Convert to nucleotide
        if(params.pep2nuc) {

            // Tuple [ id, aln, nuc ]
            runMSA.out.alignments
                .map { ids, files ->

                    def lst = []
                    files.each { f ->
                        // Convert Unix path to string
                        String p = f
                        n = p.substring(p.lastIndexOf('/') + 1)
                        n = n.substring(0, n.lastIndexOf('_msa.fasta'))

                        // Find sample ID
                        i = ids.find { it == n }
                        lst.add([ i, f ])
                    }
                    return lst
                }
                .flatMap { return it}
                .set { ch_alignment }
            
            // Import nucleotide files
            nucleotide_files = params.nucleotide_dir + '/' + params.nucleotide_ext
            Channel
                .fromFilePairs(nucleotide_files, size: 1)
                .ifEmpty { exit 1, "No complementary nucleotide fasta files at ${params.nucleotide_dir}"}
                .set { ch_nucleotide }

            // Join into tuple
            ch_alignment.join(ch_nucleotide, by: [0]).set { ch_input_p2n }

            // Create three channels: IDs, aln, nuc
            ch_input_p2n.map { id, aln, nuc -> return id }.collect().map { return it.join(' ') }.set { ids }
            ch_input_p2n.map { id, aln, nuc -> return aln }.collect().set { aln }
            ch_input_p2n.map { id, aln, nuc -> return nuc }.collect().set { nuc }

            ids.view()
            aln.view()
            nuc.view()
        }

}