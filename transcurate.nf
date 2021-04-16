/*
Curate Transcripts pipeline
    * Optionally cluster redundant transcripts using CD-HIT
    * Identify CDS regions using TransDecoder
    * Optionally extract complete ORFs
*/

// Import utility functions
include {checkTransCurateArgs; printTransCurateArgs} from '../lib/utils'
checked = checkTransCurateArgs(params)
printTransCurateArgs(checked, params.pipeline)

// Import pipeline modules
include { downloadDB } from '../nf-modules/general/transcurate_downloadDB'
include { cdhit } from '../nf-modules/cdhit/4.8.1/cdhit'
include { transdecoderLO } from '../nf-modules/transdecoder/5.5.0/transdecoderLO'
include { transdecoderP } from '../nf-modules/transdecoder/5.5.0/transdecoderP'
include { blast } from '../nf-modules/blast/2.11.0/blast'
include { hmmer } from '../nf-modules/hmmer/3.3.2/hmmer'
include { gffread } from '../nf-modules/gffread/0.12.1/gffread'

// Sub-workflow
workflow TRANSCURATE {
    main:
        // Data channel - Fasta files
        files_path = params.files_dir + '/' + params.files_ext
        Channel
            .fromFilePairs(files_path, size: 1)
            .ifEmpty { exit 1, "Can't import files at ${files_path}"}
            .set { ch_transcripts }

        // Database
        if (params.database_dir) {
            Channel
                .fromPath(checked.database_dir, type: 'dir')
                .ifEmpty { exit 1, "Database directory doesn't exist: ${params.database_dir}"}
                .set { ch_database_dir }
        } else {
            downloadDB(params.outdir)
            downloadDB.out.database_files.set { ch_database_dir }
        }

        // Remove transcript redundancy
        cdhit(ch_transcripts, checked.cdhit_pid, params.outdir)
        
        // TransDecoder step 1
        transdecoderLO(cdhit.out.fasta, params.outdir)

        // Homology search
        blast(transdecoderLO.out.longstORF, ch_database_dir, params.outdir)
        hmmer(transdecoderLO.out.longstORF, ch_database_dir, params.outdir)

        // Join channels
        transdecoderLO.out.all.join(blast.out.join(hmmer.out, by: [0]), by: [0]).set { ch_all }

        // Transdecoder step 2
        transdecoderP(ch_all, params.outdir)

        // Complete ORF sequences
        if (checked.completeORFs) {
            gffread(transdecoderP.out.complete, params.outdir)
        }
}