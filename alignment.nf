/*
Alignment pipeline
    * Parse sample-sheet of sample-reference pairs
    * Align data using either BWA2/Minimap2
    * Sort BAM file using Sambamba
    * Deduplicate
    * Statistics
        * Flagstat
        * Mosdepth
*/

// Import pipeline functions
include { faidx } from '../nf-modules/samtools/1.15/faidx'
include { minimap2_sr } from '../nf-modules/minimap2/2.24/minimap2_sr'
include { bwa_mem2_index } from '../nf-modules/bwa-mem2/2.2.1/bwa-mem2-index'
include { bwa_mem2_mem } from '../nf-modules/bwa-mem2/2.2.1/bwa-mem2-mem'
include { markdup } from '../nf-modules/sambamba/0.8.2/markdup'
include { mosdepth } from '../nf-modules/mosdepth/0.3.3/mosdepth'
include { flagstat } from '../nf-modules/samtools/1.15/flagstat'
include { multiqc } from '../nf-modules/multiqc/1.14/multiqc'

// https://github.com/nextflow-io/nextflow/discussions/1975
def inner_join(ch_a, ch_b) {
    return ch_b.cross(ch_a).map { [it[0][0], *it[1][1..-1], *it[0][1..-1]] }
}

// Sub-workflow
workflow ALIGNMENT {
    main:

    // Outprefix is used as a parent directory here, so adjust the outpath
    def outdir = [params.outdir, params.out_prefix].join('/')

    // Intro text for multiqc reqport
    def intro = 'intro_text: QC of samples aligned to their respective genomes.'

    // Optional arguments
    def mq =  params.containsKey('mapq') ?
        params.mapq :
        0

    // Get sequence reads [ reads.basename, [ R1, R2] ]
    Channel
        .fromFilePairs(
            [params.seqdir.path, params.seqdir.pattern].join('/'),
            size: params.seqdir.nfiles,
        )
        .ifEmpty { exit 1, "Can't find read files." }
        .set { ch_reads }

    // Parse CSV [ read.basename, reference ]
    Channel
        .fromPath(params.sheet)
        .splitCsv()
        .map { tuple(it[0], it[1])}
        .set { ch_csv }

    // Reads object [ bn_ref, ref, bn_reads, [R1, R2] ]
    ch_csv
        .join(ch_reads)
        .map { 
            def f = file(it[1])
            tuple(f.simpleName,it[0],it[2])
        }
        .set { ch_reads_ref }

        // Get FAI files if they exist
        def Map fai_idx_status = Tools.hasIndex(params.sheet, 'fai')
        if (fai_idx_status[(true)]) {
            Channel
                .fromList(fai_idx_status[(true)])
                .map { ref ->
                    def f = file(ref)
                    tuple(f.simpleName,file(ref + '.fai'))
                }
                .set { ch_fai_exists }
        } else { Channel.empty().set { ch_fai_exists }}

        // Create FAI files if they don't
        if (fai_idx_status[(false)]) {
            // Get reference directory path
            Channel.fromList(fai_idx_status[(false)])
                .map { 
                    def f = file(it)
                    tuple(f.simpleName,f)
                }
                .set { ch_faidx_ref }

            faidx(ch_faidx_ref)
            faidx.out.set { ch_fai_created }
        } else { Channel.empty().set { ch_fai_created }}

        // Join all the FAI channels into a single entry
        //  [ id, file(fai) ]
        ch_fai_created.concat(ch_fai_exists).set { ch_fai }
    
    // Align reads using aligner of choice
    switch(params.aligner) {
        case "bwa2":

            // Get BWA2 index files if they exists
            def Map bwa_idx_status = Tools.hasIndex(params.sheet, 'bwa')
            if (bwa_idx_status[(true)]) {
                Channel
                    .fromList(bwa_idx_status[(true)])
                    .map { ref ->
                        // File f = new File(ref)
                        def f = file(ref)
                        tuple(
                            f.simpleName,
                            f,
                            [ 
                                file(ref + '.0123'),
                                file(ref + '.ann'),
                                file(ref + '.amb'),
                                file(ref + '.bwt.2bit.64'),
                                file(ref + '.pac')
                            ]
                        )
                    }
                    .set { ch_bwa_idx_exists }
            } else { Channel.empty().set { ch_bwa_idx_exists } }

            // Create BWA2 indicies if don't exist
            if (bwa_idx_status[(false)]) {
                Channel
                    .fromList(bwa_idx_status[(false)])
                    .map { ref ->
                        def f = file(ref)
                        tuple(f.simpleName,f)
                    }
                    .set { ch_bwa_ref }

                bwa_mem2_index(ch_bwa_ref)
                bwa_mem2_index.out.set { ch_bwa_idx }
            } else { Channel.empty().set { ch_bwa_idx } }

            // Create input channel to BWA-mem2
            // NOTE: Still have to wait for these files to be made, but 
            //       frees up a fair bit of disk space
            ch_bwa_idx.concat(ch_bwa_idx_exists).set { ch_index }
            inner_join(ch_reads_ref, ch_index).set { ch_tmp }
            inner_join(ch_tmp, ch_fai)
                .map {
                    // [ id, [R1, R2], asm, fai, [bwa2-indicies] ]
                    tuple(it[1], it[2], it[3], it[5], it[4])
                }
                .set { ch_input }

            // Run BWA2 alignment
            bwa_mem2_mem(
                ch_input,
                params.platform,
                mq,
                outdir + '/qc/flagstat/pre-filter'
            )
            bwa_mem2_mem.out.bam.set { ch_bam }
            bwa_mem2_mem.out.multiqc.set { ch_flagstat_pre }
            break;
        case "minimap2":
            // Minimap2 makes the ref index pretty quickly at runtime and doesn't
            // store it on disk. Saves having to make it ahead of time.
            ch_reads.join(ch_csv).set { ch_input }

            // Run alignment
            minimap2_sr(
                ch_input,
                params.platform,
                mq,
                outdir + '/qc/flagstat/pre-filter'
            )
            
            // Output channel
            minimap2_sr.out.bam.set { ch_bam }
            minimap2_sr.out.multiqc.set { ch_flagstat_pre }
            break;
    }

    // mark duplicates
    markdup(
        ch_bam,
        outdir + '/alignments',
        outdir + '/qc/sambamba'
    )

    // Depth statistics
    mosdepth(
        markdup.out.bam,
        outdir + '/qc/mosdepth'
    )

    // Alignment statistics
    flagstat(
        markdup.out.bam.map { val -> val[1]},
        outdir + '/qc/flagstat/post-filter'
    )

    // Aggregate all summary channels into a single data channel
    ch_flagstat_pre.collect().combine(markdup.out.multiqc.collect().combine(flagstat.out.multiqc.collect().combine(mosdepth.out.multiqc.collect())))
        .collect()
        .set { ch_multiqc }

    multiqc(
        ch_multiqc,
        projectDir + '/conf/multiqc-configs/mqc-alignment.yml',
        intro,
        outdir + '/qc/multiqc'
    )
}