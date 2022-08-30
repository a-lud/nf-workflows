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
include { minimap2_sr } from '../nf-modules/minimap2/2.24/minimap2_sr'
include { bwa_mem2_index } from '../nf-modules/bwa-mem2/2.2.1/bwa-mem2-index'
include { bwa_mem2_mem } from '../nf-modules/bwa-mem2/2.2.1/bwa-mem2-mem'
include { markdup } from '../nf-modules/sambamba/0.8.2/markdup'
include { mosdepth } from '../nf-modules/mosdepth/0.3.3/mosdepth'
include { flagstat } from '../nf-modules/samtools/1.15/flagstat'
include { multiqc } from '../nf-modules/multiqc/1.12/multiqc'

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
    
    // Align reads using aligner of choice
    switch(params.aligner) {
        case "bwa2":
        // Parse CSV [ read.basename, reference ]
        Channel
            .fromPath(
                params.sheet
            )
            .splitCsv()
            .map { row ->
                tuple(row[0], row[1])
            }
            .set { ch_csv }

        // Unique reference files [ ref.basename, [ref] ]
        ch_csv
            .map {it[1]}
            .map {val ->
                File f = new File(val)
                tuple(
                    f.name.take(f.name.lastIndexOf('.')),
                    file(val)
                )
            }
            .unique()
            .set { ch_ref }

        // Index reference files [ ref.basename, [fai], [bwa2 idx files] ]
        bwa_mem2_index(
            ch_ref,
            outdir
        )

        // Join CSV data with index reference
        ch_csv
            .map {
                File f = new File(it[1])
                tuple(
                    it[0],
                    f.name.take(f.name.lastIndexOf('.'))
                )
            }
            .join(
                ch_reads
            )
            .map {tuple(it[1], it[0], it[2])} // [asm.bn, id, [reads]]
            .combine(
                bwa_mem2_index.out,
                by: 0
            )
            .map {
                tuple(
                    it[1], // id
                    it[2], // reads
                    it[3], // asm
                    it[4], // asm fai
                    it[5]  // bwa2 idx files
                )
            }
            .set { ch_input } // [ id, [reads], [asm], [fai], [bwa2 idx files] ]

            bwa_mem2_mem(
                ch_input,
                params.platform,
                mq
            )
            bwa_mem2_mem.out.set { ch_bam }
            break;
        case "minimap2":
            // Parse CSV [ read.basename, reference ]
            Channel
                .fromPath(
                    params.sheet
                )
                .splitCsv()
                .map { row ->
                    tuple(row[0], file(row[1]))
                }
                .set { ch_csv }
            ch_reads.join(ch_csv).set { ch_input }

            minimap2_sr(
                ch_input,
                params.platform,
                mq
            )
            minimap2_sr.out.set { ch_bam }
            break;
    }

    // mark duplicates
    markdup(
        ch_bam,
        outdir
    )

    // Depth statistics
    mosdepth(
        markdup.out.bam,
        outdir + '/qc/mosdepth'
    )

    // Alignment statistics
    flagstat(
        markdup.out.bam.map { val -> val[1]},
        outdir + '/qc/alignment-statistics'
    )

    // Aggregate all summary channels into a single data channel
    markdup.out.multiqc.collect().combine(mosdepth.out.multiqc.collect().combine(flagstat.out.multiqc.collect())).collect().set { ch_multiqc }

    multiqc(
        ch_multiqc,
        intro,
        outdir + '/qc/mosdepth'
    )
}