/*
Genome assembly pipeline
    * Convert BAM file to fastq file
    * KMC K-mer analysis + GenomeScope2.0 for genome size estimation
	* Genome assembly (hifiasm/flye)
	* QUAST
	* BUSCO
	* KAT-compare
	* countGapular
*/

// Import pipeline modules - change to './nf-mod...' after it's all tested
include { kmc } from '../nf-modules/kmc/3.2.1/kmc'
include { genomescope } from '../nf-modules/genomescope/2.0/genomescope'
include { hifiasm } from '../nf-modules/hifiasm/0.16.1//hifiasm'
include { hifiasm_hic } from '../nf-modules/hifiasm/0.16.1/hifiasm-hic'
include { bwa_mem2_index } from "../nf-modules/bwa-mem2/2.2.1/bwa-mem2-index"
include { bwa_mem2_mem } from "../nf-modules/bwa-mem2/2.2.1/bwa-mem2-mem"
include { pin_hic } from '../nf-modules/pin_hic/3.0.0/pin_hic'
// include { flye } from '../nf-modules/assembly/flye'
// include { countgapular } from '../nf-modules/1.0/countgapular'

// include { katcompare } from '../nf-modules/assembly/katcompare'
// include { quast } from '../nf-modules/assembly/quast'
// include { busco } from '../nf-modules/assembly/busco'

// Sub-workflow
workflow ASSEMBLY {
    main:

    // HiFi Data
    Channel
        .fromFilePairs(
            [ params.input.path, params.input.pattern].join('/'), 
            size: params.input.nfiles
        )
        .ifEmpty { exit 1, "HiFi Fastq files are empty" }
        .set { ch_hifi }

    // Hi-C data + run assembly pipeline
    if (params.containsKey("hic")) {
        Channel
            .fromFilePairs(
                [params.hic.path, params.hic.pattern].join('/'),
                size: params.hic.nfiles
            )
            .set { ch_hic }

        // Run the assembly process
        hifiasm_hic(
            ch_hifi,
            ch_hic,
            params.outdir
        )

        // Create channel for each haplotype
        hifiasm_hic.out.hap_fa.flatten().map { val ->
            return tuple(val.baseName, val)
        }
        .set { ch_haplotypes }

        // Index reference files
        bwa_mem2_index(ch_haplotypes, params.outdir)
        
        // Join haplotype channels with index information
        ch_haplotypes.join(bwa_mem2_index.out.fai.join(bwa_mem2_index.out.bwa_idx)) .set { ch_tmp }

        // Combine idx with hic-read channel
        ch_tmp.combine(ch_hic).set { ch_hap_idx_hic }
        
        // Align data to reference genomes
        bwa_mem2_mem(
            ch_hap_idx_hic,
            params.outdir
        )

        // Hi-c scaffolding
        pin_hic(bwa_mem2_mem.out.bam, params.outdir)

    } else {
        // Run the assembly process
        hifiasm(
            ch_hifi,
            params.outdir
        )
    }

    // Genome size estimation
    // TODO: Check this
    kmc(ch_hifi, params.outdir)
    genomescope(kmc.out.histo, params.outdir)

}