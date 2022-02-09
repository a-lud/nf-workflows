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
include { busco as busco_contig } from '../nf-modules/busco/5.2.2/busco'
include { busco as busco_scaffold } from '../nf-modules/busco/5.2.2/busco'
include { quast } from '../nf-modules/quast/5.0.2/quast'
include { kat_compare } from '../nf-modules/kat/2.4.2/kat_compare'
// include { countgapular } from '../nf-modules/1.0/countgapular'
// include { flye } from '../nf-modules/assembly/flye'

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

        // BUSCO - contig haplotype assemblies
        busco_contig(ch_haplotypes, params.busco_db, 'contig', params.outdir)

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

        // KAT - Compare genome to reads
        pin_hic.out.scaffolds.take(2).combine(ch_hifi).set { ch_scaff_hifi }
        kat_compare(ch_scaff_hifi, params.outdir)

        // BUSCO - Scaffold haplotype assemblies
        busco_scaffold(pin_hic.out.scaffolds.take(2), params.busco_db, 'scaffold', params.outdir)

        // QUAST - scaffold haplotype assemblies
        quast(pin_hic.out.scaffolds.take(2), params.outdir)

    } else {
        // Run the assembly process
        hifiasm(
            ch_hifi,
            params.outdir
        )

        // BUSCO - contig primary assembly
        busco_contig(hifiasm.out.fa, params.busco_db, 'contig', params.outdir)

        // KAT - Compare genome to reads
        hifiasm.out.fa.combine(ch_hifi).set { ch_contig_hifi }
        kat_compare(ch_contig_hifi, params.outdir)

        // QUAST - contig primary assembly
        quast(hifiasm.out.fa, params.outdir)
    }

    // Genome size estimation
    kmc(ch_hifi, params.outdir)
    genomescope(kmc.out.histo, params.outdir)
}