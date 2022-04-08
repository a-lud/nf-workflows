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
include { seqkit_fq2fa } from '../nf-modules/seqkit/2.1.0/seqkit_fq2fa'
include { kmc } from '../nf-modules/kmc/3.2.1/kmc'
include { genomescope } from '../nf-modules/genomescope/2.0/genomescope'
// include { hifiasm } from '../nf-modules/hifiasm/0.16.1//hifiasm'
include { hifiasm_hic } from '../nf-modules/hifiasm/0.16.1/hifiasm-hic'
include { bwa_mem2_index } from "../nf-modules/bwa-mem2/2.2.1/bwa-mem2-index"
include { arima_map_filter_combine } from "../nf-modules/arima/1.0.0/arima_map_filter_combine"
include { arima_dedup_sort } from '../nf-modules/arima/1.0.0/arima_dedup_sort'
include { bwa_mem2_mem } from "../nf-modules/bwa-mem2/2.2.1/bwa-mem2-mem"
include { pin_hic } from '../nf-modules/pin_hic/3.0.0/pin_hic'
include { salsa2 } from '../nf-modules/salsa2/2.3/salsa2'
include { tgsgapcloser } from '../nf-modules/tgs-gapcloser/1.1.1/tgsgapcloser'
include { busco as busco_contig } from '../nf-modules/busco/5.2.2/busco'
include { busco as busco_scaffold } from '../nf-modules/busco/5.2.2/busco'
include { busco_plot } from '../nf-modules/busco/5.2.2/busco_plot'
include { quast } from '../nf-modules/quast/5.0.2/quast'
include { merqury } from '../nf-modules/merqury/1.3/merqury'
include { mosdepth } from '../nf-modules/mosdepth/0.3.3/mosdepth'
include { minimap2_pb_hifi } from '../nf-modules/minimap2/2.24/minimap2_pb_hifi'

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

    // FQ2FA - hifi reads
    seqkit_fq2fa(ch_hifi, params.outdir)

    // Hi-C data + run assembly pipeline
    if (params.containsKey("hic")) {
        Channel
            .fromFilePairs(
                [params.hic.path, params.hic.pattern].join('/'),
                size: params.hic.nfiles
            )
            .set { ch_hic }      
        
        // Hifiasm - assemble reads into contigs
        hifiasm_hic(
            ch_hifi,
            ch_hic,
            params.outdir
        )

        // Create channel for each haplotype
        hifiasm_hic.out.hap_fa.flatten().map { val ->
            return tuple(val.baseName, val)
        }
        .set { ch_contigs }

        // BUSCO - contig haplotype assemblies
        busco_contig(ch_contigs, params.busco_db, 'contig', params.outdir)

        // Index reference files
        bwa_mem2_index(ch_contigs, params.outdir)
        
        // Join haplotype channels with index information
        ch_contigs.join(bwa_mem2_index.out.fai.join(bwa_mem2_index.out.bwa_idx)).set { ch_tmp }

        // Combine idx with hic-read channel
        ch_tmp.combine(ch_hic).set { ch_hap_idx_hic }
        
        // Arima Hi-C processing
        arima_map_filter_combine(ch_hap_idx_hic, params.outdir)
        arima_dedup_sort(arima_map_filter_combine.out.bam, params.outdir)

        // Combine processed Hi-C reads with genome again - join on genome ID
        ch_contigs.join(arima_dedup_sort.out.bam).set { ch_ref_hic }

        // Hi-c scaffolding
        switch(params.scaffolder) {
            case 'all':
                // Scaffold
                pin_hic(ch_ref_hic, params.outdir)
                salsa2(ch_ref_hic, params.outdir)
                break;
            case 'salsa2':
                salsa2(ch_ref_hic, params.outdir)

                break;
            case 'pin_hic':
                pin_hic(ch_ref_hic, params.outdir)
                break
        }

        // Gap closing: TGS-GapCloser
        pin_hic.out.scaffolds.combine(seqkit_fq2fa.out).set { ch_scaff_hifi_fa }
        tgsgapcloser(ch_scaff_hifi_fa, params.outdir)

        // Mosdepth - HIFI alignment and coverage
        tgsgapcloser.out.scaff_detailed.combine(ch_hifi).set { ch_filled_hifi_fq }
        minimap2_pb_hifi(ch_filled_hifi_fq, params.outdir)
        mosdepth(minimap2_pb_hifi.out, params.outdir)

        // K-mer assessment - Merqury
        tgsgapcloser.out.scaff.collect().toList().combine(ch_hifi).set { ch_filled_hifi }
        merqury(ch_filled_hifi, params.outdir)

        // BUSCO - Scaffold haplotype assemblies
        busco_scaffold(tgsgapcloser.out.scaff_detailed, params.busco_db, 'gapfilled', params.outdir)

        // QUAST - scaffold haplotype assemblies
        quast(tgsgapcloser.out.scaff.collect(), params.outdir)

        // BUSCO plot - Generate a summary plot of the BUSCO results (contig and scaffold)
        busco_contig.out.summary.concat(busco_scaffold.out.summary).collect().set { ch_busco_short }
        busco_plot(ch_busco_short, params.outdir)
        
    } else {
        // TODO: TEST THIS PORTION OF THE PIPELINE
        // Run the assembly process
        hifiasm(
            ch_hifi,
            params.outdir
        )

        // BUSCO - contig primary assembly
        busco_contig(hifiasm.out.fa, params.busco_db, 'contig', params.outdir)
        
        // BUSCO plot - Generate a summary plot of the BUSCO results
        busco_plot(busco_contig.out.summary, params.outdir)

        // K-mer assessment - Compare genome to reads
        hifiasm.out.fa.map { id, val ->
            return file(val)
        }
        .set { ch_contig }

        ch_contig.combine(ch_hifi).set { ch_contig_hifi }

        merqury(ch_contig_hifi, params.outdir)

        // QUAST - contig primary assembly
        quast(ch_contig, params.outdir)

        // Mosdepth - HIFI alignment and coverage
        hifiasm.out.fa.combine(ch_hifi).set { ch_contig_hifi_fq }
        minimap2_pb_hifi(ch_contig_hifi_fq, params.outdir)
        mosdepth(minimap2_pb_hifi.out, params.outdir)
    }

    // Genome size estimation
    kmc(ch_hifi, params.outdir)
    genomescope(kmc.out.histo, params.outdir)
}