/*
Genome assembly pipeline
    * Remove adapters from HiFi Fastq
    * Convert FASTQ to FASTA (Used by TgsGapcloser in next pipeline)
    * KMC K-mer analysis + GenomeScope2.0 for genome size estimation
	* Genome assembly (hifiasm)
	* Hi-C read processing
    * Scaffolding (salsa2/pin_hic)
    * Genearte Juicebox visualisation files
	* BUSCO (contigs/scaffolds)

NOTE: Use the 'assembly_assessment' pipeline to finalise the genome with gap-closing,
      followed by a range of QC assessments
*/

// Import pipeline modules
include { hifiadapterfilt } from '../nf-modules/hifiadapterfilt/2.0.0/hifiadapterfilt'
include { seqkit_fq2fa } from '../nf-modules/seqkit/2.1.0/seqkit_fq2fa'
include { kmc } from '../nf-modules/kmc/3.2.1/kmc'
include { genomescope } from '../nf-modules/genomescope/2.0/genomescope'
include { hifiasm } from '../nf-modules/hifiasm/0.16.1//hifiasm'
include { hifiasm_hic } from '../nf-modules/hifiasm/0.16.1/hifiasm-hic'
include { bwa_mem2_index } from "../nf-modules/bwa-mem2/2.2.1/bwa-mem2-index"
include { arima_map_filter_combine } from "../nf-modules/arima/1.0.0/arima_map_filter_combine"
include { arima_dedup_sort } from '../nf-modules/arima/1.0.0/arima_dedup_sort'
include { pin_hic } from '../nf-modules/pin_hic/3.0.0/pin_hic'
include { salsa2 } from '../nf-modules/salsa2/2.3/salsa2'
include { matlock_bam2 } from '../nf-modules/matlock/20181227/matlock_bam2'
include { assembly_visualiser as assembly_visualiser_pin_hic } from '../nf-modules/3d_dna/201008/assembly_visualiser'
include { assembly_visualiser as assembly_visualiser_salsa2 } from '../nf-modules/3d_dna/201008/assembly_visualiser'
include { busco as busco_contig } from '../nf-modules/busco/5.2.2/busco'
include { busco as busco_salsa2 } from '../nf-modules/busco/5.2.2/busco'
include { busco as busco_pin_hic } from '../nf-modules/busco/5.2.2/busco'

// Sub-workflow
workflow ASSEMBLY {
    main:

    // Output directory paths
    out_filteredreads = params.outdir + '/adapter-removed-reads'
    out_hifiasm = params.outdir + '/assembly-contigs/' + params.out_prefix
    out_busco = params.outdir + '/qc/busco'
    out_scaffold = params.outdir + '/assembly-scaffold'
    out_genomescope = params.outdir + '/genome-size'

    // HiFi Data
    Channel
        .fromFilePairs(
            [ params.hifi.path, params.hifi.pattern].join('/'), 
            size: params.hifi.nfiles
        )
        .ifEmpty { exit 1, "Long-read fastq file channel is empty. Can't find the files..." }
        .set { ch_hifi }
    
    // HifiAdapterFilt: Remove adapters from HiFi
    hifiadapterfilt(
        ch_hifi,
        out_filteredreads
    )

    // // FQ2FA: hifi reads
    seqkit_fq2fa(
        hifiadapterfilt.out.clean, 
        out_filteredreads
    )

    // Hi-C data + run assembly pipeline
    if (params.containsKey("hic")) {
        Channel
            .fromFilePairs(
                [params.hic.path, params.hic.pattern].join('/'),
                size: params.hic.nfiles
            )
            .set { ch_hic }
        
        // Hifiasm: assemble reads into contigs
        hifiasm_hic(
            hifiadapterfilt.out.clean,
            ch_hic,
            params.out_prefix,
            out_hifiasm
        )
        
        // String to match on - not pretty but does the trick
        switch(params.assembly) {
            case 'all':
                pattern = ['p_ctg', 'hap1', 'hap2' ]
                break;
            case 'primary':
                pattern = [ 'p_ctg' ]
                break;
            case 'haplotype1':
                pattern = [ 'hap1' ]
                break;
            case 'haplotype2':
                pattern = [ 'hap2' ]
                break;
            case 'haplotypes':
                pattern = [ 'hap1' ,'hap2' ]
                break;
        }

        // Filter Hifiasm output channels for only the assemblies we're after
        hifiasm_hic.out.contigs.flatten().filter{
            pattern.any { val -> it.baseName.contains(val)}
        }
        .map { ctg ->
            return tuple(ctg.baseName, ctg)
        }
        .set { ch_contigs }

        // BUSCO: contig assemblies
        busco_contig(
            ch_contigs,
            params.busco_db,
            'contig',
            out_busco
        )

        // BWA2: Index reference files
        bwa_mem2_index(
            ch_contigs
        )

        // Combine idx with hic-read channel
        bwa_mem2_index.out.combine(ch_hic).set { ch_hap_idx_hic }
        
        // Arima Hi-C processing: Align + filter hic reads
        arima_map_filter_combine(ch_hap_idx_hic)
        arima_dedup_sort(arima_map_filter_combine.out.bam)

        // Maaaaaatlooooooock: Convert Hi-C to contig BAM file to 'merged_nodups.txt' file used by juicebox
        matlock_bam2(arima_dedup_sort.out.hic_to_ctg_bam)

        // Combine processed Hi-C reads with genome again - join on genome ID (first field)
        ch_contigs.join(arima_dedup_sort.out.bam).set { ch_ref_hic }

        // Hi-c scaffolding
        switch(params.scaffolder) {
            case 'all':
                // Scaffold
                pin_hic(ch_ref_hic, out_scaffold)
                salsa2(ch_ref_hic, out_scaffold)

                // BUSCO on scaffolds ---
                busco_salsa2(
                    salsa2.out.scaffolds,
                    params.busco_db,
                    'scaffold-salsa2',
                    out_busco
                )
                busco_pin_hic(
                    pin_hic.out.scaffolds,
                    params.busco_db,
                    'scaffold-pin_hic',
                    out_busco
                )

                // Join scaffold agp with matlock output
                pin_hic.out.juicebox.join(matlock_bam2.out).set { ch_pin_hic_juicebox }
                salsa2.out.juicebox.join(matlock_bam2.out).set { ch_salsa2_juicebox }

                // Create Juicebox files
                assembly_visualiser_pin_hic(
                    ch_pin_hic_juicebox, 
                    'pin_hic',
                    out_scaffold
                )
                assembly_visualiser_salsa2(
                    ch_salsa2_juicebox,
                    'salsa2', 
                    out_scaffold
                )
                break;
            case 'salsa2':
                salsa2(ch_ref_hic, out_scaffold)
                
                // BUSCO on scaffolds ---
                busco_salsa2(
                    salsa2.out.scaffolds,
                    params.busco_db,
                    'scaffold-salsa2',
                    out_busco
                )
                
                // Generate Juicebox input files ---
                salsa2.out.juicebox.join(matlock_bam2.out).set { ch_salsa2_juicebox }
                assembly_visualiser_salsa2(
                    ch_salsa2_juicebox,
                    'salsa2',
                    out_scaffold
                )
                break;
            case 'pin_hic':
                pin_hic(ch_ref_hic, out_scaffold)
                
                // BUSCO on scaffolds ---
                busco_pin_hic(
                    pin_hic.out.scaffolds,
                    params.busco_db,
                    'scaffold-pin_hic',
                    out_busco
                )
                
                // Generate Juicebox input files ---
                pin_hic.out.juicebox.join(matlock_bam2.out).set { ch_pin_hic_juicebox }
                assembly_visualiser_pin_hic(
                    ch_pin_hic_juicebox,
                    'pin_hic',
                    out_scaffold
                )
                break;
        }

        /*
        PIPELINE BREAK: The pipeline should break here for users to edit their genome assemblies.
        - Might look at adding functionality to just continue with the pipeline? but feel
          like users should always check and edit their assemblies in Juicebox to ensure no
          obvious misjoins.
        */

    } else {
        // Hifiasm: assemble contigs
        hifiasm(
            ch_hifi,
            params.out_prefix,
            out_hifiasm
        )

        // BUSCO: contig primary assembly
        busco_contig(
            hifiasm.out.contigs,
            params.busco_db,
            'contig',
            out_busco
        )
        
        // BUSCO plot: Generate a summary plot of the BUSCO results
        // busco_plot(busco_contig.out.summary, params.outdir)

        // Merqury: K-mer assessment - Compare genome to reads
        // hifiasm.out.contigs.map { id, val -> return file(val)}
        // .set { ch_contig }

        // ch_contig.combine(ch_hifi).set { ch_contig_hifi }

        // merqury(ch_contig_hifi, params.outdir)

        // // QUAST: contig primary assembly
        // quast(ch_contig, params.outdir)

        // // Mosdepth - HIFI alignment and coverage
        // hifiasm.out.fa.combine(ch_hifi).set { ch_contig_hifi_fq }
        // minimap2_pb_hifi(ch_contig_hifi_fq, params.outdir)
        // mosdepth(minimap2_pb_hifi.out, params.outdir)
    }

    // Genome size estimation
    kmc(ch_hifi, params.out_prefix)
    genomescope(kmc.out.histo, out_genomescope)
}
