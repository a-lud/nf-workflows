/*
Ortholog detection pipeline

This pipeline is designed to find orthologs between species with genome
assemblies and annotation files (in GFF3). The pipeline aims to use best practice
methods to obtain a high-quality ortholog geneset out the other end. The pipeline is
as follows:
    1. GFF statistics (All isoforms + longest isoform performed by default) - AGAT agat_sp_statistics.pl
    2. Filter for longest isoform per gene - AGAT agat_sp_keep_longest_isoform.pl
    3. Extract CDS for longest isoform - AGAT agat_sp_extract_sequences.pl
        * Protein
        * Nucleotide
    4. Ortholog detection - OrthoFinder
    5. MSA to codon alignment
    6. MSA statistics pre-cleaning - msaSummary
    7. Clean MSA - Gblocks/CIalign
    8. MSA statistics post-cleaning - msaSummary

Potential addition?
    * Leave-one-out analysis to see which sample affects the results the most
*/

// Modules required by sub-workflow
include { agat_statistics } from '../nf-modules/agat/0.9.2/agat_statistics'
include { agat_longest_iso } from '../nf-modules/agat/0.9.2/agat_longest_iso'
include { agat_extract_seq } from '../nf-modules/agat/0.9.2/agat_extract_seq'
include { orthofinder } from '../nf-modules/orthofinder/2.5.2/orthofinder'
include { prot_to_codon_msa } from '../nf-modules/general/prot_to_codon_msa'
include { msaSummary as msaSummary_pre } from '../nf-modules/general/msaSummary'
include { msaSummary as msaSummary_post } from '../nf-modules/general/msaSummary'
include { clipkit } from '../nf-modules/clipkit/1.3.0/clipkit'

// Sub-workflow
workflow ORTHOFINDER {
    main:

    // Outprefix is used as a parent directory here, so adjust the outpath
    def outdir = [params.outdir, params.out_prefix].join('/')

    // Get GFF files
    Channel
        .fromFilePairs(
            [params.gffs.path, params.gffs.pattern].join('/'),
            size: params.gffs.nfiles,
        )
        .ifEmpty { exit 1, "Can't find GFF3 files." }
        .set { ch_gffs }

    // Get Genome files
    Channel
        .fromFilePairs(
            [params.genomes.path, params.genomes.pattern].join('/'),
            size: params.genomes.nfiles,
        )
        .ifEmpty { exit 1, "Can't find genome files." }
        .set { ch_asm }

    // Optional arguments for OrthoFinder
    def tree = (params.containsKey('tree') && params.tree) ? "-s " + params.tree : ''
    def trim_msa = (params.containsKey('trim_msa') && params.trim_msa) ? '' : '-z'
    def stop_early = (params.containsKey('stop_early') && params.stop_early ) ? '-oa' : ''

    // Statistics channel
    ch_gffs.join(ch_asm).set { ch_statistics }
    agat_statistics(ch_statistics, outdir)

    // Keep longest isoform
    agat_longest_iso(ch_gffs, outdir)
    agat_longest_iso.out.longest.join(ch_asm).set { ch_longest }

    // Extract longest seq
    agat_extract_seq(ch_longest, outdir)
    agat_extract_seq.out.protein.collect().set { ch_proteins }
    
    // Find orthologs with OrthoFinder
    orthofinder(
        ch_proteins, 
        params.search_prog,
        tree,
        stop_early, 
        trim_msa,
        outdir
    )

    // Generate codon MSA files from protein alignments
    agat_extract_seq.out.nucleotide.collect().set { ch_nucleotide }
    prot_to_codon_msa(
        ch_nucleotide,
        orthofinder.out.msa,
        outdir
    )

    // Summarise gaps pre-trimming
    msaSummary_pre(prot_to_codon_msa.out.msa_codon, 'pre', outdir)

    // Clean alignments
    clipkit(prot_to_codon_msa.out.msa_codon, outdir)

    // Summarise gaps post-trimming
    msaSummary_post(clipkit.out.clean, 'post', outdir)
}