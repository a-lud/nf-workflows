/*
Variant calling pipeline

This is an attempt to make a 'pretty' generic variant calling pipeline
using BCFtools. The user should have most of the options required to
call variants on single samples or jointly, along with filtering the VCF
as required.
*/

// Import pipeline functions
include { coverage } from '../nf-modules/samtools/1.15/coverage'
include { joint_mpileup_call } from '../nf-modules/bcftools/1.15.1/joint_mpileup_call'
include { standard_mpileup_call } from '../nf-modules/bcftools/1.15.1/standard_mpileup_call'
include { joint_filter } from '../nf-modules/bcftools/1.15.1/joint_filter'
include { standard_filter } from '../nf-modules/bcftools/1.15.1/standard_filter'
include { concat } from '../nf-modules/bcftools/1.15.1/concat'

// Sub-workflow
workflow VARIANT {
    main:

    /*
    Optional arguments
        * Setting some useful default values for optional variables
    */

    // Optional argument values
    def nodp = params.containsKey('nodp') ? 
        params.nodp :
        false

    def mapq =  params.containsKey('mapq') ?
        params.mapq :
        10
    
    def baseq = params.containsKey('baseq') ?
        params.baseq :
        10

    def genoq = params.containsKey('genoq') ?
        params.genoq :
        10

    def ploidy = params.containsKey('ploidy') ?
        params.ploidy :
        2

    def mpileupOpt = params.containsKey('mpileup_opt') ?
        params.mpileup_opt :
        ''
    
    def callOpt = params.containsKey('call_opt') ?
        params.call_opt :
        ''
    
    def filterOpt = params.containsKey('filter_opt') ?
        params.filter_opt :
        ''

    def viewOpt = params.containsKey('view_opt') ?
        params.view_opt :
        ''

    def normOpt = params.containsKey('norm_opt') ?
        params.norm_opt :
        ''
    
    def sortOpt = params.containsKey('sort_opt') ?
        params.sort_opt :
        ''

    /*
    Output directories
    */

    def outdir = [params.outdir, params.out_prefix].join('/')
    def outtmp = [outdir, 'junk'].join('/')
    def outbml = [outdir, 'bamlists'].join('/')
    def outcov = [outdir, 'coverage-stats'].join('/')
    def outvcf = [outdir, 'vcf'].join('/')
    def outflt = [outdir, 'vcf-filtered'].join('/')
    def outcat = [outdir, 'vcf-final'].join('/')

    /*
    Create a junk 'no_regions' file to prevent broken symlinks breaking output channels
    Need to think of a way to remove this once the pipeline is done...
    */

    file(outtmp).mkdirs()
    def removeThis = new File("${outtmp}")
    def junkFile = new File("${outtmp}/no_regions")
    junkFile.createNewFile()

    /*
    Sort out the input data channels
    Parse sample-sheet CSV [ bam.baseName, [BAM], [reference], [regions] ]
    */
    Channel
        .fromPath(
            params.sheet
        )
        .splitCsv()
        .map { row -> 
            bam = file(row[0])
            bai = file("${row[0]}.bai") // I'm sure this will come back to bite me by not finding the files...
            ref = file(row[1])
            
            tuple( tuple(file(bam), file(bai)), file(ref) )
        }
        .set { ch_tmp }
    
    // Check regions - this dicates the parallel nature below
    //      - Output channel: [ID, ref, <regions>]
    if ( params.containsKey('regions') ) {

        // File channel of regions files
        Channel
            .fromFilePairs(
                [params.regions.path, params.regions.pattern].join('/'),
                size: params.regions.nfiles,
            )
            .set { ch_r }

        // Reference channel
        ch_tmp
            .map { tuple(it[1].baseName,it[1])}
            .unique()
            .set { ch_ref }

        // Joint method should be parallelised based on chromosomes in region file
        if (params.method == 'joint') {

            // Get chromsomes from regions file
            ch_r
                .map { tuple(it[0], it[1][0])}
                .splitCsv(elem: 1, sep: "\t")
                .map {
                    tuple(it[0], it[1][0])
                }
                .unique()
                .groupTuple()
                .set { ch_chr }

            // Join and transpose so we end up with [ refid, reference, chr ] for every combination
            ch_ref
                .join(ch_chr, by: 0, remainder: true)
                .transpose()
                .set { ch_regions }

            /*
            Create a dummy channel so there is something to join on.
            Prevents the dropping of channel elements when there is no key to join on.
            */
            ch_ref
                .map { tuple(it[0], null)}
                .concat(ch_r)
                .groupTuple()
                .map { id, bed ->
                    if (bed.size() > 1) {
                        b = bed - null
                        return tuple(id, b[0][0])
                    }
                    return tuple(id, file("${outtmp}/no_regions"))
                }
                .set { ch_dummy }

            ch_regions
                .combine(ch_dummy, by: 0)
                .set { ch_regions }

        } else {
            // Join with data channels - set as no_regions if no regions
            ch_ref
                .join(ch_r, by: 0, remainder: true)
                .map {
                    it[2] = it[2] ?: file("${outtmp}/no_regions")
                    return it
                }
                .set { ch_regions }
        }
    } else {
        // Reference channel
        ch_tmp
            .map { tuple(it[1].baseName,it[1])}
            .unique()
            .set { ch_ref }

        // Set no regions (null or file depending on method)
        if (params.method == 'joint') {
            ch_ref
                .map {
                    tuple(it[0], it[1], null, file("${outtmp}/no_regions"))
                }
                .set { ch_regions }

        } else {
            ch_ref
                .map {
                    tuple(it[0], it[1], file("${outtmp}/no_regions"))
                }
                .set { ch_regions }
        }
    }

    /*
    The design here is having a channel of the form:
        
        - Joint genotyping
            * [ ID, [BAM/S], BAMLIST, REFERENCE, CHR, REG.BED ]
        - Standard genotyping
            * [ ID, BAM, REF, REG.BED]
    
    Where ID will either be 
        - BAM basename if in standard mode
        - REF basename if in join mode

    In either case, we just need a unique identifier to have for the output.
    Each of these IDs make sense for their own variant calling approach.
    */

    // Manage input channels depending on variant calling type
    if (params.method == 'joint') {

        // Group BAM files by reference genomes - [ [BAM1, BAM2, BAM..., BAMN], reference ]
        ch_tmp
            .groupTuple(by: 1)
            .map { bams, ref ->
                tuple(ref.baseName, bams)
            }
            .set { ch_id_bam }

        // Create output directory
        file(outbml).mkdirs()
        
        // Make a bam-list file for each reference
        ch_id_bam
            .map { refid, bams ->
                bamListFile = file("${outbml}/${refid}.bamlist") //  [outbml, refid + '.bamlist'].join('/') 
                bams.each {
                    bamListFile << it[0].getName() + '\n'
                }
            }
            .map { 
                // Create the termination file 
                file( "$outbml/zzz.bamlist") << '\n'
            }

        Channel
            .watchPath("${outbml}/*.bamlist", 'create,modify')
            .until { file -> file.baseName == 'zzz' }
            .map {
                tuple(it.baseName, it)
            }
            .set { ch_bamlist }
        
        // Join bamlist to reference - bams channel (on reference id key)
        ch_id_bam
            .join(ch_bamlist)
            .combine(ch_regions, by: 0)
            .set { ch_data }

    // Standard mode - calling variants for individual samples
    } else {
        ch_tmp
            .map { bam, ref ->
                tuple(ref.baseName, bam)
            }
            .combine(ch_regions, by: 0)
            .map {
                tuple(it[1][0].baseName, it[1][0], it[1][1], it[2], it[3])
            }
            .set { ch_data }
    }

    /*
    Call variants and filter
    */

    // Make BAM tuple (bam and idx) into separate tuples
    if (params.method == 'standard') {

        /*
        Coverage information per sample
        */
        ch_data
            .map { 
                tuple(it[0], it[1], it[3])
            }
            .set { ch_cov }

        coverage(ch_cov, outcov)

        // Genotype and call variants
        standard_mpileup_call(
            ch_data,
            params.vcftype,
            mapq,
            baseq,
            ploidy,
            mpileupOpt,
            callOpt,
            outvcf
        )

        // Join each samples respective coverage information
        if (nodp) {
            standard_mpileup_call.out.vcf
                .map { tuple(it[0], it[1], it[2], it[3], null, null) }
                .set { ch_flt }
        } else {
            standard_mpileup_call.out.vcf
                .join(coverage.out.covdepth)
                .set { ch_flt }
        }

        // Filter variants
        standard_filter(
            ch_flt,
            filterOpt,
            viewOpt,
            normOpt,
            sortOpt,
            outflt
        )
    } else {
        
        /*
        Coverage information per-samples
        */
        ch_data
            .map {
                // TODO: This causes java.util.ConcurrentModificationException.
                def lst = []
                it[1].each { bam_idx ->
                    lst << bam_idx[0]
                }
                tuple(lst, it[3])
            }
            .transpose()
            .unique()
            .map {tuple(it[0].baseName, it[0], it[1])}
            .set { ch_cov }

        coverage(ch_cov, outcov)

        // One final channel manipulation to get BAM files into a single
        // 'list' and their BAI files into another. Required for linking
        // into the working directory (avoids the 'input1.bam' from tmp issue).
        ch_data
            .map {
                def bams = []
                def idx = []
                it[1].each { tup ->
                    bams << file(tup[0])
                    idx << file(tup[1])
                }
                tuple(
                    it[0],
                    bams,
                    idx,
                    it[2],
                    it[3],
                    it[4],
                    it[5]
                )
            }
            .set { ch_data }

        // Genotype and call variants - each channel input is a separate chromosome
        joint_mpileup_call(
            ch_data,
            params.vcftype,
            mapq,
            baseq,
            ploidy,
            mpileupOpt,
            callOpt,
            outvcf
        )

        // Filter the chromosome VCF
        joint_filter(
            joint_mpileup_call.out.vcf,
            filterOpt,
            viewOpt,
            normOpt,
            sortOpt,
            outflt
        )

        // Concatenate the chromosome VCF files together into a single output
        joint_filter.out.vcf
            .groupTuple()
            .map {
                vcfs = []
                idxs = []
                it[1].each { i ->
                    vcfs << i[0]
                    idxs << i[1]
                }
                tuple(it[0], vcfs, idxs)
            }
            .set { ch_vcfs }
        
        concat(
            ch_vcfs,
            outcat
        )
    }

    // Remove temp directory on complete
    workflow.onComplete {
        removeThis.deleteDir()
    }
}