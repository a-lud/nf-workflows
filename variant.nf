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
    def caller = params.caller == 'new' ? '-m' : '-c'

    def force = params.containsKey('force') ?
        params.force :
        false

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
    def outreg = [outdir, 'regions'].join('/')
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

    file("${outtmp}/no_regions") << '\n'

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

    /*
    Create two reference channels
        * ch_ref_fa = [ ref.basename, ref.fasta ]
        * ch_ref = [ ref.basename ]
    */
    ch_tmp
        .map { tuple(it[1].baseName,it[1])}
        .unique()
        .tap  { ch_ref_fa }
        .map { it[0] }
        .tap { ch_ref }
    
    // Check regions - this dicates the parallel nature below
    //      - Output channel: [ID, ref, <regions>]
    if ( params.containsKey('regions') ) {

        /*
        Default regions channel - [ bed.basename, file ]
        */
        Channel
            .fromFilePairs(
                [params.regions.path, params.regions.pattern].join('/'),
                size: params.regions.nfiles,
            )
            .set { ch_regions }

        /*
        CheckRegionsSubset will check if
            1. the regions sub-directory exists. If it doesn't it returns false
            2. If there are files within the regions directory. Returns false if there are none.
        */

        def checkRegionSubset = Tools.checkDirAndFiles(outreg)

        // Joint method should be parallelised based on chromosomes in region file
        if (params.method == 'joint') {

            /*
            Subset BED file by chromosome if the files don't exist/force isn't used
                * The whole point of doing this in the script is to enable resume
                  functionality. Filtering for the chromosomes in the process
                  messess up the resume flag.
            */

            // Split on chromosome because force/the files don't exist already
            //  * Channel returned [ ref.basename, ref chr, bed.filepath ]
            if ( (!checkRegionSubset) || (force) ) {

                // Create output directory
                file(outreg).mkdirs()

                // Split region BED files up by chromosome
                ch_regions
                    .join(ch_ref)
                    .map { tuple(it[0], it[1][0]) }
                    .splitCsv(elem: 1, sep: "\t")
                    .groupTuple()
                    .map { refid, regions ->
                    // Write regions to respecitve BED files
                        regions.each {
                            regFile = file("${outreg}/${refid}:${it[0]}.bed")
                            regFile << it.join('\t') + "\n"
                        }
                    }
                    .map { file( "${outreg}/DONE.bed") << '\n' }

                // Read in the now split BED files
                Channel
                    .watchPath("${outreg}/*.bed", 'create,modify')
                    .until { file -> file.baseName == 'DONE' }
                    .map {
                        def sp = it.baseName.split(":")
                        tuple(sp[0], sp[1], file(it) ) 
                    }
                    .set { ch_regions_by_chrom }

            } else {

                // Split BED files exist from a previous run - simply create a channel for them
                Channel
                    .fromFilePairs("${outreg}/*.bed", size: 1)
                    .filter { it[0] != "DONE" }
                    .map { id, bed ->
                        def sp = id.split(':')
                        tuple(sp[0], sp[1], file(bed[0]))
                    }
                    .combine(ch_ref, by: 0)
                    .set { ch_regions_by_chrom }
            }

            // Reintroduce reference files that might not have associated region files
            ch_regions_by_chrom.groupTuple(by:0).set { temp }
            ch_ref_fa
                .join(temp, remainder: true)
                .map { 
                    if (!it[2]) {
                        return tuple(it[0], it[1], [null], [file("${outtmp}/no_regions")])
                    }
                    return it
                }
                .transpose()
                .set { ch_ref_fa_regions }

        // Standard - [ ref.basename, reference, bed ]
        } else {
            // Join with data channels - set as no_regions if no regions
            ch_ref_fa
                .join(ch_regions, by: 0, remainder: true)
                .map {
                    it[2] = it[2] ?: [ file("${outtmp}/no_regions") ]
                    return it
                }
                .unique()
                .set { ch_ref_fa_regions }
        }

    // No regions - empty channels to fill in the gaps
    } else {
 
        // Set no regions (null or file depending on method)
        if (params.method == 'joint') {
            ch_ref_fa
                .map {
                    tuple(it[0], it[1], null, file("${outtmp}/no_regions"))
                }
                .set { ch_ref_fa_regions }

        } else {
            ch_ref_fa
                .map {
                    tuple(it[0], it[1], [ file("${outtmp}/no_regions") ])
                }
                .set { ch_ref_fa_regions }
        }
    }

    /*
    Channel structures going into the next section:
        * Joint with regions - [ ref.bn, ref, chr, bed ]
        * Joint no regions - [ ref.bn, ref, null, no_regions ]
        * Standard with regions - [ ref.bn, ref, bed ]
        * Standard no regions - [ ref.bn, ref, no_regions ]
    */

    // Manage input channels depending on variant calling type
    if (params.method == 'joint') {

        // Check if bamlist files exist already - use existing to preserve 'resume' functionality
        def checkBamlist = Tools.checkDirAndFiles(outbml)

        // If Bam-list files don't exist, or force is provided
        if ( (!checkBamlist) || (force) ) {
            // Create output directory for bamlist files
            file(outbml).mkdirs()

            // Group bam files by their respective reference genomes
            ch_tmp
                .groupTuple(by: 1)
                .map { bams, ref ->
                    tuple(ref.baseName, bams)
                }
                .set { ch_id_bam }
            
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
                    file( "$outbml/DONE.bamlist") << '\n'
                }
            
            Channel
                .watchPath("${outbml}/*.bamlist", 'create,modify')
                .until { file -> file.baseName == 'DONE' }
                .map {
                    tuple(it.baseName, it)
                }
                .set { ch_bamlist }
        } else {
            Channel
                .fromFilePairs("${outbml}/*.bamlist", size: 1)
                .filter { it[0] != "DONE" }
                .set { ch_bamlist }
        }

        // Join reference + bams channel with bamlist and regions channel
        ch_tmp
            .groupTuple(by: 1)
            .map { bams, ref ->
                tuple(ref.baseName, bams)
            }
            .join(ch_bamlist)
            .combine(ch_ref_fa_regions, by: 0)
            .unique()
            .set { ch_data }

    // Standard mode - calling variants for individual samples
    } else {
        ch_tmp
            .map { bam, ref ->
                tuple(ref.baseName, bam)
            }
            .combine(ch_ref_fa_regions, by: 0)
            .map {
                tuple(it[1][0].baseName, it[1][0], it[1][1], it[2], it[3][0])
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
            caller,
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
            .set { ch_input }

        // Genotype and call variants - each channel input is a separate chromosome
        joint_mpileup_call(
            ch_input,
            caller,
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
                def vcfs = []
                def idxs = []
                it[1].each { i ->
                    vcfs << i[0]
                    idxs << i[1]
                }
                tuple(it[0], vcfs, idxs)
            }
            .set { ch_vcfs }
        
        // TODO: Only concat if regions == True
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