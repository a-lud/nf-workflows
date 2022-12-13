/*
PSMC pipeline
    * 1. Generate consensus sequence
        - BAM input -> use SAMtools consensus
        - VCF input -> use BCFtools/VCFutils (older but traditional method)
    * 2. Convert consensus FASTQ to PSMC-FASTA format
    * 3. Run PSMC
    * 4. Optionally - Run bootstrap 
        - Combine results
    * Generate plots?
*/

// Import pipeline functions
include { samtools_consensus } from '../nf-modules/samtools/1.16.1/samtools_consensus'
include { bcftools_vcfutils } from '../nf-modules/bcftools/1.15.1/bcftools_vcfutils'
include { fq2psmcfa } from '../nf-modules/psmc/0.6.5/fq2psmcfa'
include { psmc } from '../nf-modules/psmc/0.6.5/psmc'
include { splitfa } from '../nf-modules/psmc/0.6.5/splitfa'
include { psmc_bootstrap } from '../nf-modules/psmc/0.6.5/psmc_bootstrap'
include { psmc_write } from '../nf-modules/psmc/0.6.5/psmc_write'

// Sub-workflow
workflow PSMC {
    main:

    // Output directories
    def outdir = [params.outdir, params.out_prefix].join('/')
    def outdir_cov = [ outdir, "coverage" ].join("/")
    def outdir_psmc = [ outdir, "psmc" ].join("/")

    // Read in files (BAM or VCF)
    Channel
        .fromFilePairs(
            [params.filedir.path, params.filedir.pattern].join('/'),
            size: params.filedir.nfiles,
        )
        .ifEmpty { exit 1, "Can't find read files." }
        .set { ch_files }

    // Split clocks`
    Channel
        .fromList(params.clock.tokenize(' '))
        .set { ch_clocks }

    // Parse regions file
    if (params.containsKey('regions')) {
        Channel
            .fromPath(params.regions)
            .splitCsv()
            .map { row -> 
                id = row[0]
                bed = file(row[1])

                tuple(id, bed)
            }
            .set { ch_regions }
    }
    
    // Process data for type of pipeline
    if (params.filetype == "bam") {
        // Assume BAI exists in BAM directory
        ch_files
            .map { tuple(it[0], it[1], file([it[1][0], "bai"].join('.'))) }
            .set { ch_temp }

        // Join with regions channel
        ch_bam = params.containsKey("regions") ? 
            ch_temp
                .join(ch_regions, remainder: true)
                .map { 
                    it[3] = it[3] ?: file("EMPTY")
                    return it
                } : 
            ch_temp.combine([file("EMPTY")])

        // Double check the join has worked - second field will be null if CSV key doesn't match
        // a file.
        ch_bam.map {
            assert it[1] != null : "ERROR: ${it[0]} in CSV doesn't have a matching key."
        }

        // Create consensus FASTQ
        samtools_consensus(ch_bam, params.minMQ)
        samtools_consensus.out.set { ch_consensus }
    } else {
        // Join VCF channel with regions
        ch_vcf = params.containsKey("regions") ? 
            ch_files
                .join(ch_regions, remainder: true)
                .map { 
                    it[2] = it[2] ?: file("EMPTY")
                    return it
                 } :
                ch_files.combine([file("EMPTY")])
        
        // Double check the join has worked - second field will be null if CSV key doesn't match
        // a file.
        ch_vcf.map {
            assert it[1] : "ERROR: ${it[0]} in CSV doesn't have a matching key."
        }

        // Consensus sequence from VCF/BAM
        bcftools_vcfutils(ch_vcf)
        bcftools_vcfutils.out.set { ch_consensus }
    }

    // Convert VCF-FASTQ to PSMCFA
    fq2psmcfa(ch_consensus)

    // Make PSMCFA-clock combinations
    fq2psmcfa.out.combine(ch_clocks).set { ch_psmc_in }

    // Run PSMC
    psmc(ch_psmc_in)

    // Bootstrap
    if (params.bootstrap) {
        // Split fasta into chunks
        splitfa(fq2psmcfa.out)
        splitfa.out.combine(ch_clocks).set { ch_splitfa_clocks } 

        // Run PSMC bootstrap
        psmc_bootstrap(ch_splitfa_clocks)

        // Combine bootstrap output with original run
        psmc.out
            .join(psmc_bootstrap.out, by: [0,1])
            .map {
                // [ id, clock, [psmc files] ] - first file in 'psmc files' is the full run (-0.psmc)
                return tuple(
                    it[0],
                    it[1],
                    [ it[2] ] + it[3]
                )
            }
            .set { ch_psmc }
    } else {
        psmc.out.set { ch_psmc }
    }

    // Write psmc outputs to file
    psmc_write(ch_psmc, outdir_psmc)
}