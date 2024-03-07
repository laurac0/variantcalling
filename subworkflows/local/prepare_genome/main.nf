//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { TABIX_TABIX as TABIX_KNOWN_INDELS          } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DBSNP                 } from '../../../modules/nf-core/tabix/tabix/main'
include { GATK4_CREATESEQUENCEDICTIONARY             } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GET_CHROM_SIZES                            } from '../../../modules/local/get_chrom_sizes'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME    } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {
    take:
    fasta             // channel: [mandatory] fasta 
    dbsnp             // channel: [optional]  dbsnp
    known_indels      // channel: [optional]  known_indels

    main:
    versions = Channel.empty()
    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [ file1, file2 ] becomes [ [ meta1, file1 ], [ meta2, file2 ] ]
    // outputs are collected to maintain a single channel for relevant TBI files

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    SAMTOOLS_FAIDX_GENOME(fasta, [[],[]])
    GET_CHROM_SIZES( SAMTOOLS_FAIDX_GENOME.out.fai )

    TABIX_DBSNP(dbsnp.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_KNOWN_INDELS(known_indels.flatten().map{ it -> [ [ id:it.baseName ], it ] } )

    // Gather versions of all tools used

    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(TABIX_KNOWN_INDELS.out.versions)
    versions = versions.mix(TABIX_DBSNP.out.versions)
    versions = versions.mix(GET_CHROM_SIZES.out.versions)


    emit:
    dict                  = GATK4_CREATESEQUENCEDICTIONARY.out.dict                               // path: genome.fasta.dict
    known_indels_tbi      = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
    dbsnp_tbi             = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()               // path: dbsnb.vcf.gz.tbi
    genome_chrom_sizes    = GET_CHROM_SIZES.out.sizes.collect()

    versions // channel: [ versions.yml ]
}