/*
 * FILTER_LOWCONF
 * Purpose: apply hard filters to produce the final preprocessed VCF per family.
 * Filters:
 *   1) Keep biallelic sites only
 *   2) Exclude sites where all three samples are reference-homozygous (0/0 or 0|0)
 *   3) Enforce MIN(FMT/GQ) and MIN(FMT/DP) across the trio
 *   4) Exclude centromeres, segmental duplications and HLA/KIR regions
 *
 * All bcftools arguments are configured in conf/modules.config via withName selectors.
 */

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_BIALLELIC         } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_REFHOMO_EXCLUDE    } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_QUAL_MIN           } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_EXCL_ALL   } from '../../../modules/nf-core/bcftools/view/main'

workflow FILTER_LOWCONF {
    
    take:
    vcfs_ch    // channel: [meta, vcf, tbi]
    
    main:
    ch_versions = Channel.empty()
    
    // Create empty channel for unused inputs
    empty_ch = Channel.value([])
    
    // Helper function to combine VCF and index
    def combineVcfIndex = { vcf_ch, index_ch ->
        return vcf_ch.join(index_ch).map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }
    }

    // Filter 1: Keep only biallelic variants
    BCFTOOLS_VIEW_BIALLELIC(vcfs_ch, empty_ch, empty_ch, empty_ch)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_BIALLELIC.out.versions)
    
    // Filter 2: Apply quality filters (GQ and DP)
    BCFTOOLS_VIEW_QUAL_MIN(combineVcfIndex(BCFTOOLS_VIEW_BIALLELIC.out.vcf, BCFTOOLS_VIEW_BIALLELIC.out.tbi), empty_ch, empty_ch, empty_ch)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_QUAL_MIN.out.versions)
    
    // Filter 3: Exclude reference homozygous variants
    BCFTOOLS_VIEW_REFHOMO_EXCLUDE(combineVcfIndex(BCFTOOLS_VIEW_QUAL_MIN.out.vcf, BCFTOOLS_VIEW_QUAL_MIN.out.tbi), empty_ch, empty_ch, empty_ch)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_REFHOMO_EXCLUDE.out.versions)
    
    // Filter 4: Exclude problematic genomic regions (centromeres, segmental duplications and HLA/KIR regions)
    BCFTOOLS_VIEW_EXCL_ALL(combineVcfIndex(BCFTOOLS_VIEW_REFHOMO_EXCLUDE.out.vcf, BCFTOOLS_VIEW_REFHOMO_EXCLUDE.out.tbi), empty_ch, empty_ch, empty_ch)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_EXCL_ALL.out.versions)
    
    // Final output
    final_vcfs = combineVcfIndex(BCFTOOLS_VIEW_EXCL_ALL.out.vcf, BCFTOOLS_VIEW_EXCL_ALL.out.tbi)
    
    emit:
    vcfs     = final_vcfs            // channel: [meta, vcf, tbi]
    versions = ch_versions           // channel: versions.yml
}