//
// Filter structural variants using BED regions
//
// Steps:
//   1. Split samples into those with SVs and those without.
//   2. For samples with SVs: concatenate and merge BED files into a unified region mask.
//   3. Apply the mask to the input VCF using bcftools view.
//   4. Pass through unmodified VCFs for samples without SVs.
//   5. Emit combined results.
//

include { CONCAT_BED } from '../../../modules/local/concat_bed/main'
include { BEDTOOLS_MERGE } from '../../../modules/nf-core/bedtools/merge/main'
include { BCFTOOLS_VIEW as VIEW_MASK } from '../../../modules/nf-core/bcftools/view/main'

workflow SV_MASK_BED {
    
    take:
    sv_paths_ch    // channel: [meta, sv_p, sv_m, sv_f, vcf, tbi]
    
    main:
    
    ch_versions = Channel.empty()

    // Split samples depending on whether they have SVs or not
    sv_paths_ch.branch { meta, sv_p, sv_m, sv_f, vcf, tbi ->
        with_sv: [sv_p, sv_m, sv_f].any { it != '-' }
        no_sv: [sv_p, sv_m, sv_f].every { it == '-' }
    }.set { sv_branches }

    // ------------------
    // Process samples with SVs
    // ------------------

    beds_ch = sv_branches.with_sv.map { meta, sv_p, sv_m, sv_f, vcf, tbi ->
        def beds = [sv_p, sv_m, sv_f].findAll { it != '-' }.collect { file(it) }
        tuple(meta, beds, vcf, tbi)
    }

    // Concatenate multiple BEDs into one file per trio
    CONCAT_BED(beds_ch.map { meta, beds, vcf, tbi -> tuple(meta, beds) })
    concat_bed_ch = CONCAT_BED.out  
    
    // Merge overlapping regions to create a unified mask
    BEDTOOLS_MERGE(concat_bed_ch)
    merged_bed_ch = BEDTOOLS_MERGE.out.bed
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

    mask_combined_ch = beds_ch
        .map { meta, beds, vcf, tbi -> tuple(meta, vcf, tbi) }      
        .join(merged_bed_ch.map { meta, bed -> tuple(meta, bed) })
        .map { meta, vcf, tbi, bed ->
            tuple(meta, vcf, tbi, bed)
        }

    vcf_input_ch   = mask_combined_ch.map { meta, vcf, tbi, bed -> tuple(meta, vcf, tbi) }
    regions_only_ch = mask_combined_ch.map { meta, vcf, tbi, bed -> bed }

    empty_ch = Channel.value([])

    // Apply the BED mask to the VCF using bcftools view
    VIEW_MASK(vcf_input_ch, regions_only_ch, empty_ch, empty_ch) 

    with_sv_vcf_ch = VIEW_MASK.out.vcf
    with_sv_tbi_ch = VIEW_MASK.out.tbi
    with_sv_joined_ch = with_sv_vcf_ch.join(with_sv_tbi_ch)
        .map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }
    ch_versions = ch_versions.mix(VIEW_MASK.out.versions)
    
    // ------------------
    // Process samples without SVs 
    // ------------------
    no_sv_vcf_ch = sv_branches.no_sv.map { meta, sv_p, sv_m, sv_f, vcf, tbi -> tuple(meta, vcf, tbi) }

    // Combine both channels
    all_vcfs_ch = with_sv_joined_ch.mix(no_sv_vcf_ch)
    
    emit:
    vcfs     = all_vcfs_ch           // channel: [meta, vcf, tbi]
    versions = ch_versions           // channel: versions.yml
}