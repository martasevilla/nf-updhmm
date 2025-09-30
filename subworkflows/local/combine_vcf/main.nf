//
// Perform VCF intersection and merge operations
//

include { BCFTOOLS_ISEC  } from '../../../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_MERGE } from '../../../modules/nf-core/bcftools/merge/main'

workflow COMBINE_VCF {
    
    take:
    vcfs_ch    // channel: [meta, vcfs, tbis]
    
    main:
    ch_versions = Channel.empty()

    // VCF intersection
    BCFTOOLS_ISEC(vcfs_ch)
    ch_versions = ch_versions.mix(BCFTOOLS_ISEC.out.versions)

    merge_ch = BCFTOOLS_ISEC.out.results.map { meta, dir ->
        def vcfs = file("${dir}/000*.vcf.gz")
        def tbis = file("${dir}/*.tbi")
        tuple(meta, vcfs, tbis)
    }

    empty_tuple_ch = Channel.value([[:], []])

    // Merge intersected VCFs
    BCFTOOLS_MERGE(merge_ch, empty_tuple_ch, empty_tuple_ch, empty_tuple_ch)
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    merged_with_sv_paths = BCFTOOLS_MERGE.out.vcf
        .join(BCFTOOLS_MERGE.out.index)
        .map { meta, vcf, tbi ->
            tuple(meta, meta.sv_p, meta.sv_m, meta.sv_f, vcf, tbi)
        }
    
    emit:
    merged_vcf   = BCFTOOLS_MERGE.out.vcf      // channel: [meta, vcf]
    merged_index = BCFTOOLS_MERGE.out.index    // channel: [meta, tbi]
    sv_paths     = merged_with_sv_paths        // channel: [meta, sv_p, sv_m, sv_f, vcf, tbi]
    versions     = ch_versions                 // channel: versions.yml

}