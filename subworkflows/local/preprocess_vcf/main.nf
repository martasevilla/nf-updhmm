//
// Preprocess VCF files: validate samplesheet, remove annotations, combine VCFs, filter SVs, and apply quality filters
//

include { SAMPLESHEET_CHECK } from '../../../modules/local/samplesheet_check'
include { REMOVE_ANNOTATIONS } from '../../../subworkflows/local/remove_annotations'
include { COMBINE_VCF } from '../../../subworkflows/local/combine_vcf'
include { SV_MASK_BED } from '../../../subworkflows/local/sv_mask_bed/main'
include { FILTER_LOWCONF } from '../../../subworkflows/local/filter_lowconf/main'

workflow PREPROCESS_VCF {
    
    take:
    samplesheet_path    // path: samplesheet CSV file
    
    main:
    
    // Create channel from samplesheet
    samplesheet_ch = Channel.fromPath(samplesheet_path, checkIfExists: true)

    // Validate samplesheet
    SAMPLESHEET_CHECK(samplesheet_ch)
    validated_samplesheet = SAMPLESHEET_CHECK.out.csv

    // Build channel of VCFs and metadata
    samples_ch = validated_samplesheet
        .splitCsv(header:true)
        .map { row ->
            def meta = [
                id  : row.fam_id,
                sv_p: row.path_sv_proband ?: '-',
                sv_m: row.path_sv_mother  ?: '-',
                sv_f: row.path_sv_father  ?: '-'
            ]
            def vcfs = [
                file(row.path_vcf_proband),
                file(row.path_vcf_mother),
                file(row.path_vcf_father)
            ]
            tuple(meta, vcfs)
        }

    // Step 1: Remove annotations (Keep only essential fields)
    REMOVE_ANNOTATIONS(samples_ch)
    
    // Step 2: Intersection and merge (Combine the three VCFs into one)
    COMBINE_VCF(REMOVE_ANNOTATIONS.out.vcfs)

    // Step 3: Filter Structural Variants
    SV_MASK_BED(COMBINE_VCF.out.sv_paths)
    
    // Step 4: Apply low confidence filters
    FILTER_LOWCONF(SV_MASK_BED.out.vcfs)
    
    emit:
    vcfs     = FILTER_LOWCONF.out.vcfs    // channel: [meta, vcf, tbi] - Final processed VCFs
}