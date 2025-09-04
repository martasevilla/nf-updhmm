include { CHECK_SAMPLESHEET } from '../../modules/local/check_samplesheet'   // nf-core local default (one row per family)
include { BCFTOOLS_MERGE    } from '../../modules/nf-core/bcftools/merge/main'
include { SV_MASK_BED       } from './sv_mask_bed'
include { FILTER_LOWCONF    } from './filter_lowconf'

workflow PREPROCESS_VCF {
  take:
    // Channel of Maps (one per family) from the CSV samplesheet
    ROWS

  main:
    // 1) Samplesheet check → ensure required files/columns exist
    ROWS | CHECK_SAMPLESHEET

    // 2) Prepare meta (carry SV BED paths early) and SNV VCF list → MERGE
    CHECK_SAMPLESHEET.out
      .map { row ->
        def meta = [
          id  : row.fam_id,
          sv_f: row.path_sv_father  ?: '-',
          sv_m: row.path_sv_mother  ?: '-',
          sv_p: row.path_sv_proband ?: '-'
        ]
        def vcfs = [
          file(row.path_vcf_father),
          file(row.path_vcf_mother),
          file(row.path_vcf_proband)
        ]
        tuple(meta, vcfs)
      } \
      | BCFTOOLS_MERGE

    // 3) Build BED list from meta, then apply SV mask (if any BED present)
    BCFTOOLS_MERGE.out
      .map { meta, merged_vcf ->
        def beds = []
        if (meta.sv_f && meta.sv_f != '-') beds << file(meta.sv_f)
        if (meta.sv_m && meta.sv_m != '-') beds << file(meta.sv_m)
        if (meta.sv_p && meta.sv_p != '-') beds << file(meta.sv_p)
        tuple(meta, merged_vcf, beds)
      } \
      | SV_MASK_BED

    // 4) Low-confidence filters (biallelic, trio ref-homo, GQ/DP, centromeres, segdups, HLA/KIR)
    SV_MASK_BED.masked_vcf | FILTER_LOWCONF

  emit:
    // tuple(meta, vcf_path_indexed)
    preprocessed_vcf = FILTER_LOWCONF.preprocessed_vcf
}
