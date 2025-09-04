include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

include { CONCAT_BED } from '../../modules/local/concat_bed'
include { BCFTOOLS_VIEW as VIEW_MASK } from '../../modules/nf-core/bcftools/view/main'
include { TABIX_INDEX } from '../../modules/nf-core/htslib/tabix/main'

workflow SV_MASK_BED {
  take:
    // Tuple: [ meta(Map: [id: fam_id]), merged_vcf(Path), sv_beds(List<Path or empty]) ]
    MERGED_VCF_SNV_BEDS

  main:
    MERGED_VCF_SNV_BEDS
      .branch { t -> (t[2] && t[2].size()>0) ? 'with_sv' : 'no_sv' }

    // With SVs: concat → bcftools view -T ^mask
    MERGED_VCF_SNV_BEDS.with_sv
      .map { meta, vcf, beds -> tuple(meta, beds) }
      | CONCAT_BED
      .map { meta, mask_bed -> tuple(meta, mask_bed) }
      .combine( MERGED_VCF_SNV_BEDS.with_sv.map{ meta, vcf, beds -> tuple(meta, vcf) } )
      .map { (m1, bed), (m2, vcf) -> tuple(m1, vcf, bed) }
      | VIEW_MASK( extra_args: ["-T", "^"], bed_file: true )
      | TABIX_INDEX

    // No SVs: pass-through + index
    MERGED_VCF_SNV_BEDS.no_sv
      | TABIX_INDEX

  emit:
    // both branches should emit: tuple(meta, vcf_path)
    masked_vcf = merge( VIEW_MASK.out, TABIX_INDEX.out )
}
