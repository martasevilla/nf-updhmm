/*
 * FILTER_LOWCONF
 * Purpose: apply hard filters to produce the final preprocessed VCF per family.
 * Filters:
 *   1) Keep biallelic sites only
 *   2) Exclude sites where all three samples are reference-homozygous (0/0 or 0|0)
 *   3) Enforce MIN(FMT/GQ) and MIN(FMT/DP) across the trio
 *   4) Exclude centromeres
 *   5) Exclude segmental duplications
 *   6) Exclude HLA/KIR regions
 *
 * All bcftools arguments are configured in conf/modules.config via withName selectors.
 */

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_BIALLELIC         } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_REFHOMO_EXCLUDE    } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_QUAL_MIN           } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_EXCL_CENTROMERES   } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_EXCL_SEGDUPS       } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_EXCL_HLAKIR        } from '../../modules/nf-core/bcftools/view/main'
include { TABIX_INDEX                                      } from '../../modules/nf-core/tabix/tabix/main'

workflow FILTER_LOWCONF {
  take:
    // tuple(meta, vcf_path)
    MASKED_VCF

  main:
    MASKED_VCF \
      | BCFTOOLS_VIEW_BIALLELIC \
      | BCFTOOLS_VIEW_REFHOMO_EXCLUDE \
      | BCFTOOLS_VIEW_QUAL_MIN \
      | BCFTOOLS_VIEW_EXCL_CENTROMERES \
      | BCFTOOLS_VIEW_EXCL_SEGDUPS \
      | BCFTOOLS_VIEW_EXCL_HLAKIR \
      | TABIX_INDEX

  emit:
    // tuple(meta, vcf_path_indexed)
    preprocessed_vcf = TABIX_INDEX.out
}