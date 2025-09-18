//
// This process takes multiple BED files containing structural variant (SV) regions 
// from different sources (e.g. proband, mother, father), concatenates them into a 
// single file, keeps only valid BED fields (chromosome, start, end), sorts the entries 
// by genomic coordinates, removes duplicates, and produces a clean, unified BED file.
//

process CONCAT_BED {
  tag { meta.id }
  publishDir "${params.outdir}/PREPROCESS/SV_MASK_BED", mode: 'copy'

  input:
    tuple val(meta), path(beds)

  output:
    tuple val(meta), path("${meta.id}_svs_concatenated.bed")

  when:
    beds && beds.size() > 0

  script: 
  def bed_files = (beds instanceof List) ? beds.join(' ') : beds
  """
  set -euo pipefail
  cat ${bed_files} \
    | awk 'BEGIN{OFS="\\t"}{ if(NF>=3){print \$1,\$2,\$3} }' \
    | sort -k1,1 -k2,2n | uniq > ${meta.id}_svs_concatenated.bed
  """
}
