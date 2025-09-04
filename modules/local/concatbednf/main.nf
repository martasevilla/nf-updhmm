process CONCAT_BED {
  tag { meta.id }
  publishDir "${params.outdir}/preprocess/sv_mask", mode: 'copy'

  input:
    tuple val(meta), path(beds)

  output:
    tuple val(meta), path("${meta.id}.sv.mask.sorted.bed")

  when:
    beds && beds.size() > 0

  shell:
  '''
  set -euo pipefail
  cat ${beds} | awk 'BEGIN{OFS="\\t"}{ if(NF>=3){print $1,$2,$3} }' \
    | sort -k1,1 -k2,2n | uniq > ${meta.id}.sv.mask.sorted.bed
  '''
}
