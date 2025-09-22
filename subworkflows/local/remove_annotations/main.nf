//
// Remove annotations from individual VCF files and regroup by family
//

include { BCFTOOLS_ANNOTATE as BCFTOOLS_DELETE_ANNOTATIONS } from '../../../modules/nf-core/bcftools/annotate/main'

workflow REMOVE_ANNOTATIONS {
    
    take:
    samples_ch    // channel: [meta, vcfs] 
    
    main:
    
    // Flatten to process each VCF individually
    individual_vcfs_ch = samples_ch.flatMap { meta, vcfs ->  
        def individuals = []
        def roles = ['proband', 'mother', 'father']
        
        vcfs.eachWithIndex { vcf, idx ->
            def individual_meta = meta.clone()
            individual_meta.role = roles[idx]
            individual_meta.original_family_id = meta.id
            individual_meta.id = "${meta.id}_${roles[idx]}"
            
            individuals << tuple(individual_meta, vcf, [], [], [])
        }
        return individuals
    }

    // Apply delete annotations to each individual VCF
    BCFTOOLS_DELETE_ANNOTATIONS(individual_vcfs_ch, [], [], [])
    
    // Regroup the processed VCFs back into family trios
    regrouped_ch = BCFTOOLS_DELETE_ANNOTATIONS.out.vcf
        .join(BCFTOOLS_DELETE_ANNOTATIONS.out.tbi)
        .map { meta, vcf, tbi ->
            tuple(meta.original_family_id, meta.role, meta, vcf, tbi)
        }
        .groupTuple(by: 0)  // Group by original family ID
        .map { family_id, roles, metas, vcfs, tbis ->
            // Reconstruct the original metadata
            def family_meta = [
                id: family_id,
                sv_p: metas[0].sv_p,
                sv_m: metas[0].sv_m,
                sv_f: metas[0].sv_f  
            ]
            
            // Sort VCFs and indices by role to maintain order (proband, mother, father)
            def role_order = ['proband', 'mother', 'father']
            def sorted_data = [roles, vcfs, tbis].transpose().sort { 
                it[0] in role_order ? role_order.indexOf(it[0]) : 999 
            }
            def sorted_vcfs = sorted_data.collect { it[1] }
            def sorted_tbis = sorted_data.collect { it[2] }
            
            tuple(family_meta, sorted_vcfs, sorted_tbis)
        }
    
    emit:
    vcfs     = regrouped_ch           // channel: [meta, vcfs, tbis]
}