//
// This process validates the provided sample sheet to ensure it follows the expected 
// format and contains all required fields. It runs the custom script 
// `check_samplesheet.py` to perform the validation. 
//

process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    //conda "conda-forge::python=3.8.3"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path samplesheet

    output:
    path 'samplesheet.valid.csv', emit: csv
    path "versions.yml"         , emit: versions

    script: // nf-core/chipseq/bin/
    """
    check_samplesheet.py \
        $samplesheet \
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
