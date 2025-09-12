//
// Check input samplesheet and get read, sample, and case channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
    canal_sp = SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' ) // canal donde cada elemento es una fila del csv

    canal_sp
        .map {create_bam_bai_bed_channel(it)}
        .set { reads_bam_bai_bed }

    canal_sp
        .map {create_bam_channel(it)}
        .set { reads_bam }
    
    canal_sp
        .map {create_bam_bai_channel(it)}
        .set { reads_bam_bai }

    canal_sp
        .map {create_bed_channel(it)}
        .set { reads_bed }

    //reads_bam = reads.flatten().first().concat(reads.flatten().filter(~/.*.bam/).toList()).toList()

    emit:
    reads_bam_bai_bed // channel: [ meta, [ bam, bai, bed] ]
    reads_bam // channel: [ meta, [ bam ] ]
    reads_bam_bai // channel: [ meta, [ bam, bai ] ]
    reads_bed // channel: [ meta, [ bed ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

}

//Function to get list of [ meta, [ bam, bai, bed] ]

def create_bam_bai_bed_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the bam/bai files to the meta map
    def bam_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
    }

    if (!file(row.bai).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam.bai file does not exist!\n${row.bai}"
    }

    if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bed file does not exist!\n${row.bed}"
    }

    bam_meta = [ meta, file(row.bam), file(row.bai), file(row.bed) ]
    //bam_meta = [ meta, [ file(row.bam)] ]

    return bam_meta
}

def create_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the bam/bai files to the meta map
    def bam_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
    }

    bam_meta = [ meta, file(row.bam) ]
    //bam_meta = [ meta, [ file(row.bam)] ]

    return bam_meta
}

def create_bam_bai_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the bam/bai files to the meta map
    def bam_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
    }

    if (!file(row.bai).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam.bai file does not exist!\n${row.bai}"
    }

    bam_bai_meta = [ meta, file(row.bam), file(row.bai) ]

    return bam_bai_meta
}

def create_bed_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the bam/bai files to the meta map

    if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bed file does not exist!\n${row.bed}"
    }

    bed_meta = [ meta, file(row.bed) ]
    //bam_meta = [ meta, [ file(row.bam)] ]

    return bed_meta
}
