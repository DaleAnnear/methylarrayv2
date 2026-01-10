process FIND_DMP {
    tag "${samplesheet_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ params.methylarray_deps_container }"

    input:
    tuple val(samplesheet_name), path(bVALS_SNPPROBES)
    path(extensive_metadata)

    output:
    tuple val(samplesheet_name), path("dmp_champ.*.csv"), emit: all
    tuple val(samplesheet_name), path("dmp_minfi.csv")  , emit: minfi
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'find_dmp_4.R'
}
