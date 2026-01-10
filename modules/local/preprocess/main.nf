process PREPROCESS {
    tag "${idats.size()}"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ params.methylarray_deps_container }"

    input:
    path(idats)
    path(sample_mapping)

    output:
    tuple val("${idats.size()}"), path("*.csv")          , emit: csv
    tuple val("${idats.size()}"), path("mSetSqFlt.RData"), emit: rdata
    tuple val("${idats.size()}"), path("rgSet.RData")    , emit: rdata_rgSet
    path  "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "preprocess_1.R"
}
