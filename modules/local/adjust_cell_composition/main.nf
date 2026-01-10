process ADJUST_CELL_COMPOSITION {
    tag "${samplesheet_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ params.methylarray_deps_container }"

    input:
    tuple val(samplesheet_name), path(bVALS_SNPPROBES)

    output:
    tuple val(samplesheet_name), path("cmVals.cell_comp.csv"), emit: mVals
    tuple val(samplesheet_name), path("cbVals.cell_comp.csv"), emit: bVals
    path  "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'adjust_cell_composition_o3.R'
}
