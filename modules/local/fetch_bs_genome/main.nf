process FETCH_BS_GENOME {
    tag "${params.bs_genome_path}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ params.methylarray_deps_container }"

    input:
    tuple val(samplesheet_name), path(RData_PREPROCESSING)
    path(genome_path)

    output:
    path("genome_bs/${params.bs_genome_version}"), emit: bs_genome

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget https://zenodo.org/records/4088020/files/genome_bs.zip?download=1 -O genome_bs.zip
    unzip genome_bs.zip
    """

}
