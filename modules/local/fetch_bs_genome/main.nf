process FETCH_BS_GENOME {
    tag "Version: ${params.bs_genome_version}, fetching from zenodo/4088020/files/genome_bs.zip"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ params.methylarray_deps_container }"

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
