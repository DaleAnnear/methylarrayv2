/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// methylarray
include { PREPROCESS                    } from "../modules/local/preprocess/main"
include { FETCH_BS_GENOME               } from "../modules/local/fetch_bs_genome/main"
include { XREACTIVE_PROBES_FIND_REMOVE  } from "../modules/local/xreactive_probes_find_remove/main"
include { REMOVE_SNP_PROBES             } from "../modules/local/remove_snp_probes/main"
include { REMOVE_SEX_CHROMOSOMES        } from "../modules/local/remove_sex_chromosomes/main"
include { REMOVE_CONFOUNDING_PROBES     } from "../modules/local/remove_confounding_probes/main"
include { ADJUST_CELL_COMPOSITION       } from "../modules/local/adjust_cell_composition/main"
include { ADJUST_BATCH_EFFECT           } from "../modules/local/adjust_batch_effect/main"
include { FIND_DMP                      } from "../modules/local/find_dmp/main"
include { FIND_DMR                      } from "../modules/local/find_dmr/main"
include { FIND_BLOCKS                   } from "../modules/local/find_blocks/main"


// nf-core
include { MULTIQC                } from "../modules/nf-core/multiqc/main"
include { paramsSummaryMap       } from "plugin/nf-schema"
include { paramsSummaryMultiqc   } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { softwareVersionsToYAML } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { methodsDescriptionText } from "../subworkflows/local/utils_nfcore_methylarray_pipeline"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLARRAY {

    take:
    ch_input

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_preprocessed_files = channel.empty()
    extensive_metadata = params.sample_metadata ? channel.fromPath(params.sample_metadata) : channel.empty()

    //
    // MODULE: Run PREPROCESS
    //
    ch_input
        .map{ meta, red, green -> red.getName().replaceFirst(/_Red\.idat(?:\.gz)?$/,"") + "," + meta.id }
        .collectFile(storeDir: params.outdir, name: "sample_mapping.csv", sort: true, newLine: true)
        .set{ ch_sample_mapping }

    PREPROCESS (
        ch_input.map{meta, red, green -> [red, green]}.collect(),
        ch_sample_mapping
    )
    ch_versions = ch_versions.mix(PREPROCESS.out.versions.first())

    //
    // MODULE: Run XREACTIVE_PROBES_FIND_REMOVE
    //
    // Download from: https://github.com/pjhop/DNAmCrosshyb at https://doi.org/10.5281/zenodo.4088019 and point to the version
    // TODO: If bs_genome_path is not provided then the pipeline might resolve it itself
    if (params.remove_xreactive) {
        FETCH_BS_GENOME()
        XREACTIVE_PROBES_FIND_REMOVE (
            PREPROCESS.out.rdata,
            params.bs_genome_path ? file(params.bs_genome_path) : FETCH_BS_GENOME.out.bs_genome
        )
        ch_versions = ch_versions.mix(XREACTIVE_PROBES_FIND_REMOVE.out.versions.first())
    }

    //
    // MODULE: Run REMOVE_SNP_PROBES
    //
    if (params.remove_snp_probes) {
        REMOVE_SNP_PROBES (
            params.remove_xreactive ? XREACTIVE_PROBES_FIND_REMOVE.out.rdata : PREPROCESS.out.rdata
        )
        ch_versions = ch_versions.mix(REMOVE_SNP_PROBES.out.versions.first())
    }

    //
    // Optional steps of methylarray
    //

    //
    // Output channel following optional steps
    //
    final_bVals_ch = Channel.empty()
    current_bVals_ch = Channel.empty()

    if (params.run_optional_steps) {

        current_bVals_ch = params.remove_snp_probes ? REMOVE_SNP_PROBES.out.csv_bVals : XREACTIVE_PROBES_FIND_REMOVE.out.csv.filter { it == "bVals_noXprob.csv" }
        current_mVals_ch = params.remove_snp_probes ? REMOVE_SNP_PROBES.out.csv_mVals : XREACTIVE_PROBES_FIND_REMOVE.out.csv.filter { it == "mVals_noXprob.csv" }
        current_mSetSqFlt_ch = params.remove_snp_probes ? REMOVE_SNP_PROBES.out.rdata : XREACTIVE_PROBES_FIND_REMOVE.out.rdata.filter { it == "mSetSqFlt_noXprob.RData" }

        if (params.remove_sex_chromosomes) {
            current_bVals_ch = params.remove_snp_probes ? REMOVE_SNP_PROBES.out.rdata : XREACTIVE_PROBES_FIND_REMOVE.out.rdata
            //
            // MODULE: Run REMOVE_SEX_CHROMOSOMES
            //
            REMOVE_SEX_CHROMOSOMES (
                current_bVals_ch,
                PREPROCESS.out.rdata_rgSet
            )
            current_bVals_ch = REMOVE_SEX_CHROMOSOMES.out.bVals_csv
            current_mVals_ch = REMOVE_SEX_CHROMOSOMES.out.mVals_csv
            current_mSetSqFlt_ch = REMOVE_SEX_CHROMOSOMES.out.mSetSqFlt
            ch_versions = ch_versions.mix(REMOVE_SEX_CHROMOSOMES.out.versions.first())
        }

        if (params.remove_confounding_probes) {
            //
            // MODULE: Run REMOVE_CONFOUNDING_PROBES
            //
            REMOVE_CONFOUNDING_PROBES (
                current_mVals_ch,
                current_bVals_ch,
                current_mSetSqFlt_ch,
                extensive_metadata
            )
            current_bVals_ch = REMOVE_CONFOUNDING_PROBES.out.bVals
            ch_versions = ch_versions.mix(REMOVE_CONFOUNDING_PROBES.out.versions.first())
        }

        if (params.adjust_cell_composition) {
            //
            // MODULE: Run ADJUST_CELL_COMPOSITION
            //
            ADJUST_CELL_COMPOSITION (
                current_bVals_ch
            )
            current_bVals_ch = ADJUST_CELL_COMPOSITION.out.bVals
            ch_versions = ch_versions.mix(ADJUST_CELL_COMPOSITION.out.versions.first())
        }

        if (params.adjust_batch_effect) {
            //
            // MODULE: Run REMOVE_CONFOUNDING_PROBES
            //
            ADJUST_BATCH_EFFECT (
                current_bVals_ch,
                extensive_metadata
            )
            current_bVals_ch = ADJUST_BATCH_EFFECT.out.bVals
            ch_versions = ch_versions.mix(ADJUST_BATCH_EFFECT.out.versions.first())
        }
    }

    //
    // Update final bVals channel
    //
    final_bVals_ch = current_bVals_ch

    //
    // Generate group annotation
    //
    if (params.find_dmps || params.find_dmrs || params.find_blocks) {
        ch_input.map { meta, red, green -> "${meta.id},${meta.group}" }
            .collectFile(storeDir: params.outdir, name: "sample_annotation.csv", sort: true, newLine: true)
            .set{ ch_sample_annotation }
    }

    //
    // MODULE: Run FIND_DMP
    //
    if (params.find_dmps) {
        FIND_DMP (
            final_bVals_ch,
            ch_sample_annotation
        )
        ch_versions = ch_versions.mix(FIND_DMP.out.versions.first())
    }

    //
    // MODULE: Run FIND_DMR
    //
    if (params.find_dmrs) { // Will not be able to find DMRs with test data
        FIND_DMR (
            final_bVals_ch,
            ch_sample_annotation
        )
        ch_versions = ch_versions.mix(FIND_DMR.out.versions.first())
    }
    //
    // MODULE OPTIONAL: Run FIND_BLOCKS
    //
    if (params.find_blocks) {
        //
        // MODULE: Run FIND_BLOCKS
        //
        FIND_BLOCKS (
            final_bVals_ch,
            ch_sample_annotation
        )
        ch_versions = ch_versions.mix(FIND_BLOCKS.out.versions.first())
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: "nf_core_"  +  "methylarray_software_"  + "mqc_"  + "versions.yml",
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml"))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: "methods_description_mqc.yaml",
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
