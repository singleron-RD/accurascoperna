/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_accurascoperna_pipeline'

process CONVERT {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::pyfastx=2.1.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyfastx:2.1.0--py39h3d4b85c_0' :
        'biocontainers/pyfastx:2.1.0--py39h3d4b85c_0' }"

    input:
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    tuple val(meta), path(reads)
    path star_genome
    path assets_dir
    val protocol
    val starsolo_common_args
    val starsolo_p3_args
    val starsolo_p3p5_args

    output:
    tuple val(meta), path("${meta.id}_p3_*.fastq"), emit: p3_reads
    tuple val(meta), path("${meta.id}_p5_*.fastq"), emit: p5_reads
    path "${meta.id}.p3.starsolo_cmd.txt" , emit: p3_starsolo_cmd
    path "${meta.id}.p3p5.starsolo_cmd.txt" , emit: p3p5_starsolo_cmd
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    convert.py \\
        --sample ${meta.id} \\
        --genomeDir ${star_genome} \\
        --fq1 ${forward.join( "," )} \\
        --fq2 ${reverse.join( "," )} \\
        --assets_dir ${assets_dir} \\
        --protocol ${protocol} \\
        --thread $task.cpus \\
        --common_args \"${starsolo_common_args}\"  \\
        --p3_args \"${starsolo_p3_args}\"  \\
        --p3p5_args \"${starsolo_p3p5_args}\"  

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyfastx: \$(pyfastx --version | sed -e "s/pyfastx version //g")
    END_VERSIONS
    """
}

process FILTER_GTF {
    tag "$gtf"
    label 'process_single'

    conda 'conda-forge::python==3.12'
    container "biocontainers/python:3.12"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    path gtf
    val attributes

    output:
    path "*.filtered.gtf", emit: filtered_gtf
    path "gtf_filter.log", emit: log_file

    script:
    def args = task.ext.args ?: ''

    """
    filter_gtf.py ${gtf} \"${attributes}\"
    """
}

process STAR_GENOME {
    tag "$genome_name"
    label 'process_high'

    conda "bioconda::star==2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    path fasta
    path gtf
    val genome_name

    output:
    path "$genome_name" , emit: star_genome
    path "versions.yml"            , emit: versions

    script:
    def args        = task.ext.args ?: ''
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def include_gtf = gtf ? "--sjdbGTFfile $gtf" : ''
    def fasta_sa = ( Math.log(fasta.size()) / Math.log(2) ) / 2 - 1
    def sa = Math.floor( Math.min(14, fasta_sa) )
    """
    mkdir ${genome_name}
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir ${genome_name}/ \\
        --genomeFastaFiles $fasta \\
        $include_gtf \\
        --runThreadN $task.cpus \\
        --genomeSAindexNbases ${sa} \\
        $memory \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}

process STARSOLO_P3 {
    tag "${meta.id}.p3"
    label 'process_high'

    conda "bioconda::star==2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(reads)
    path star_genome
    path assets_dir
    val starsolo_cmd

    output:
    tuple val(meta), path('*d.out.bam')               , emit: bam
    tuple val(meta), path('*.Solo.out/*')             , emit: solo_out
    path('*Log.final.out')                            , emit: log_mapping
    path  "versions.yml"                              , emit: versions
    path "*.Solo.out/GeneFull_Ex50pAS/${meta.id}.p3.Summary.csv" , emit: summary

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*out.mate')               , optional:true, emit: unmap
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab


    script:
    def prefix = "${meta.id}.p3."

    // do not indent
"""
${starsolo_cmd}

if [ -d ${prefix}Solo.out ]; then
    # Backslashes still need to be escaped (https://github.com/nextflow-io/nextflow/issues/67)
    find ${prefix}Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip -f {} \\;
fi

mv ${prefix}Solo.out/GeneFull_Ex50pAS/Summary.csv ${prefix}Solo.out/GeneFull_Ex50pAS/${prefix}Summary.csv

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    star: \$(STAR --version | sed -e "s/STAR_//g")
END_VERSIONS
"""
}

process STARSOLO_P3P5 {
    tag "${meta.id}.p3p5"
    label 'process_high'

    conda "bioconda::star==2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(p3_reads)
    tuple val(meta), path(p5_reads)
    path star_genome
    path assets_dir
    val starsolo_cmd

    output:
    tuple val(meta), path('*d.out.bam')               , emit: bam
    tuple val(meta), path('*.Solo.out/*')             , emit: solo_out
    path('*Log.final.out')                            , emit: log_mapping
    path  "versions.yml"                              , emit: versions

    path "*.Solo.out/GeneFull_Ex50pAS/${meta.id}.Summary.csv" , optional:true, emit: summary
    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*out.mate')               , optional:true, emit: unmap
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab


    script:
    def prefix = "${meta.id}"

    // do not indent
"""
${starsolo_cmd}

if [ -d ${prefix}.Solo.out ]; then
    # Backslashes still need to be escaped (https://github.com/nextflow-io/nextflow/issues/67)
    find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip -f {} \\;
fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    star: \$(STAR --version | sed -e "s/STAR_//g")
END_VERSIONS
"""
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ACCURASCOPERNA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC

    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // STAR genome
    def ch_star_genome = null
    if (params.star_genome) {
        ch_star_genome = params.star_genome
    } else {
        FILTER_GTF(
            params.gtf,
            params.keep_attributes,
        )
        STAR_GENOME(
            params.fasta,
            FILTER_GTF.out.filtered_gtf,
            params.genome_name
        )
        ch_versions = ch_versions.mix(STAR_GENOME.out.versions.first())
        ch_star_genome = STAR_GENOME.out.star_genome
    }

    // convert p5 barcode to p3
    CONVERT (
        ch_samplesheet,
        ch_star_genome,
        "${projectDir}/assets/",
        params.protocol,
        params.starsolo_common_args,
        params.starsolo_p3_args,
        params.starsolo_p3p5_args
    )

    // starsolo p3
    STARSOLO_P3 (
        CONVERT.out.p3_reads,
        ch_star_genome,
        "${projectDir}/assets/",
        CONVERT.out.p3_starsolo_cmd.map { it.text }
    )
    ch_multiqc_files = ch_multiqc_files.mix(STARSOLO_P3.out.log_mapping.collect()).mix(STARSOLO_P3.out.summary.collect())
    ch_versions = ch_versions.mix(STARSOLO_P3.out.versions.first())

    // starsolo p3p5
    STARSOLO_P3P5 (
        CONVERT.out.p3_reads,
        CONVERT.out.p5_reads,
        ch_star_genome,
        "${projectDir}/assets/",
        CONVERT.out.p3p5_starsolo_cmd.map { it.text }
    )


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
