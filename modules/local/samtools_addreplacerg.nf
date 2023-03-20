process SAMTOOLS_ADDREPLACERG {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(input), path(index)

    output:
    tuple val(meta), path("*.{cram,bam}"), emit: bam
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args  ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    """
    samtools addreplacerg -w -r  ${meta.read_group} $input -o ${prefix}.${file_type}
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
