process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bcftools=1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1':
        'quay.io/biocontainers/bcftools:1.16--hfe4b78e_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path regions
    path targets
    path samples

    output:
    tuple val(meta), path("*stats.txt"), emit: stats
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_file = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file =  samples ? "--samples-file ${samples}" : ""
    def is_clairs = meta.variantcaller == "clairs" ? "bcftools +fill-tags $vcf -Ob -o ${vcf.baseName}_tagged.vcf.gz -- -t all,\"INFO/DP:1=int(sum(FORMAT/DP))\"\ncat ${vcf.baseName}_tagged.vcf.gz" : ""
    def vcf_filename = meta.variantcaller == "clairs" ? "${vcf.baseName}_tagged.vcf.gz" : "${vcf.baseName}.gz"
    """
    ${is_clairs}
    bcftools stats \\
        --verbose -s - \\
        $args \\
        $regions_file \\
        $targets_file \\
        $samples_file \\
        ${vcf_filename} > ${prefix}.bcftools_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}