import groovy.json.JsonBuilder

process SNIFFLES2 {
    tag "$meta.id"
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:shac591518dd32ecc3936666c95ff08f6d7474e9728"

    input:
        tuple val(meta), path(xam), path(xam_idx)
        path tr_bed
        path ref //, path(ref_idx), path(ref_cache), env(REF_PATH) 
        // val genome_build
    output:
        tuple val(meta), path("*.sniffles.vcf"), emit: vcf, optional: true
        path "${xam}.wf_sv.snf", emit: snf, optional: true
        path "versions.yml", emit: versions
    script:
        // if tr_arg is not provided and genome_build is set
        // automatically pick the relevant TR BED from the SV image
        def tr_arg = ""
        if (tr_bed.name != 'OPTIONAL_FILE'){
            tr_arg = "--tandem-repeats ${tr_bed}"
        }

        genome_build = 'hg38' // TODO: make this into params
        tr_arg = "--tandem-repeats \${WFSV_TRBED_PATH}/${genome_build}.trf.bed"

        // TODO: make these into params
        params.sniffles_args = false
        params.min_sv_length = false
        params.phased = false
        params.cluster_merge_pos = false

        def sniffles_args = params.sniffles_args ?: ''
        def min_sv_len = params.min_sv_length ? "--minsvlen ${params.min_sv_length}" : ""
        // Perform internal phasing only if snp not requested; otherwise, use joint phasing.
        def phase = params.phased ? "--phase" : ""
        def cluster_merge_pos = params.cluster_merge_pos ? params.cluster_merge_pos : -1
    """
    sniffles \
        --threads $task.cpus \
        --sample-id ${meta.id} \
        --output-rnames \
        ${min_sv_len} \
        --cluster-merge-pos $cluster_merge_pos \
        --input $xam \
        --reference $ref \
        --input-exclude-flags 2308 \
        --snf ${xam}.wf_sv.snf \
        $tr_arg \
        $sniffles_args \
        $phase \
        --vcf ${xam}.sniffles.vcf
    sed '/.:0:0:0:NULL/d' ${xam}.sniffles.vcf > tmp.vcf
    mv tmp.vcf ${xam}.sniffles.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles2: \$(echo \$(sniffles --version) | sed 's/ /,/')
    END_VERSIONS
    """
}

process filterCalls {
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:shac591518dd32ecc3936666c95ff08f6d7474e9728"

    input:
        tuple val(meta), path(vcf)
        tuple val(meta), path(mosdepth_summary)
        path target_bed
    output:
        tuple val(meta), path("*.filtered.vcf"), emit: vcf
    script:
    // Programmatically define chromosome codes.
    // note that we avoid interpolation (eg. "${chr}N") to ensure that values
    // are Strings and not GStringImpl, ensuring that .contains works.
    ArrayList chromosome_codes = []
    ArrayList chromosomes = [1..22] + ["X", "Y", "M", "MT"]
    for (N in chromosomes.flatten()){
        chromosome_codes += ["chr" + N, "" + N]
    }
    String ctgs = chromosome_codes.join(',')
    def ctgs_filter = "--contigs ${ctgs}"
    """
    # Filter contigs requre the input VCF to be compressed and indexed
    bcftools view -O z $vcf > input.vcf.gz && tabix -p vcf input.vcf.gz

    # Create filtering script
    get_filter_calls_command.py \
        --bcftools_threads $task.cpus \
        --target_bedfile $target_bed \
        --vcf input.vcf.gz \
        --depth_summary $mosdepth_summary \
        --min_read_support "auto" \
        --min_read_support_limit 2 \
        ${ctgs_filter} > filter.sh

    # Run filtering
    bash filter.sh > ${meta.id}.filtered.vcf
    """
}

// NOTE This is the last touch the VCF has as part of the workflow,
//  we'll rename it with its desired output name here
process sortVCF {
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:shac591518dd32ecc3936666c95ff08f6d7474e9728"

    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), path("${meta.id}.wf_sv.vcf.gz"), emit: vcf_gz
        tuple val(meta), path("${meta.id}.wf_sv.vcf.gz.tbi"), emit: vcf_tbi
    script:
    """
    bcftools sort -m 2G -T ./ -O z $vcf > ${meta.id}.wf_sv.vcf.gz
    tabix -p vcf ${meta.id}.wf_sv.vcf.gz
    """
}