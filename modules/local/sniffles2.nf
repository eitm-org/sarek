import groovy.json.JsonBuilder

process SNIFFLES2 {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::sniffles" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:2.5.3--pyhdfd78af_0':
        'quay.io/biocontainers/sniffles:2.5.3--pyhdfd78af_0' }"

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
        genome_build = 'hg38'
        if (tr_bed.name != 'OPTIONAL_FILE'){
            tr_arg = "--tandem-repeats ${tr_bed}"
        } else {
            tr_arg = "--tandem-repeats \${WFSV_TRBED_PATH}/${genome_build}.trf.bed"
        }

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
        --snf ${xam}.wf_sv.snf \
        $tr_arg \
        --mosaic \
        --mosaic-include-germline \
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
        path filter_bed
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
    # bedtools intersect -a $target_bed -b $filter_bed -v > filtered_target.bed

    # Filter contigs requre the input VCF to be compressed and indexed
    bcftools view -O z $vcf > input.vcf.gz && tabix -p vcf input.vcf.gz

    # Create filtering script
    get_filter_calls_command.py \
        --bcftools_threads $task.cpus \
        --target_bedfile $target_bed \ # filtered_target.bed
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

process getVersions {
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:shac591518dd32ecc3936666c95ff08f6d7474e9728"
    output:
        path "versions.txt"
    script:
    """
    trap '' PIPE # suppress SIGPIPE without interfering with pipefail
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    truvari version | sed 's/ /,/' >> versions.txt
    sniffles --version | head -n 1 | sed 's/ Version //' >> versions.txt
    bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    echo `seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d ':' -f 2 | sed 's/ /seqtk,/'` >> versions.txt
    """
}

process getParams {
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:shac591518dd32ecc3936666c95ff08f6d7474e9728"
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process report {
    label 'process_low'
    container "ontresearch/wf-common:shaabceef445fb63214073cbf5836fdd33c04be4ac7"
    input:
        tuple val(xam_meta), path(vcf)
        file eval_json
        file versions
        path "params.json"
    output:
        path "${xam_meta.id}.report.html", emit: html, optional: true
        path "${xam_meta.id}.svs.json", emit: json
    script:
        def report_name = "${xam_meta.id}.report.html"
        def evalResults = eval_json.name != 'OPTIONAL_FILE' ? "--eval_results ${eval_json}" : ""

    """
    report_sv.py \
        $report_name \
        --vcf $vcf \
        --params params.json \
        --params-hidden 'help,schema_ignore_params,${params.schema_ignore_params}' \
        --versions $versions \
        --revision ${workflow.revision} \
        --commit ${workflow.commitId} \
        --output_json "${xam_meta.id}.svs.json" \
        --workflow_version ${workflow.manifest.version} \
        $evalResults
    """
}
