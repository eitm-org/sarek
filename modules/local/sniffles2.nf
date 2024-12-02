import groovy.json.JsonBuilder

process SNIFFLES2 {
    tag "$meta.id"
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:shac591518dd32ecc3936666c95ff08f6d7474e9728"

    input:
        tuple val(meta), path(xam), path(xam_idx)
        file tr_bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH) 
        val genome_build
    output:
        tuple val(meta), path("*.sniffles.vcf"), emit: vcf
        path "${xam}.wf_sv.snf", emit: snf
        path "versions.yml", emit: versions
    script:
        // if tr_arg is not provided and genome_build is set
        // automatically pick the relevant TR BED from the SV image
        def tr_arg = ""
        if (tr_bed.name != 'OPTIONAL_FILE'){
            tr_arg = "--tandem-repeats ${tr_bed}"
        }
        else if (genome_build) {
            log.warn "Automatically selecting TR BED: ${genome_build}.trf.bed"
            tr_arg = "--tandem-repeats \${WFSV_TRBED_PATH}/${genome_build}.trf.bed"
        }
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
        sniffles2: \$(echo \$(sniffles --version 2>&1) | sed 's/^.*(sniffles) v//; s/ .*\$//')
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
        val chromosome_codes
    output:
        tuple val(meta), path("*.filtered.vcf"), emit: vcf
    script:
    String ctgs = chromosome_codes.join(',')
    def ctgs_filter = params.include_all_ctgs ? "" : "--contigs ${ctgs}"
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

process getGenome {
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:shac591518dd32ecc3936666c95ff08f6d7474e9728"

    input:
        tuple val(meta), path(xam), path(xam_idx)
    output:
        env genome_build, emit: genome_build, optional: false
     script:
        // set flags for subworkflows that have genome build restrictions
        def str_arg = ""
        def cnv_arg = ""
        def qdnaseq_arg = ""
        """
        # use view -H rather than idxstats, as idxstats will still cause a scan of the whole CRAM (https://github.com/samtools/samtools/issues/303)
        samtools view -H ${xam} --no-PG | grep '^@SQ' | sed -nE 's,.*SN:([^[:space:]]*).*LN:([^[:space:]]*).*,\\1\\t\\2,p' > ${xam}_genome.txt
        get_genome.py --chr_counts ${xam}_genome.txt -o output.txt ${str_arg} ${cnv_arg} ${qdnaseq_arg}
        genome_build=`cat output.txt`
        """
}