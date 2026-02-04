// filter VCF such that only INS, DEL and DUP remain for benchmarking
// possibly make this optional in future but there is only one canonical data set for this benchmarking currently
process filterBenchmarkVcf {
    //tag "$meta.id"
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:sha8134f9fef5e19605c7fb4c1348961d6771f1af79"

    input:
        tuple val(xam_meta), path(calls_vcf)
    output:
        tuple val(xam_meta), path("benchmarkCalls.vcf.gz"), path("benchmarkCalls.vcf.gz.tbi")
    script:
    """
    zcat $calls_vcf \
    | bcftools view -i '(SVTYPE = \"INS\" || SVTYPE = \"DEL\" || SVTYPE = \"DUP\")' \
    | bgziptabix benchmarkCalls.vcf.gz
    """
}

// ensure to output the exact intersection of A and B to only benchmark
// against regions that have been analysed by the workflow
process intersectBedWithTruthset {
    //tag "$meta.id"
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:sha8134f9fef5e19605c7fb4c1348961d6771f1af79"
    
    input:
        path target_bed
        path user_truthset_bed
    output:
        path "target_truthset.bed", emit: intersected_bed
    script:
    // use the bundled benchmark BED if the user_truthset_bed is the dummy
    def tru_bed_arg = user_truthset_bed.name.startsWith("OPTIONAL_FILE") ? "\${WFSV_EVAL_DATA_PATH}/benchmark.bed" : user_truthset_bed

    """
    bedtools intersect \
        -a ${tru_bed_arg} \
        -b $target_bed \
        > target_truthset.bed
    if [ ! -s target_truthset.bed ]
    then
        echo "No overlaps found between truth and target"
        echo "Chr names in your target or reference and truthset may differ"
        exit 1
    fi
    """
}

// --dup-to-ins included to handle the hg2 NIST data set where DUP are coded as INS
process truvari {
    //tag "$meta.id"
    label 'process_low'
    container "ontresearch/wf-human-variation-sv:sha8134f9fef5e19605c7fb4c1348961d6771f1af79"
    
    input:
        path(ref) //, path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple val(xam_meta), path(calls_vcf), path(calls_vcf_tbi)
        tuple path(user_truthset_vcf), path(user_truthset_tbi)
        file include_bed
    output:
        path "*.truvari.json", emit: truvari_json
    script:
    // use the bundled benchmark data if the user_truthset_vcf is the dummy
    def tru_vcf_arg = user_truthset_vcf.name.startsWith("OPTIONAL_FILE") ? "\${WFSV_EVAL_DATA_PATH}/benchmark.vcf.gz" : user_truthset_vcf
    """
    
    bcftools norm -f $ref -m -any $calls_vcf -Oz -o ${calls_vcf}.split.vcf
    bgzip ${calls_vcf}.split.vcf
    bcftools index -t ${calls_vcf}.split.vcf.gz

    truvari benc \
        --passonly \
        --dup-to-ins \
        --pick ac \
        -b ${tru_vcf_arg} \
        -c ${calls_vcf}.split.vcf.gz \
        -f ${ref} \
        -o ${xam_meta.id} \
        --includebed $include_bed
/*
    truvari refine \
        --recount \
        --use-region-coords \
        --use-original-vcfs \
        --align mafft \
        ${xam_meta.id}
*/ 
    mv ${xam_meta.id}/summary.txt ${xam_meta.id}.truvari.json
    """
}