
include {
    SNIFFLES2;
    filterCalls;
    sortVCF;
} from "../../../modules/local/sniffles2.nf"
include {
    filterBenchmarkVcf;
    intersectBedWithTruthset;
    truvari;
} from "../../../modules/local/sniffles2_eval.nf"
include { 
    MOSDEPTH           
} from '../../../modules/nf-core/mosdepth/main'

workflow BAM_VARIANT_CALLING_STRUCTURAL_SNIFFLES2 {
    take:
        bam_channel                                 // channel: [mandatory] [meta, bam, bai]
        reference                                   // channel: [mandatory] [fa]
        target                                      // channel: [mandatory] [bed]
        mosdepth_stats
        optional_file

    main:

        ch_versions = Channel.empty()

        sniffles2 = SNIFFLES2(bam_channel.map{ meta, xam, xai -> [meta, xam, xai] }, optional_file, reference)
        filterCalls(sniffles2.vcf, mosdepth_stats, target)
        sorted_vcf = sortVCF(filterCalls.out.vcf)

        hg002_seq = false
        if (hg002_seq) {
            benchmark_result = runBenchmark(sorted_vcf.vcf, reference, target)
        } else {
            benchmark_result = Channel.fromPath(optional_file)
        }

        ch_versions = ch_versions.mix(sniffles2.versions)

    emit:
        sniffles2_vcf = sorted_vcf.vcf_gz
        versions = ch_versions // software versions
}


workflow runBenchmark {
    take:
        vcf
        reference
        target
    main:
        // for benchmarking we bundle a dataset in the SV container in $WFSV_EVAL_DATA_PATH
        // rather than coupling that dataset to the workflow by referring to it here
        //   in a value channel (or similar), we'll instead interpret use of dummy files
        //   as a flag to load from the bundled dataset inside the process scope
        // note we're not using the usual `optional_file` as this will cause an input collision error
        //   instead we just reference some OPTIONAL_FILE.ext that we know don't exist
        //   we can get away with this as the files will never be opened (so don't need to exist)

        // reconcile workflow target BED and benchmark truthset BED
        //   recall if user does not input a BED, one covering all genomic
        //   intervals in the ref is generated by getAllChromosomesBed
        params.sv_benchmark_bed = false
        if (params.sv_benchmark_bed) {
            truthset_bed = Channel.fromPath(params.sv_benchmark_bed, checkIfExists: true)
        }
        else {
            truthset_bed = file("OPTIONAL_FILE.bed") // this will trigger process to use bundled benchmark bed
        }
        
        intersected = intersectBedWithTruthset(target, truthset_bed)

        // load user-provided benchmark data
        if (params.sv_benchmark_vcf) {
            // truvari assumes index is [vcf].tbi
            truthset_vcf = Channel.fromPath(params.sv_benchmark_vcf, checkIfExists: true)
            truthset_tbi = Channel.fromPath(params.sv_benchmark_vcf + '.tbi', checkIfExists: true)
        }
        else {
            // we'll create some non-existent optional files to stage
            // again this will trigger the process to use the bundled benchmark data
            // we use channels here so we can concat them later
            truthset_vcf = Channel.fromPath("OPTIONAL_FILE.vcf.gz", checkIfExists: false)
            truthset_tbi = Channel.fromPath("OPTIONAL_FILE.vcf.gz.tbi", checkIfExists: false)
        }

        // run benchmark
        filtered = filterBenchmarkVcf(vcf)
        truvari(
            reference,
            filtered,
            truthset_vcf.concat(truthset_tbi).toList(),
            intersected.intersected_bed)
    emit:
        json = truvari.out.truvari_json
}