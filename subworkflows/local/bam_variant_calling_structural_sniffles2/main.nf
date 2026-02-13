
include {
    SNIFFLES2;
    filterCalls;
    sortVCF;
    getVersions;
    getParams;
    report;
} from "../../../modules/local/sniffles2.nf"
include { 
    MOSDEPTH           
} from '../../../modules/nf-core/mosdepth/main'

workflow BAM_VARIANT_CALLING_STRUCTURAL_SNIFFLES2 {
    take:
        xam_channel                                 // channel: [mandatory] [meta, xam, xai]
        reference                                   // channel: [mandatory] [fa]
        target                                      // channel: [mandatory] [bed]
        filter                                      // channel: [mandatory] [bed]
        mosdepth_stats
        optional_file
        sv_benchmark_bed
        sv_benchmark_vcf
        hg002_seq                                   // set as true if you have hg002 data from sequencer

    main:

        ch_versions = Channel.empty()

        sniffles2 = SNIFFLES2(xam_channel.map{ meta, xam, xai -> [meta, xam, xai] }, target, reference)  // run Sniffles2!
        filterCalls(sniffles2.vcf, mosdepth_stats, target, filter) // removes the pathogenic regions + applies read support filters
        sorted_vcf = sortVCF(filterCalls.out.vcf) // sorts VCF file by chromosome and position and then indexes

        final_vcf = sorted_vcf.vcf_gz.join(sorted_vcf.vcf_tbi) // join output vcf with respective index file

        // if (hg002_seq) {
        //    benchmark_result = runBenchmark(sorted_vcf.vcf_gz, reference, target, sv_benchmark_bed, sv_benchmark_vcf).first()
        // } else {
        //    benchmark_result = Channel.fromPath(optional_file).first()
        // }

        // run report.py to get nice, clean stats/counts from VCFs
        // report = runReport(
        //    sorted_vcf.vcf_gz.groupTuple(),
        //    benchmark_result
        //)
        // sv_stats_json = report.json
        // report = report.html.concat(
        //    final_vcf.map{meta, vcf, tbi -> [vcf, tbi]}
        // )

        // remap to add variantcaller info
        sniffles2_vcf = Channel.empty().mix(sorted_vcf.vcf_gz).map{ meta, vcf ->
                [[
                    id:             meta.sample,
                    num_intervals:  meta.num_intervals,
                    patient:        meta.patient,
                    sample:         meta.sample,
                    sex:            meta.sex,
                    status:         meta.status,
                    variantcaller:  "sniffles2"
                ],vcf]
            }

        ch_versions = ch_versions.mix(sniffles2.versions)
        
    emit:
        // report = report
        // sv_stats_json = sv_stats_json
        sniffles2_vcf = sniffles2_vcf
        versions = ch_versions // software versions
}
