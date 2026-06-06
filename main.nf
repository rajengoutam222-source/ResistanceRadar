nextflow.enable.dsl=2

include { FASTP } from './modules/fastp/main.nf'
include { BOWTIE2_HOST_REMOVAL } from './modules/bowtie2/main.nf'

params.samplesheet = "tests/samplesheets/erp015409_samples.csv"

workflow {
    reads_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample_id, file(row.read1), file(row.read2))
        }

    FASTP(reads_ch)

    bowtie2_input_ch = FASTP.out.map { sample_id, clean_r1, clean_r2, fastp_json, fastp_html ->
        tuple(sample_id, clean_r1, clean_r2)
    }

    BOWTIE2_HOST_REMOVAL(bowtie2_input_ch)
}
