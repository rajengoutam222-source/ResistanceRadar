/*
 * ResistanceRadar
 * Stage 3.3 - Host and Contamination Removal
 *
 * Author: Rajen K. Goutam, Ph.D.
 * Role: Bioinformatics Lead
 * Organization: TriAxis Biosciences
 *
 * Bowtie2 2.5.3
 */

process BOWTIE2_HOST_REMOVAL {

    tag "${sample_id}"

    cpus 8

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}_bowtie2_clean_R1.fastq.gz"),
          path("${sample_id}_bowtie2_clean_R2.fastq.gz")
    path("${sample_id}_host_removal_log.tsv")

    script:
    """
INDEX=/rsstu/users/s/sleblan/MismatchRepair/Project/references/decontamination/v1.0/decontamination_v1.0

       STAGE="host_removal"
    TOOL="bowtie2"
    TOOL_VERSION=\$(bowtie2 --version | head -1 | awk '{print \$3}')
    BOWTIE2_MODE="paired_end"
    BOWTIE2_PRESET="very_sensitive_local"
    DECONTAMINATION_REFERENCE="GRCh38p14_PhiX174_pUC19_pBR322_UniVec_v1.0"
 bowtie2 \
        --very-sensitive-local \
        -x \$INDEX \
        -1 ${read1} \
        -2 ${read2} \
        --threads ${task.cpus} \
        --un-conc-gz ${sample_id}_clean.fastq.gz \
        2> ${sample_id}_bowtie2.log
    mv ${sample_id}_clean.fastq.1.gz ${sample_id}_bowtie2_clean_R1.fastq.gz
    mv ${sample_id}_clean.fastq.2.gz ${sample_id}_bowtie2_clean_R2.fastq.gz   


    INPUT_READ_PAIRS=\$(awk 'NR==1 {print \$1}' ${sample_id}_bowtie2.log)

    CONCORDANTLY_UNMAPPED_PAIRS=\$(grep "aligned concordantly 0 times" ${sample_id}_bowtie2.log | head -1 | awk '{print \$1}')

    CONCORDANTLY_MAPPED_ONCE=\$(grep "aligned concordantly exactly 1 time" ${sample_id}_bowtie2.log | awk '{print \$1}')

    CONCORDANTLY_MAPPED_MULTI=\$(grep "aligned concordantly >1 times" ${sample_id}_bowtie2.log | awk '{print \$1}')

    OVERALL_ALIGNMENT_RATE=\$(grep "overall alignment rate" ${sample_id}_bowtie2.log | awk '{print \$1}' | sed 's/%//')

    FINAL_CLEAN_READ_PAIRS=\$CONCORDANTLY_UNMAPPED_PAIRS

        printf "sample_id\tstage\ttool\ttool_version\tbowtie2_mode\tbowtie2_preset\tdecontamination_reference\tbowtie2_index\tinput_read_pairs\tconcordantly_unmapped_pairs\tconcordantly_mapped_once\tconcordantly_mapped_multi\toverall_alignment_rate_percent\tpercent_removed\tfinal_clean_read_pairs\tclean_R1\tclean_R2\tstatus\n" > ${sample_id}_host_removal_log.tsv

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${sample_id}" \
        "\$STAGE" \
        "\$TOOL" \
        "\$TOOL_VERSION" \
        "\$BOWTIE2_MODE" \
        "\$BOWTIE2_PRESET" \
        "\$DECONTAMINATION_REFERENCE" \
        "\$INDEX" \
        "\$INPUT_READ_PAIRS" \
        "\$CONCORDANTLY_UNMAPPED_PAIRS" \
        "\$CONCORDANTLY_MAPPED_ONCE" \
        "\$CONCORDANTLY_MAPPED_MULTI" \
        "\$OVERALL_ALIGNMENT_RATE" \
        "\$OVERALL_ALIGNMENT_RATE" \
        "\$FINAL_CLEAN_READ_PAIRS" \
        "${sample_id}_bowtie2_clean_R1.fastq.gz" \
        "${sample_id}_bowtie2_clean_R2.fastq.gz" \
        "completed" >> ${sample_id}_host_removal_log.tsv





    """
}
