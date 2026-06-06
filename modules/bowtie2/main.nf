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
    INDEX=/rsstu/users/s/sleblan/MismatchRepair/Project/references/decontamination/test/decontamination_test
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

    echo -e "sample_id\\tstatus" > ${sample_id}_host_removal_log.tsv
    echo -e "${sample_id}\\tcompleted" >> ${sample_id}_host_removal_log.tsv
    """
}
