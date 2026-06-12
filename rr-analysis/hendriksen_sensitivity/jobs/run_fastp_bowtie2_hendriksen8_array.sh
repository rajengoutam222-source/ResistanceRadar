#!/bin/bash
#BSUB -J HENDRIKSEN8_FP_BT[1-8]
#BSUB -q single_chassis
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=32GB]"
#BSUB -W 12:00
#BSUB -o rr-analysis/hendriksen_sensitivity/fastp_bowtie2_8sample/logs/HENDRIKSEN8_FP_BT_%I_%J.out
#BSUB -e rr-analysis/hendriksen_sensitivity/fastp_bowtie2_8sample/logs/HENDRIKSEN8_FP_BT_%I_%J.err

set -euo pipefail

REPO="/rsstu/users/s/sleblan/MismatchRepair/Project/pipeline/rr-pipeline"
SAMPLESHEET="${REPO}/tests/samplesheets/erp015409_hendriksen_8_samples.csv"

OUTDIR="${REPO}/rr-analysis/hendriksen_sensitivity/fastp_bowtie2_8sample"
QC_DIR="${OUTDIR}/qc_reports"
BT_DIR="${OUTDIR}/bowtie2"
CLEAN_DIR="${OUTDIR}/clean_fastq"

FASTP="/rsstu/users/s/sleblan/MismatchRepair/Project/envs/rr_meta/bin/fastp"
BOWTIE2="/rsstu/users/s/sleblan/MismatchRepair/Project/envs/rr_meta/bin/bowtie2"

INDEX="/rsstu/users/s/sleblan/MismatchRepair/Project/references/decontamination/v1.0/decontamination_v1.0"

mkdir -p "${QC_DIR}" "${BT_DIR}" "${CLEAN_DIR}"

cd "${REPO}"

LINE=$(awk -F',' -v n="${LSB_JOBINDEX}" 'NR==n+1 {print}' "${SAMPLESHEET}")

if [[ -z "${LINE}" ]]; then
    echo "No samplesheet line found for LSB_JOBINDEX=${LSB_JOBINDEX}"
    exit 1
fi

IFS=',' read -r SAMPLE_ID READ1 READ2 <<< "${LINE}"

echo "========================================="
echo "Hendriksen FASTP + Bowtie2 array task"
echo "Date: $(date)"
echo "Array index: ${LSB_JOBINDEX}"
echo "Sample ID: ${SAMPLE_ID}"
echo "Read1: ${READ1}"
echo "Read2: ${READ2}"
echo "========================================="

if [[ ! -s "${READ1}" ]]; then
    echo "Missing READ1: ${READ1}"
    exit 1
fi

if [[ ! -s "${READ2}" ]]; then
    echo "Missing READ2: ${READ2}"
    exit 1
fi

echo
echo "Running FASTP..."

"${FASTP}" \
  -i "${READ1}" \
  -I "${READ2}" \
  -o "${CLEAN_DIR}/${SAMPLE_ID}_fastp_clean_R1.fastq.gz" \
  -O "${CLEAN_DIR}/${SAMPLE_ID}_fastp_clean_R2.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --length_required 50 \
  --correction \
  --overlap_len_require 10 \
  --thread 4 \
  --json "${QC_DIR}/${SAMPLE_ID}_fastp.json" \
  --html "${QC_DIR}/${SAMPLE_ID}_fastp.html"

echo
echo "Running Bowtie2 host/contaminant removal..."

"${BOWTIE2}" \
  --very-sensitive-local \
  -x "${INDEX}" \
  -1 "${CLEAN_DIR}/${SAMPLE_ID}_fastp_clean_R1.fastq.gz" \
  -2 "${CLEAN_DIR}/${SAMPLE_ID}_fastp_clean_R2.fastq.gz" \
  --threads 8 \
  --un-conc-gz "${CLEAN_DIR}/${SAMPLE_ID}_clean.fastq.gz" \
  > /dev/null 2> "${BT_DIR}/${SAMPLE_ID}_bowtie2.log"

mv "${CLEAN_DIR}/${SAMPLE_ID}_clean.fastq.1.gz" "${CLEAN_DIR}/${SAMPLE_ID}_bowtie2_clean_R1.fastq.gz"
mv "${CLEAN_DIR}/${SAMPLE_ID}_clean.fastq.2.gz" "${CLEAN_DIR}/${SAMPLE_ID}_bowtie2_clean_R2.fastq.gz"

echo
echo "Creating host-removal summary TSV..."

INPUT_READ_PAIRS=$(awk 'NR==1 {print $1}' "${BT_DIR}/${SAMPLE_ID}_bowtie2.log")
CONCORDANTLY_UNMAPPED_PAIRS=$(grep "aligned concordantly 0 times" "${BT_DIR}/${SAMPLE_ID}_bowtie2.log" | head -1 | awk '{print $1}')
CONCORDANTLY_MAPPED_ONCE=$(grep "aligned concordantly exactly 1 time" "${BT_DIR}/${SAMPLE_ID}_bowtie2.log" | awk '{print $1}')
CONCORDANTLY_MAPPED_MULTI=$(grep "aligned concordantly >1 times" "${BT_DIR}/${SAMPLE_ID}_bowtie2.log" | awk '{print $1}')
OVERALL_ALIGNMENT_RATE=$(grep "overall alignment rate" "${BT_DIR}/${SAMPLE_ID}_bowtie2.log" | awk '{print $1}' | sed 's/%//')
TOOL_VERSION=$("${BOWTIE2}" --version | head -1 | awk '{print $3}')

{
    printf "sample_id\tstage\ttool\ttool_version\tbowtie2_mode\tbowtie2_preset\tdecontamination_reference\tbowtie2_index\tinput_read_pairs\tconcordantly_unmapped_pairs\tconcordantly_mapped_once\tconcordantly_mapped_multi\toverall_alignment_rate_percent\tpercent_removed\tfinal_clean_read_pairs\tclean_R1\tclean_R2\tstatus\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${SAMPLE_ID}" \
        "host_removal" \
        "bowtie2" \
        "${TOOL_VERSION}" \
        "paired_end" \
        "very_sensitive_local" \
        "GRCh38p14_PhiX174_pUC19_pBR322_UniVec_v1.0" \
        "${INDEX}" \
        "${INPUT_READ_PAIRS}" \
        "${CONCORDANTLY_UNMAPPED_PAIRS}" \
        "${CONCORDANTLY_MAPPED_ONCE}" \
        "${CONCORDANTLY_MAPPED_MULTI}" \
        "${OVERALL_ALIGNMENT_RATE}" \
        "${OVERALL_ALIGNMENT_RATE}" \
        "${CONCORDANTLY_UNMAPPED_PAIRS}" \
        "${CLEAN_DIR}/${SAMPLE_ID}_bowtie2_clean_R1.fastq.gz" \
        "${CLEAN_DIR}/${SAMPLE_ID}_bowtie2_clean_R2.fastq.gz" \
        "completed"
} > "${BT_DIR}/${SAMPLE_ID}_host_removal_log.tsv"

echo
echo "Completed sample: ${SAMPLE_ID}"
echo "Date: $(date)"
