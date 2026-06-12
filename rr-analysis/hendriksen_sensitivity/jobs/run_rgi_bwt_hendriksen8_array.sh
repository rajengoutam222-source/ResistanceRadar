#!/bin/bash
#BSUB -J HENDRIKSEN8_RGI_BWT[1-8]
#BSUB -q single_chassis
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=64GB]"
#BSUB -W 36:00
#BSUB -o rr-analysis/hendriksen_sensitivity/rgi_bwt/logs/HENDRIKSEN8_RGI_BWT_%I_%J.out
#BSUB -e rr-analysis/hendriksen_sensitivity/rgi_bwt/logs/HENDRIKSEN8_RGI_BWT_%I_%J.err

set -euo pipefail

REPO="/rsstu/users/s/sleblan/MismatchRepair/Project/pipeline/rr-pipeline"
SAMPLESHEET="${REPO}/rr-analysis/hendriksen_sensitivity/metadata/hendriksen_8_bowtie2_clean_samplesheet.csv"

RGI="/rsstu/users/s/sleblan/MismatchRepair/Project/envs/rr_rgi/bin/rgi"
CARD_DIR="/rsstu/users/s/sleblan/MismatchRepair/Project/references/card"

OUTDIR="${REPO}/rr-analysis/hendriksen_sensitivity/rgi_bwt"
mkdir -p "${OUTDIR}"

LINE=$(awk -F',' -v n="${LSB_JOBINDEX}" 'NR==n+1 {print}' "${SAMPLESHEET}")

if [[ -z "${LINE}" ]]; then
    echo "No samplesheet line found for LSB_JOBINDEX=${LSB_JOBINDEX}"
    exit 1
fi

IFS=',' read -r SAMPLE_ID READ1 READ2 <<< "${LINE}"

echo "========================================="
echo "Hendriksen Stage 3.4 CARD/RGI BWT"
echo "Date: $(date)"
echo "Array index: ${LSB_JOBINDEX}"
echo "Sample ID: ${SAMPLE_ID}"
echo "Read1: ${READ1}"
echo "Read2: ${READ2}"
echo "Output prefix: ${OUTDIR}/${SAMPLE_ID}_rgi_bwt"
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
echo "Checking CARD local database..."
cd "${CARD_DIR}"
"${RGI}" database -v --local

echo
echo "Running RGI BWT..."
"${RGI}" bwt \
  -1 "${READ1}" \
  -2 "${READ2}" \
  -a kma \
  -n 8 \
  -o "${OUTDIR}/${SAMPLE_ID}_rgi_bwt" \
  --local \
  --clean

echo
echo "Checking expected RGI outputs..."
ls -lh "${OUTDIR}/${SAMPLE_ID}_rgi_bwt.overall_mapping_stats.txt"
ls -lh "${OUTDIR}/${SAMPLE_ID}_rgi_bwt.gene_mapping_data.txt"

echo
echo "Completed RGI BWT sample: ${SAMPLE_ID}"
echo "Date: $(date)"
