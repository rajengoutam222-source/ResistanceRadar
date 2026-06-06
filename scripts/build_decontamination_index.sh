#!/bin/bash
#
# ResistanceRadar
# Stage 3.3 - Host and Contamination Removal
#
# Author: Rajen K. Goutam, Ph.D.
# Role: Bioinformatics Lead
# Organization: TriAxis Biosciences
#
# Description:
# Builds the Bowtie2 decontamination index used for host and
# contaminant read removal in the ResistanceRadar pipeline.
#
# References included by specification:
#   - GRCh38 primary assembly
#   - PhiX174
#   - pUC19
#   - pBR322
#   - UniVec
#
# Version: 1.0
# Created: 2026-06-06
#

set -euo pipefail

REF_DIR="/rsstu/users/s/sleblan/MismatchRepair/Project/references/decontamination/v1.0"
FASTA="${REF_DIR}/decontamination_v1.0.fa"
INDEX_PREFIX="${REF_DIR}/decontamination_v1.0"

mkdir -p "${REF_DIR}"

echo "========================================="
echo "ResistanceRadar Stage 3.3"
echo "Build Bowtie2 decontamination index"
echo "Target references:"
echo "  - GRCh38 primary assembly"
echo "  - PhiX174"
echo "  - pUC19"
echo "  - pBR322"
echo "  - UniVec"
echo "========================================="

echo "Reference directory:"
echo "${REF_DIR}"

echo
echo "NOTE:"
echo "This script currently expects the combined FASTA to exist at:"
echo "${FASTA}"

echo
if [[ ! -s "${FASTA}" ]]; then
    echo "ERROR: Combined FASTA not found:"
    echo "${FASTA}"
    echo
    echo "Next step: download/prepare GRCh38, PhiX174, pUC19, pBR322, and UniVec,"
    echo "then concatenate them into decontamination_v1.0.fa."
    exit 1
fi

echo "Building Bowtie2 index..."
bowtie2-build "${FASTA}" "${INDEX_PREFIX}"

echo "Index complete:"
ls -lh "${REF_DIR}"/decontamination_v1.0*.bt2*
