#!/bin/bash
#
# ResistanceRadar
# Stage 3.3 - Test Decontamination Reference Builder
#
# Author: Rajen K. Goutam, Ph.D.
# Role: Bioinformatics Lead
# Organization: TriAxis Biosciences
#
# Description:
# Builds a small PhiX-only Bowtie2 test decontamination reference
# to validate Stage 3.3 workflow wiring before constructing the full
# production GRCh38 + PhiX174 + pUC19 + pBR322 + UniVec index.
#
# Version: 1.1
# Created: 2026-06-06
#

set -euo pipefail

REF_DIR="/rsstu/users/s/sleblan/MismatchRepair/Project/references/decontamination/test"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

echo "========================================="
echo "ResistanceRadar Stage 3.3"
echo "Build PhiX-only TEST decontamination reference"
echo "========================================="

rm -f decontamination_test.fa decontamination_test*.bt2* phix174.fa phix174.fa.gz pUC19.fa pBR322.fa

echo "Downloading PhiX174..."
curl -fL \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz \
-o phix174.fa.gz

gunzip -f phix174.fa.gz

cp phix174.fa decontamination_test.fa

echo
echo "Reference created:"
ls -lh decontamination_test.fa

echo
echo "Building Bowtie2 index..."

bowtie2-build \
decontamination_test.fa \
decontamination_test

echo
echo "Index files:"
ls -lh decontamination_test*.bt2*
