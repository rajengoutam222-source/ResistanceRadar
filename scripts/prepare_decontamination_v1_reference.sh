#!/bin/bash
#
# ResistanceRadar
# Stage 3.3 - Production Decontamination Reference Preparation
#
# Builds combined FASTA for:
#   - GRCh38.p14 primary assembly
#   - PhiX174
#   - pUC19
#   - pBR322
#   - UniVec
#
# Output:
#   /rsstu/users/s/sleblan/MismatchRepair/Project/references/decontamination/v1.0/decontamination_v1.0.fa
#

set -euo pipefail

REF_DIR="/rsstu/users/s/sleblan/MismatchRepair/Project/references/decontamination/v1.0"
FASTA="${REF_DIR}/decontamination_v1.0.fa"
MANIFEST="${REF_DIR}/decontamination_v1.0_manifest.tsv"
SHA_FILE="${REF_DIR}/decontamination_v1.0.sha256"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

echo "========================================="
echo "ResistanceRadar Stage 3.3"
echo "Prepare production decontamination reference v1.0"
echo "========================================="
echo "Reference directory: ${REF_DIR}"
echo

rm -f "${FASTA}" "${MANIFEST}" "${SHA_FILE}"

echo -e "reference\taccession_or_source\tfile\tstatus" > "${MANIFEST}"

echo "Downloading GRCh38.p14 primary assembly..."
curl -fL \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
  -o GRCh38.p14_genomic.fna.gz

gunzip -f GRCh38.p14_genomic.fna.gz
cat GRCh38.p14_genomic.fna >> "${FASTA}"
echo -e "GRCh38.p14\tGCF_000001405.40\tGRCh38.p14_genomic.fna\tincluded" >> "${MANIFEST}"

echo "Downloading PhiX174..."
curl -fL \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz \
  -o PhiX174_genomic.fna.gz

gunzip -f PhiX174_genomic.fna.gz
cat PhiX174_genomic.fna >> "${FASTA}"
echo -e "PhiX174\tGCF_000819615.1\tPhiX174_genomic.fna\tincluded" >> "${MANIFEST}"

echo "Downloading pUC19..."
curl -fL \
  "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?id=M77789.2&db=nuccore&report=fasta&retmode=text" \
  -o pUC19.fa

cat pUC19.fa >> "${FASTA}"
echo -e "pUC19\tM77789.2\tpUC19.fa\tincluded" >> "${MANIFEST}"

echo "Downloading pBR322..."
curl -fL \
  "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?id=J01749.1&db=nuccore&report=fasta&retmode=text" \
  -o pBR322.fa

cat pBR322.fa >> "${FASTA}"
echo -e "pBR322\tJ01749.1\tpBR322.fa\tincluded" >> "${MANIFEST}"

echo "Downloading UniVec..."
curl -fL \
  https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec \
  -o UniVec.fa

cat UniVec.fa >> "${FASTA}"
echo -e "UniVec\tNCBI_UniVec\tUniVec.fa\tincluded" >> "${MANIFEST}"

echo
echo "Creating SHA256 checksum..."
sha256sum "${FASTA}" > "${SHA_FILE}"

echo
echo "Combined FASTA:"
ls -lh "${FASTA}"

echo
echo "Manifest:"
cat "${MANIFEST}"

echo
echo "SHA256:"
cat "${SHA_FILE}"

echo
echo "Reference preparation complete."
