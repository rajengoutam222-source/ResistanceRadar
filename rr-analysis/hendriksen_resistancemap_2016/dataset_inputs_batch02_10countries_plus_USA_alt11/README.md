# ResistanceRadar RM2016 Batch02 compact dataset inputs

This directory is the standardized, GitHub-safe release for the
RM2016 Batch02 cohort: 10 countries plus one alternate USA sample.

## Cohort

- Samples: 11
- Sample IDs: ERR1713332, ERR1713334, ERR1713337, ERR1713335, ERR1713339, ERR1713349, ERR1713352, ERR1713362, ERR1713379, ERR1713343, ERR1713403
- FASTP, Bowtie2, CARD/RGI BWT, DIAMOND/CARD, MEGAHIT,
  geNomAD, and CARD/RGI Main stages completed.

## Directory organization

- `FASTP/<sample>/`: FASTP JSON quality-control report.
- `Bowtie2/<sample>/`: Bowtie2 decontamination log.
- `CARD_RGI_BWT/<sample>/`: read-level RGI BWT mapping statistics
  and gene-mapping data.
- `DIAMOND_CARD/<sample>/`: compact CARD hit and sample summaries.
- `Assembly/<sample>/`: assembly statistics calculated from the
  final MEGAHIT contigs.
- `geNomAD/<sample>/contig_classification.tsv`: standardized copy
  of the geNomAD aggregated contig classification.
- `CARD_RGI_Main/<sample>/`: compact contig-level RGI Main table.
- `metadata/`: sample manifest.

## Files intentionally excluded

The release does not include raw FASTQ reads, FASTP-cleaned FASTQ
reads, Bowtie2-cleaned FASTQ reads, full MEGAHIT contig FASTA files,
MEGAHIT temporary working files, raw per-read DIAMOND tables, geNomAD
sequence/protein files, or RGI JSON files.

The compact RGI Main tables omit the sequence-heavy columns
`Predicted_DNA`, `Predicted_Protein`, and `CARD_Protein_Sequence`.

## Integrity and provenance

- `PROVENANCE.tsv` maps every released analytical file to its source
  and transformation.
- `FILE_INDEX.tsv` lists released files and sizes.
- `SHA256SUMS.txt` provides SHA-256 checksums.
