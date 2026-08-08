# ResistanceRadar RM2016 Batch03 global10 compact dataset inputs

This directory is the standardized, GitHub-safe release for the
RM2016 Batch03 global10 cohort.

## Cohort

- Samples: 10
- Sample IDs: ERR1713346, ERR1713370, ERR1713382, ERR1713385, ERR1713386, ERR1725938, ERR1725952, ERR1725965, ERR1725980, ERR1725981
- FASTP, Bowtie2, CARD/RGI BWT, DIAMOND/CARD, MEGAHIT,
  geNomAD, CARD/RGI Main, and ARG-mobility overlap stages completed.

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
- `ARG_Mobility_Overlap/`: ARG-to-mobile-element overlap, strict/perfect subset, sample summaries, and drug-class summary.
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
