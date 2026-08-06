Kathmandu KBase 18-Sample Standardized Analysis Release
==========================================================

Generated: 2026-08-06T17:48:03Z
Unique samples: 18
Portable stage files: 270
Symbolic links: 0

Important FASTP note
--------------------
FASTP JSON reports are available for all 18 samples. The I_5-6 report
was regenerated on 2026-08-06 by streaming the original KBase
interleaved FASTQ ZIP through fastp v1.3.3 in report-only mode.
Filtered FASTQ output was discarded, and the existing Bowtie2 and
downstream results were not rerun or changed.
For E_1-2, the interleaved FASTP result was selected instead of the older
paired-file result to prevent duplicate delivery.

Release contents
----------------
FASTP
  FASTP JSON reports for all 18 samples.

Bowtie2
  Bowtie2 decontamination summary files for all 18 samples.

CARD_RGI_BWT
  RGI BWT gene-mapping and overall-mapping statistics for all 18 samples.

DIAMOND_CARD
  Streaming compact DIAMOND/CARD hit and sample summaries for all 18
  samples. The original 500–600 MB BLASTX tables are not included.
  Summary statistics include hit counts, ARO/CARD identities, mean
  percent identity, mean bitscore, and maximum bitscore.

  Unlike the Hendriksen cohort releases, Kathmandu does not include
  six-column per-read DIAMOND tables. A tested gzip-compressed version
  exceeded GitHub's 100 MiB per-file limit for at least one sample
  (E_5-6: 149.44 MiB). Derived summaries are therefore used as the
  portable Git release format, while provenance records identify the
  original complete DIAMOND source tables.

Assembly
  Metrics calculated from resolved final MEGAHIT contigs for all 18
  samples. Full FASTA assemblies are not included.

geNomAD
  Plasmid and virus summary tables for all 18 samples.

CARD_RGI_Main
  Compact RGI Main tables for all 18 samples. All rows are retained.
  Predicted_DNA, Predicted_Protein, and CARD_Protein_Sequence were removed.

ARG_Mobility_Overlap
  Five per-sample ARG-mobility tables for all 18 samples.

Metadata
--------
sample_manifest.tsv
source_selection.tsv
dataset_inputs_index.tsv
file_checksums.sha256
rgi_main_compact_provenance.tsv
diamond_card_summary_provenance.tsv

Stage file counts
-----------------
FASTP: 18
Bowtie2: 18
CARD_RGI_BWT: 36
DIAMOND_CARD: 36
Assembly: 18
geNomAD: 36
CARD_RGI_Main: 18
ARG_Mobility_Overlap: 90
