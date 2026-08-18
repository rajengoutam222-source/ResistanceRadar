# Virginia Annual Wastewater 2017-2025

Standardized ResistanceRadar outputs for nine annual Virginia wastewater
metagenomic samples spanning 2017-2025.

Each sample contains six analysis-stage directories and eight compact files:

- fastp: QC JSON
- bowtie2: decontamination log
- CARD_RGI_BWT: overall and gene-level ARG mapping summaries
- ARG_mobility_overlap: strict/perfect ARG-mobility overlap
- geNomAD: plasmid and virus summaries
- Assembly: MEGAHIT assembly metrics

Large runtime files such as FASTQ, BAM, contig FASTA, full RGI Main outputs,
and full ARG-mobility overlap tables are intentionally excluded.

Samples:
2017_SRR13208898
2018_SRR13208892
2019_SRR28226461
2020_SRR33084338
2021_SRR28226424
2022_SRR28226418
2023_SRR37365438
2024_ERR17001293
2025_ERR17359858
