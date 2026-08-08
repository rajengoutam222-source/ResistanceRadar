#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import statistics
import sys
from collections import defaultdict
from pathlib import Path


REPO = Path.cwd()

BASE = (
    REPO
    / "rr-analysis"
    / "hendriksen_resistancemap_2016"
)

TARGET = BASE / "dataset_inputs_USA_crossstudy_additional8"
BUILD = BASE / ".dataset_inputs_USA_crossstudy_additional8.building"

MEGA_STAGE = BASE / "stage3_5_megahit_USA_crossstudy_additional8"

MEGA_SHEET = (
    MEGA_STAGE
    / "samplesheets"
    / "USA_crossstudy_additional8_megahit_samples.tsv"
)

FASTP_ROOT = BASE / "fastp" / "USA_crossstudy_additional8"

BOWTIE2_ROOT = (
    BASE
    / "results"
    / "bowtie2_USA_crossstudy_additional8"
)

RGI_MAIN_ROOT = (
    BASE
    / "stage3_5_rgi_main_contigs_USA_crossstudy_additional8"
    / "results"
)

GENOMAD_ROOT = (
    BASE
    / "stage3_5_genomad_USA_crossstudy_additional8"
    / "results"
)

DIAMOND_ROOT = (
    BASE
    / "stage3_4c_diamond_card_USA_crossstudy_additional8"
    / "results_cleanfasta"
)

RGI_BWT_ROOT_CANDIDATES = [
    BASE / "results" / "rgi_bwt_USA_crossstudy_additional8",
]

DROP_RGI_MAIN_COLUMNS = {
    "Predicted_DNA",
    "Predicted_Protein",
    "CARD_Protein_Sequence",
}

EXPECTED_SAMPLES = 8
GITHUB_LIMIT = 100 * 1024 * 1024

ARG_MOBILITY_ROOT = (
    BASE
    / "stage3_5_arg_mobility_overlap_USA_crossstudy_additional8"
)

ARG_MOBILITY_PREFIX = "rm2016_USA_crossstudy_additional8"


def require_file(path: Path, allow_empty: bool = False) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")

    if not allow_empty and path.stat().st_size == 0:
        raise ValueError(f"Empty file: {path}")

    return path


def repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO.resolve()))
    except ValueError:
        return str(path.resolve())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def destination_path(relative_path: Path) -> Path:
    path = BUILD / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


provenance: list[dict[str, str]] = []


def add_provenance(
    destination: Path,
    sources: list[Path],
    transformation: str,
) -> None:
    provenance.append(
        {
            "destination": str(destination),
            "source": ";".join(repo_relative(path) for path in sources),
            "transformation": transformation,
        }
    )


def copy_release_file(
    source: Path,
    relative_destination: Path,
    transformation: str = "Copied without modification",
    allow_empty: bool = False,
) -> None:
    require_file(source, allow_empty=allow_empty)

    destination = destination_path(relative_destination)
    shutil.copy2(source, destination)

    add_provenance(
        relative_destination,
        [source],
        transformation,
    )


def copy_normalized_rgi_bwt_file(
    source: Path,
    relative_destination: Path,
) -> None:
    """Copy an RGI BWT text file with Git-safe line formatting.

    This normalizes CRLF/CR line endings to LF, removes trailing ASCII
    spaces, and removes blank lines only at the end of the file. Tabs
    are preserved so tabular fields remain unchanged.
    """
    require_file(source)

    destination = destination_path(relative_destination)

    data = source.read_bytes()
    normalized = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    lines = normalized.split(b"\n")

    while lines and lines[-1] == b"":
        lines.pop()

    lines = [line.rstrip(b" ") for line in lines]

    if not lines:
        raise ValueError(f"Empty RGI BWT file after normalization: {source}")

    destination.write_bytes(b"\n".join(lines) + b"\n")

    add_provenance(
        relative_destination,
        [source],
        (
            "Copied with LF line endings; trailing ASCII spaces and "
            "terminal blank lines removed; tabs preserved"
        ),
    )


def find_rgi_bwt_file(sample: str, suffix: str) -> Path:
    filename = f"{sample}_rgi_bwt.{suffix}"

    candidates: list[Path] = []

    for root in RGI_BWT_ROOT_CANDIDATES:
        candidates.extend(
            [
                root / filename,
                root / sample / filename,
            ]
        )

    existing = [path for path in candidates if path.is_file()]

    if len(existing) == 1:
        return existing[0]

    if len(existing) > 1:
        raise RuntimeError(
            f"Multiple RGI BWT sources found for {sample}: "
            + ", ".join(str(path) for path in existing)
        )

    attempted = "\n  ".join(str(path) for path in candidates)

    raise FileNotFoundError(
        f"Could not locate {filename}. Tried:\n  {attempted}"
    )


def compact_rgi_main(source: Path, destination: Path) -> None:
    require_file(source)

    destination.parent.mkdir(parents=True, exist_ok=True)

    with source.open("r", newline="", errors="replace") as src:
        reader = csv.reader(src, delimiter="\t")

        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty RGI Main table: {source}") from exc

        keep_indexes = [
            index
            for index, column in enumerate(header)
            if column not in DROP_RGI_MAIN_COLUMNS
        ]

        with destination.open("w", newline="") as dst:
            writer = csv.writer(
                dst,
                delimiter="\t",
                lineterminator="\n",
            )

            writer.writerow([header[index] for index in keep_indexes])

            for row in reader:
                padded = row + [""] * max(0, len(header) - len(row))
                writer.writerow([padded[index] for index in keep_indexes])


def fasta_metrics(path: Path) -> dict[str, str]:
    require_file(path)

    lengths: list[int] = []
    total_gc = 0
    total_bases = 0
    current_length = 0

    with path.open("r", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_length:
                    lengths.append(current_length)
                current_length = 0
                continue

            sequence = line.upper()
            current_length += len(sequence)
            total_bases += len(sequence)
            total_gc += sequence.count("G") + sequence.count("C")

    if current_length:
        lengths.append(current_length)

    if not lengths:
        raise ValueError(f"No FASTA sequences found in {path}")

    descending = sorted(lengths, reverse=True)
    half_total = sum(lengths) / 2
    cumulative = 0
    n50 = 0
    l50 = 0

    for index, length in enumerate(descending, start=1):
        cumulative += length

        if cumulative >= half_total:
            n50 = length
            l50 = index
            break

    return {
        "contig_count": str(len(lengths)),
        "total_length_bp": str(sum(lengths)),
        "minimum_contig_bp": str(min(lengths)),
        "maximum_contig_bp": str(max(lengths)),
        "mean_contig_bp": f"{statistics.mean(lengths):.2f}",
        "median_contig_bp": f"{statistics.median(lengths):.2f}",
        "N50_bp": str(n50),
        "L50": str(l50),
        "GC_percent": (
            f"{100 * total_gc / total_bases:.4f}"
            if total_bases
            else "0.0000"
        ),
        "source_fasta_bytes": str(path.stat().st_size),
    }


def write_assembly_metrics(
    sample: str,
    source: Path,
    destination: Path,
) -> None:
    metrics = fasta_metrics(source)

    destination.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "sample_id",
        "contig_count",
        "total_length_bp",
        "minimum_contig_bp",
        "maximum_contig_bp",
        "mean_contig_bp",
        "median_contig_bp",
        "N50_bp",
        "L50",
        "GC_percent",
        "source_fasta_bytes",
    ]

    row = {"sample_id": sample, **metrics}

    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(row)


def summarize_diamond(
    sample: str,
    r1_path: Path,
    r2_path: Path,
    hit_destination: Path,
    sample_destination: Path,
) -> None:
    require_file(r1_path, allow_empty=True)
    require_file(r2_path, allow_empty=True)

    aggregated: dict[str, dict[str, float]] = defaultdict(
        lambda: {
            "total": 0,
            "r1": 0,
            "r2": 0,
            "sum_pident": 0.0,
            "max_length": 0,
            "min_evalue": float("inf"),
            "max_bitscore": float("-inf"),
        }
    )

    totals = {
        "R1_hits": 0,
        "R2_hits": 0,
        "all_hits": 0,
        "sum_pident": 0.0,
        "max_bitscore": float("-inf"),
        "min_evalue": float("inf"),
    }

    def process(path: Path, read_label: str) -> None:
        with path.open("r", errors="replace") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\n")

                if not line:
                    continue

                fields = line.split("\t")

                if len(fields) < 6:
                    raise ValueError(
                        f"Malformed DIAMOND row: {path}:{line_number}"
                    )

                _, subject, pident, alignment_length, evalue, bitscore = (
                    fields[:6]
                )

                pident_value = float(pident)
                length_value = int(float(alignment_length))
                evalue_value = float(evalue)
                bitscore_value = float(bitscore)

                record = aggregated[subject]
                record["total"] += 1
                record[read_label.lower()] += 1
                record["sum_pident"] += pident_value
                record["max_length"] = max(
                    record["max_length"],
                    length_value,
                )
                record["min_evalue"] = min(
                    record["min_evalue"],
                    evalue_value,
                )
                record["max_bitscore"] = max(
                    record["max_bitscore"],
                    bitscore_value,
                )

                totals[f"{read_label}_hits"] += 1
                totals["all_hits"] += 1
                totals["sum_pident"] += pident_value
                totals["max_bitscore"] = max(
                    totals["max_bitscore"],
                    bitscore_value,
                )
                totals["min_evalue"] = min(
                    totals["min_evalue"],
                    evalue_value,
                )

    process(r1_path, "R1")
    process(r2_path, "R2")

    hit_destination.parent.mkdir(parents=True, exist_ok=True)

    with hit_destination.open("w", newline="") as handle:
        columns = [
            "card_subject_id",
            "total_read_hits",
            "R1_read_hits",
            "R2_read_hits",
            "mean_percent_identity",
            "maximum_alignment_length",
            "minimum_evalue",
            "maximum_bitscore",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for subject in sorted(
            aggregated,
            key=lambda key: (
                -aggregated[key]["total"],
                key,
            ),
        ):
            record = aggregated[subject]

            writer.writerow(
                {
                    "card_subject_id": subject,
                    "total_read_hits": int(record["total"]),
                    "R1_read_hits": int(record["r1"]),
                    "R2_read_hits": int(record["r2"]),
                    "mean_percent_identity": (
                        f"{record['sum_pident'] / record['total']:.4f}"
                    ),
                    "maximum_alignment_length": int(
                        record["max_length"]
                    ),
                    "minimum_evalue": (
                        f"{record['min_evalue']:.6g}"
                    ),
                    "maximum_bitscore": (
                        f"{record['max_bitscore']:.4f}"
                    ),
                }
            )

    metrics = [
        ("sample_id", sample),
        ("R1_hits", str(totals["R1_hits"])),
        ("R2_hits", str(totals["R2_hits"])),
        ("total_hits", str(totals["all_hits"])),
        ("unique_CARD_subjects", str(len(aggregated))),
        (
            "mean_percent_identity",
            (
                f"{totals['sum_pident'] / totals['all_hits']:.4f}"
                if totals["all_hits"]
                else "NA"
            ),
        ),
        (
            "minimum_evalue",
            (
                f"{totals['min_evalue']:.6g}"
                if totals["all_hits"]
                else "NA"
            ),
        ),
        (
            "maximum_bitscore",
            (
                f"{totals['max_bitscore']:.4f}"
                if totals["all_hits"]
                else "NA"
            ),
        ),
    ]

    with sample_destination.open("w", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(["metric", "value"])
        writer.writerows(metrics)


if TARGET.exists():
    raise SystemExit(
        f"Release target already exists; refusing to overwrite: {TARGET}"
    )

if not MEGA_SHEET.is_file():
    raise SystemExit(f"Missing MEGAHIT samplesheet: {MEGA_SHEET}")

if BUILD.exists():
    shutil.rmtree(BUILD)

BUILD.mkdir(parents=True)

with MEGA_SHEET.open(newline="") as handle:
    samplesheet_rows = list(
        csv.DictReader(handle, delimiter="\t")
    )

samples = [
    row["sample_id"].strip()
    for row in samplesheet_rows
]

if len(samples) != EXPECTED_SAMPLES:
    raise SystemExit(
        f"Expected {EXPECTED_SAMPLES} samples; found {len(samples)}"
    )

if len(set(samples)) != EXPECTED_SAMPLES:
    raise SystemExit("Duplicate sample IDs detected")

mega_rows = {
    row["sample_id"].strip(): row
    for row in samplesheet_rows
}

categories = [
    "CARD_RGI_BWT",
    "CARD_RGI_Main",
    "geNomAD",
    "Assembly",
    "FASTP",
    "Bowtie2",
    "DIAMOND_CARD",
    "ARG_Mobility_Overlap",
    "metadata",
]

for category in categories:
    (BUILD / category).mkdir(parents=True, exist_ok=True)

manifest_rows: list[dict[str, str]] = []

for index, sample in enumerate(samples, start=1):
    print(f"[{index:02d}/{EXPECTED_SAMPLES}] Packaging {sample}")

    # CARD/RGI BWT
    bwt_stats = find_rgi_bwt_file(
        sample,
        "overall_mapping_stats.txt",
    )

    bwt_genes = find_rgi_bwt_file(
        sample,
        "gene_mapping_data.txt",
    )

    copy_normalized_rgi_bwt_file(
        bwt_stats,
        Path("CARD_RGI_BWT")
        / sample
        / f"{sample}_rgi_bwt.overall_mapping_stats.txt",
    )

    copy_normalized_rgi_bwt_file(
        bwt_genes,
        Path("CARD_RGI_BWT")
        / sample
        / f"{sample}_rgi_bwt.gene_mapping_data.txt",
    )

    # CARD/RGI Main compact table
    rgi_main_source = (
        RGI_MAIN_ROOT
        / sample
        / f"{sample}_rgi_main_contigs.txt"
    )

    require_file(rgi_main_source)

    rgi_main_relative = (
        Path("CARD_RGI_Main")
        / sample
        / f"{sample}_rgi_main.txt"
    )

    compact_rgi_main(
        rgi_main_source,
        destination_path(rgi_main_relative),
    )

    add_provenance(
        rgi_main_relative,
        [rgi_main_source],
        (
            "Tabular RGI Main output copied after removing "
            "Predicted_DNA, Predicted_Protein, and "
            "CARD_Protein_Sequence columns"
        ),
    )

    # geNomAD aggregated contig classification
    genomad_source = (
        GENOMAD_ROOT
        / sample
        / "final.contigs_aggregated_classification"
        / "final.contigs_aggregated_classification.tsv"
    )

    copy_release_file(
        genomad_source,
        Path("geNomAD")
        / sample
        / "contig_classification.tsv",
        (
            "geNomAD aggregated contig classification copied "
            "and standardized as contig_classification.tsv"
        ),
    )

    # Assembly metrics
    mega_outdir = Path(
        mega_rows[sample]["outdir"].strip()
    )

    contigs = mega_outdir / f"{sample}.contigs.fa"

    if not contigs.is_file():
        contigs = mega_outdir / "final.contigs.fa"

    require_file(contigs)

    assembly_relative = (
        Path("Assembly")
        / sample
        / f"{sample}_assembly_metrics.tsv"
    )

    write_assembly_metrics(
        sample,
        contigs,
        destination_path(assembly_relative),
    )

    add_provenance(
        assembly_relative,
        [contigs],
        (
            "Assembly statistics calculated from final MEGAHIT "
            "contigs; full contig FASTA excluded"
        ),
    )

    # FASTP
    fastp_source = (
        FASTP_ROOT
        / sample
        / f"{sample}.fastp.json"
    )

    copy_release_file(
        fastp_source,
        Path("FASTP")
        / sample
        / f"{sample}.fastp.json",
    )

    # Bowtie2
    bowtie_source = (
        BOWTIE2_ROOT
        / sample
        / f"{sample}.bowtie2.log"
    )

    copy_release_file(
        bowtie_source,
        Path("Bowtie2")
        / sample
        / f"{sample}.bowtie2.log",
    )

    # DIAMOND/CARD compact summaries
    diamond_sample_root = DIAMOND_ROOT / sample

    r1_diamond = (
        diamond_sample_root
        / f"{sample}_R1_diamond_card.tsv"
    )

    r2_diamond = (
        diamond_sample_root
        / f"{sample}_R2_diamond_card.tsv"
    )

    hit_relative = (
        Path("DIAMOND_CARD")
        / sample
        / f"{sample}_diamond_card_hit_summary.tsv"
    )

    sample_relative = (
        Path("DIAMOND_CARD")
        / sample
        / f"{sample}_diamond_card_sample_summary.tsv"
    )

    summarize_diamond(
        sample,
        r1_diamond,
        r2_diamond,
        destination_path(hit_relative),
        destination_path(sample_relative),
    )

    diamond_transformation = (
        "Compact aggregate summaries generated from six-column "
        "R1 and R2 DIAMOND/CARD BLASTX tables; raw per-read "
        "tables excluded"
    )

    add_provenance(
        hit_relative,
        [r1_diamond, r2_diamond],
        diamond_transformation,
    )

    add_provenance(
        sample_relative,
        [r1_diamond, r2_diamond],
        diamond_transformation,
    )

    manifest_rows.append(
        {
            "sample_id": sample,
            "cohort": "RM2016_USA_crossstudy_additional8",
            "FASTP": "PASS",
            "Bowtie2": "PASS",
            "CARD_RGI_BWT": "PASS",
            "DIAMOND_CARD": "PASS",
            "Assembly": "PASS",
            "geNomAD": "PASS",
            "CARD_RGI_Main": "PASS",
            "ARG_Mobility_Overlap": "PASS",
        }
    )



# ARG-mobility overlap
arg_mobility_files = [
    f"{ARG_MOBILITY_PREFIX}_ARG_mobility_overlap.tsv",
    f"{ARG_MOBILITY_PREFIX}_ARG_mobility_overlap_strict_perfect.tsv",
    f"{ARG_MOBILITY_PREFIX}_ARG_mobility_sample_summary.tsv",
    f"{ARG_MOBILITY_PREFIX}_ARG_mobility_sample_summary_strict_perfect.tsv",
    f"{ARG_MOBILITY_PREFIX}_ARG_mobility_drug_class_summary.tsv",
]

for filename in arg_mobility_files:
    source = ARG_MOBILITY_ROOT / filename

    copy_release_file(
        source,
        Path("ARG_Mobility_Overlap") / filename,
        (
            "Validated cohort-level overlap linking contig-level "
            "CARD/RGI Main ARG calls with geNomAD plasmid/virus context"
        ),
    )

manifest_relative = (
    Path("metadata")
    / "USA_crossstudy_additional8_sample_manifest.tsv"
)

manifest_path = destination_path(manifest_relative)

manifest_columns = [
    "sample_id",
    "cohort",
    "FASTP",
    "Bowtie2",
    "CARD_RGI_BWT",
    "DIAMOND_CARD",
    "Assembly",
    "geNomAD",
    "CARD_RGI_Main",
    "ARG_Mobility_Overlap",
]

with manifest_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=manifest_columns,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(manifest_rows)

add_provenance(
    manifest_relative,
    [MEGA_SHEET],
    "Standardized 10-sample release manifest generated",
)


readme = f"""# ResistanceRadar RM2016 USA cross-study additional8 compact dataset inputs

This directory is the standardized, GitHub-safe release for the
RM2016 USA cross-study additional8 cohort.

## Cohort

- Samples: {EXPECTED_SAMPLES}
- Sample IDs: {", ".join(samples)}
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
"""

readme_path = BUILD / "README.md"
readme_path.write_text(readme)

add_provenance(
    Path("README.md"),
    [MEGA_SHEET],
    "Release documentation generated",
)


provenance_path = BUILD / "PROVENANCE.tsv"

with provenance_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "destination",
            "source",
            "transformation",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(
        sorted(
            provenance,
            key=lambda row: row["destination"],
        )
    )


# Build file index before checksums. The index intentionally excludes
# itself and SHA256SUMS.txt.
index_path = BUILD / "FILE_INDEX.tsv"

indexed_files = [
    path
    for path in BUILD.rglob("*")
    if path.is_file()
    and path.name not in {
        "FILE_INDEX.tsv",
        "SHA256SUMS.txt",
    }
]

with index_path.open("w", newline="") as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(
        [
            "relative_path",
            "category",
            "sample_id",
            "bytes",
        ]
    )

    for path in sorted(indexed_files):
        relative = path.relative_to(BUILD)
        parts = relative.parts
        category = parts[0]
        sample_id = (
            parts[1]
            if len(parts) >= 3
            and category
            in {
                "CARD_RGI_BWT",
                "CARD_RGI_Main",
                "geNomAD",
                "Assembly",
                "FASTP",
                "Bowtie2",
                "DIAMOND_CARD",
            }
            else ""
        )

        writer.writerow(
            [
                str(relative),
                category,
                sample_id,
                path.stat().st_size,
            ]
        )


checksum_path = BUILD / "SHA256SUMS.txt"

checksum_files = [
    path
    for path in BUILD.rglob("*")
    if path.is_file()
    and path.name != "SHA256SUMS.txt"
]

with checksum_path.open("w") as handle:
    for path in sorted(checksum_files):
        relative = path.relative_to(BUILD)
        handle.write(f"{sha256(path)}  {relative}\n")


expected_counts = {
    "CARD_RGI_BWT": 2 * EXPECTED_SAMPLES,
    "CARD_RGI_Main": EXPECTED_SAMPLES,
    "geNomAD": EXPECTED_SAMPLES,
    "Assembly": EXPECTED_SAMPLES,
    "FASTP": EXPECTED_SAMPLES,
    "Bowtie2": EXPECTED_SAMPLES,
    "DIAMOND_CARD": 2 * EXPECTED_SAMPLES,
    "ARG_Mobility_Overlap": 5,
    "metadata": 1,
}

print("\n===== RELEASE VALIDATION =====")

for category, expected in expected_counts.items():
    observed = sum(
        1
        for path in (BUILD / category).rglob("*")
        if path.is_file()
    )

    print(f"{category:18s} {observed:3d}/{expected}")

    if observed != expected:
        raise SystemExit(
            f"{category}: expected {expected}, found {observed}"
        )


all_files = [
    path
    for path in BUILD.rglob("*")
    if path.is_file()
]

symlinks = [
    path
    for path in BUILD.rglob("*")
    if path.is_symlink()
]

empty_files = [
    path
    for path in all_files
    if path.stat().st_size == 0
]

oversized = [
    path
    for path in all_files
    if path.stat().st_size >= GITHUB_LIMIT
]

if symlinks:
    raise SystemExit(
        "Symlinks found in release:\n"
        + "\n".join(str(path) for path in symlinks)
    )

if empty_files:
    raise SystemExit(
        "Empty files found in release:\n"
        + "\n".join(str(path) for path in empty_files)
    )

if oversized:
    raise SystemExit(
        "Files exceed GitHub's 100 MiB limit:\n"
        + "\n".join(str(path) for path in oversized)
    )

checksum_lines = [
    line
    for line in checksum_path.read_text().splitlines()
    if line.strip()
]

if len(checksum_lines) != len(all_files) - 1:
    raise SystemExit(
        "Checksum count does not match released file count"
    )

largest = max(all_files, key=lambda path: path.stat().st_size)
total_bytes = sum(path.stat().st_size for path in all_files)

print(f"\nTotal files:       {len(all_files)}")
print(f"Symlinks:          {len(symlinks)}")
print(f"Empty files:       {len(empty_files)}")
print(f"Total size:        {total_bytes / 1048576:.2f} MiB")
print(
    "Largest file:     "
    f"{largest.stat().st_size / 1048576:.2f} MiB "
    f"({largest.relative_to(BUILD)})"
)
print(f"Checksum entries:  {len(checksum_lines)}")

BUILD.rename(TARGET)

print("\n===== RELEASE BUILD COMPLETE — PASS =====")
print(f"Release: {TARGET.relative_to(REPO)}")
