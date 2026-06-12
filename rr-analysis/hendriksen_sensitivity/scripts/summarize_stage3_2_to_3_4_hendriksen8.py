#!/usr/bin/env python3
#@rgoutam
from pathlib import Path
import json
import re
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("rr-analysis/hendriksen_sensitivity")

QC_DIR = BASE / "fastp_bowtie2_8sample" / "qc_reports"
BT_DIR = BASE / "fastp_bowtie2_8sample" / "bowtie2"
RGI_DIR = BASE / "rgi_bwt"

TABLE_DIR = BASE / "tables"
FIG_DIR = BASE / "figures"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def sample_from_fastp_json(path):
    return path.name.replace("_fastp.json", "")

def sample_from_host_log(path):
    return path.name.replace("_host_removal_log.tsv", "")

def sample_from_rgi_gene(path):
    return path.name.replace("_rgi_bwt.gene_mapping_data.txt", "")

def sample_from_rgi_stats(path):
    return path.name.replace("_rgi_bwt.overall_mapping_stats.txt", "")

def to_num(x):
    try:
        return pd.to_numeric(str(x).replace(",", "").replace("%", ""), errors="coerce")
    except Exception:
        return None

def find_col(df, candidates, contains=None):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    if contains:
        for c in df.columns:
            cl = c.lower()
            if all(k.lower() in cl for k in contains):
                return c
    return None

# ------------------------------------------------------------
# Stage 3.2 FASTP summary
# ------------------------------------------------------------

fastp_rows = []

for f in sorted(QC_DIR.glob("*_fastp.json")):
    sample = sample_from_fastp_json(f)

    with open(f, "r") as handle:
        js = json.load(handle)

    summary = js.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})

    row = {
        "sample_id": sample,
        "stage3_2_fastp_json": str(f),
        "raw_total_reads": before.get("total_reads"),
        "raw_total_bases": before.get("total_bases"),
        "raw_q20_rate": before.get("q20_rate"),
        "raw_q30_rate": before.get("q30_rate"),
        "raw_gc_content": before.get("gc_content"),
        "filtered_total_reads": after.get("total_reads"),
        "filtered_total_bases": after.get("total_bases"),
        "filtered_q20_rate": after.get("q20_rate"),
        "filtered_q30_rate": after.get("q30_rate"),
        "filtered_gc_content": after.get("gc_content"),
        "duplication_rate": js.get("duplication", {}).get("rate"),
        "adapter_trimmed_reads": js.get("adapter_cutting", {}).get("adapter_trimmed_reads"),
    }

    if row["raw_total_reads"] and row["filtered_total_reads"]:
        row["fastp_read_retention_percent"] = 100 * row["filtered_total_reads"] / row["raw_total_reads"]
    else:
        row["fastp_read_retention_percent"] = None

    fastp_rows.append(row)

fastp_df = pd.DataFrame(fastp_rows)

# ------------------------------------------------------------
# Stage 3.3 Bowtie2 host-removal summary
# ------------------------------------------------------------

host_logs = sorted(BT_DIR.glob("*_host_removal_log.tsv"))
host_tables = []

for f in host_logs:
    try:
        df = pd.read_csv(f, sep="\t")
        host_tables.append(df)
    except Exception as e:
        print(f"Could not read {f}: {e}")

if host_tables:
    host_df = pd.concat(host_tables, ignore_index=True)
else:
    host_df = pd.DataFrame()

# ------------------------------------------------------------
# Stage 3.4 RGI BWT summaries
# ------------------------------------------------------------

rgi_gene_files = sorted(RGI_DIR.glob("*_rgi_bwt.gene_mapping_data.txt"))
rgi_stat_files = sorted(RGI_DIR.glob("*_rgi_bwt.overall_mapping_stats.txt"))

gene_tables = []

for f in rgi_gene_files:
    sample = sample_from_rgi_gene(f)
    try:
        df = pd.read_csv(f, sep="\t")
        df["sample_id"] = sample
        gene_tables.append(df)
    except Exception as e:
        print(f"Could not read {f}: {e}")

if gene_tables:
    genes = pd.concat(gene_tables, ignore_index=True)
else:
    genes = pd.DataFrame()

rgi_stat_rows = []

for f in rgi_stat_files:
    sample = sample_from_rgi_stats(f)
    text = f.read_text(errors="ignore")
    row = {"sample_id": sample, "stage3_4_rgi_stats_file": str(f)}

    for line in text.splitlines():
        clean = line.strip()
        if not clean:
            continue

        if "\t" in clean:
            parts = [x.strip() for x in clean.split("\t") if x.strip()]
            if len(parts) >= 2:
                key = re.sub(r"[^A-Za-z0-9]+", "_", parts[0]).strip("_").lower()
                row[key] = parts[-1]

        elif ":" in clean:
            key, val = clean.split(":", 1)
            key = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
            row[key] = val.strip()

    rgi_stat_rows.append(row)

rgi_stats_df = pd.DataFrame(rgi_stat_rows)

# ------------------------------------------------------------
# Per-sample RGI ARG summary
# ------------------------------------------------------------

if not genes.empty:
    gene_col = find_col(
        genes,
        ["Best_Hit_ARO", "Best Hit ARO", "ARO Term", "ARO Name", "Gene"],
        contains=["aro"]
    )

    drug_col = find_col(
        genes,
        ["Drug Class", "drug_class", "ARO Drug Class"],
        contains=["drug", "class"]
    )

    mech_col = find_col(
        genes,
        ["Resistance Mechanism", "resistance_mechanism", "ARO Resistance Mechanism"],
        contains=["mechanism"]
    )

    mapped_col = find_col(
        genes,
        ["All Mapped Reads", "Mapped Reads", "mapped_reads", "read_count"],
        contains=["mapped", "read"]
    )

    per_rgi = genes.groupby("sample_id").size().reset_index(name="rgi_gene_mapping_rows")

    if gene_col:
        tmp = genes.groupby("sample_id")[gene_col].nunique().reset_index(name="unique_rgi_arg_names")
        per_rgi = per_rgi.merge(tmp, on="sample_id", how="left")

    if drug_col:
        dtmp = genes[["sample_id", drug_col]].copy()
        dtmp[drug_col] = dtmp[drug_col].astype(str).str.split(";")
        dtmp = dtmp.explode(drug_col)
        dtmp[drug_col] = dtmp[drug_col].astype(str).str.strip()
        dtmp = dtmp[dtmp[drug_col].ne("")]
        tmp = dtmp.groupby("sample_id")[drug_col].nunique().reset_index(name="unique_rgi_drug_classes")
        per_rgi = per_rgi.merge(tmp, on="sample_id", how="left")

    if mech_col:
        mtmp = genes[["sample_id", mech_col]].copy()
        mtmp[mech_col] = mtmp[mech_col].astype(str).str.split(";")
        mtmp = mtmp.explode(mech_col)
        mtmp[mech_col] = mtmp[mech_col].astype(str).str.strip()
        mtmp = mtmp[mtmp[mech_col].ne("")]
        tmp = mtmp.groupby("sample_id")[mech_col].nunique().reset_index(name="unique_resistance_mechanisms")
        per_rgi = per_rgi.merge(tmp, on="sample_id", how="left")
else:
    gene_col = drug_col = mech_col = mapped_col = None
    per_rgi = pd.DataFrame()

# ------------------------------------------------------------
# Combined Stage 3.2-3.4 summary table
# ------------------------------------------------------------

summary = fastp_df.copy()

if not host_df.empty:
    summary = summary.merge(host_df, on="sample_id", how="outer")

if not per_rgi.empty:
    summary = summary.merge(per_rgi, on="sample_id", how="outer")

summary_out = TABLE_DIR / "Stage3_2_to_3_4_Hendriksen8_pipeline_summary.tsv"
summary.to_csv(summary_out, sep="\t", index=False)

fastp_df.to_csv(TABLE_DIR / "Stage3_2_FASTP_summary.tsv", sep="\t", index=False)
if not host_df.empty:
    host_df.to_csv(TABLE_DIR / "Stage3_3_Bowtie2_host_removal_summary.tsv", sep="\t", index=False)
if not rgi_stats_df.empty:
    rgi_stats_df.to_csv(TABLE_DIR / "Stage3_4_RGI_overall_mapping_summary_raw.tsv", sep="\t", index=False)
if not genes.empty:
    genes.to_csv(TABLE_DIR / "Stage3_4_RGI_gene_mapping_combined.tsv", sep="\t", index=False)
if not per_rgi.empty:
    per_rgi.to_csv(TABLE_DIR / "Stage3_4_RGI_per_sample_ARG_summary.tsv", sep="\t", index=False)

# ------------------------------------------------------------
# Plot 1: pipeline read flow, Stage 3.2 to 3.3
# ------------------------------------------------------------

plot_df = summary.sort_values("sample_id").copy()

# Try to find host-removal clean/read columns
input_pairs_col = "input_read_pairs" if "input_read_pairs" in plot_df.columns else None
clean_pairs_col = "final_clean_read_pairs" if "final_clean_read_pairs" in plot_df.columns else None
removed_col = "percent_removed" if "percent_removed" in plot_df.columns else None

for c in ["raw_total_reads", "filtered_total_reads", input_pairs_col, clean_pairs_col, removed_col]:
    if c and c in plot_df.columns:
        plot_df[c] = pd.to_numeric(plot_df[c], errors="coerce")

if "raw_total_reads" in plot_df.columns and "filtered_total_reads" in plot_df.columns:
    plt.figure(figsize=(9, 5))
    x = range(len(plot_df))
    plt.plot(x, plot_df["raw_total_reads"] / 2, marker="o", label="Raw read pairs")
    plt.plot(x, plot_df["filtered_total_reads"] / 2, marker="o", label="FASTP-filtered read pairs")

    if clean_pairs_col and clean_pairs_col in plot_df.columns:
        plt.plot(x, plot_df[clean_pairs_col], marker="o", label="Bowtie2-clean read pairs")

    plt.xticks(x, plot_df["sample_id"], rotation=45, ha="right")
    plt.ylabel("Read pairs")
    plt.xlabel("Sample")
    plt.title("Stage 3.2–3.3: Read-pair flow through QC and host removal")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Stage3_2_to_3_3_read_pair_flow.png", dpi=300)
    plt.savefig(FIG_DIR / "Stage3_2_to_3_3_read_pair_flow.pdf")
    plt.close()

# ------------------------------------------------------------
# Plot 2: FASTP Q30 and retention
# ------------------------------------------------------------

if "filtered_q30_rate" in plot_df.columns and "fastp_read_retention_percent" in plot_df.columns:
    plot_df["filtered_q30_percent"] = pd.to_numeric(plot_df["filtered_q30_rate"], errors="coerce") * 100

    plt.figure(figsize=(9, 5))
    x = range(len(plot_df))
    plt.bar(x, plot_df["fastp_read_retention_percent"], label="Read retention (%)")
    plt.plot(x, plot_df["filtered_q30_percent"], marker="o", label="Filtered Q30 (%)")
    plt.xticks(x, plot_df["sample_id"], rotation=45, ha="right")
    plt.ylabel("Percent")
    plt.xlabel("Sample")
    plt.title("Stage 3.2 FASTP: Read retention and Q30")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Stage3_2_FASTP_retention_Q30.png", dpi=300)
    plt.savefig(FIG_DIR / "Stage3_2_FASTP_retention_Q30.pdf")
    plt.close()

# ------------------------------------------------------------
# Plot 3: Bowtie2 host/contaminant removal percent
# ------------------------------------------------------------

if removed_col and removed_col in plot_df.columns:
    plt.figure(figsize=(9, 5))
    plt.bar(plot_df["sample_id"], plot_df[removed_col])
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Percent removed")
    plt.xlabel("Sample")
    plt.title("Stage 3.3 Bowtie2: Host/contaminant removal rate")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Stage3_3_Bowtie2_percent_removed.png", dpi=300)
    plt.savefig(FIG_DIR / "Stage3_3_Bowtie2_percent_removed.pdf")
    plt.close()

# ------------------------------------------------------------
# Plot 4: RGI ARG rows per sample
# ------------------------------------------------------------

if "rgi_gene_mapping_rows" in plot_df.columns:
    plt.figure(figsize=(9, 5))
    plt.bar(plot_df["sample_id"], plot_df["rgi_gene_mapping_rows"])
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("RGI gene-mapping rows")
    plt.xlabel("Sample")
    plt.title("Stage 3.4 CARD/RGI: ARG hit rows per sample")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Stage3_4_RGI_ARG_hit_rows_per_sample.png", dpi=300)
    plt.savefig(FIG_DIR / "Stage3_4_RGI_ARG_hit_rows_per_sample.pdf")
    plt.close()

# ------------------------------------------------------------
# Plot 5: Drug class summary and sample x drug-class heatmap
# ------------------------------------------------------------

if not genes.empty and drug_col:
    drug = genes.copy()
    drug[drug_col] = drug[drug_col].astype(str).str.split(";")
    drug = drug.explode(drug_col)
    drug[drug_col] = drug[drug_col].astype(str).str.strip()
    drug = drug[drug[drug_col].ne("")]

    if mapped_col and mapped_col in drug.columns:
        drug[mapped_col] = pd.to_numeric(drug[mapped_col], errors="coerce").fillna(0)
        sample_drug = drug.groupby(["sample_id", drug_col])[mapped_col].sum().reset_index()
        value_col = mapped_col
    else:
        sample_drug = drug.groupby(["sample_id", drug_col]).size().reset_index(name="n_hits")
        value_col = "n_hits"

    sample_drug.to_csv(TABLE_DIR / "Stage3_4_RGI_sample_by_drug_class.tsv", sep="\t", index=False)

    top_drug = sample_drug.groupby(drug_col)[value_col].sum().sort_values(ascending=False).head(15).reset_index()

    plt.figure(figsize=(9, 6))
    plt.barh(top_drug[drug_col][::-1], top_drug[value_col][::-1])
    plt.xlabel(value_col)
    plt.ylabel("Drug class")
    plt.title("Stage 3.4 CARD/RGI: Top drug classes across Hendriksen 8-sample pilot")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Stage3_4_RGI_top_drug_classes.png", dpi=300)
    plt.savefig(FIG_DIR / "Stage3_4_RGI_top_drug_classes.pdf")
    plt.close()

    matrix = sample_drug.pivot_table(
        index="sample_id",
        columns=drug_col,
        values=value_col,
        aggfunc="sum",
        fill_value=0
    )

    top_cols = matrix.sum(axis=0).sort_values(ascending=False).head(15).index
    matrix_top = matrix[top_cols]
    matrix_top.to_csv(TABLE_DIR / "Stage3_4_RGI_sample_by_top_drug_class_matrix.tsv", sep="\t")

    plt.figure(figsize=(11, 5.5))
    plt.imshow(matrix_top.values, aspect="auto")
    plt.xticks(range(len(matrix_top.columns)), matrix_top.columns, rotation=45, ha="right")
    plt.yticks(range(len(matrix_top.index)), matrix_top.index)
    plt.colorbar(label=value_col)
    plt.title("Stage 3.4 CARD/RGI: Sample × drug-class signal")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Stage3_4_RGI_sample_by_drug_class_heatmap.png", dpi=300)
    plt.savefig(FIG_DIR / "Stage3_4_RGI_sample_by_drug_class_heatmap.pdf")
    plt.close()

# ------------------------------------------------------------
# Plot 6: Resistance mechanism summary
# ------------------------------------------------------------

if not genes.empty and mech_col:
    mech = genes.copy()
    mech[mech_col] = mech[mech_col].astype(str).str.split(";")
    mech = mech.explode(mech_col)
    mech[mech_col] = mech[mech_col].astype(str).str.strip()
    mech = mech[mech[mech_col].ne("")]

    mech_summary = mech.groupby(mech_col).size().sort_values(ascending=False).reset_index(name="n_hits")
    mech_summary.to_csv(TABLE_DIR / "Stage3_4_RGI_resistance_mechanism_summary.tsv", sep="\t", index=False)

    top_mech = mech_summary.head(12)

    plt.figure(figsize=(9, 5.5))
    plt.barh(top_mech[mech_col][::-1], top_mech["n_hits"][::-1])
    plt.xlabel("Number of RGI rows")
    plt.ylabel("Resistance mechanism")
    plt.title("Stage 3.4 CARD/RGI: Resistance mechanisms")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "Stage3_4_RGI_resistance_mechanisms.png", dpi=300)
    plt.savefig(FIG_DIR / "Stage3_4_RGI_resistance_mechanisms.pdf")
    plt.close()

print("\nDetected RGI columns:")
print("gene_col:", gene_col)
print("drug_col:", drug_col)
print("mechanism_col:", mech_col)
print("mapped_col:", mapped_col)

print("\nWrote summary table:")
print(summary_out)

print("\nTables:")
for p in sorted(TABLE_DIR.glob("Stage3_*")):
    print(" ", p)

print("\nFigures:")
for p in sorted(FIG_DIR.glob("Stage3_*")):
    print(" ", p)

print("\nDone.")
