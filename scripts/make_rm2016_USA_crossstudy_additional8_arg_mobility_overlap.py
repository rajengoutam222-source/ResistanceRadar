#!/usr/bin/env python3

from pathlib import Path
import gzip
import pandas as pd

ROOT = Path("rr-analysis/hendriksen_resistancemap_2016")

RGI_DIR = ROOT / "stage3_5_rgi_main_contigs_USA_crossstudy_additional8/results"
GENOMAD_DIR = ROOT / "stage3_5_genomad_USA_crossstudy_additional8"

OUTDIR = ROOT / "stage3_5_arg_mobility_overlap_USA_crossstudy_additional8"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUTDIR / "rm2016_USA_crossstudy_additional8_ARG_mobility_overlap.tsv"
OUT_ALL_GZ = OUTDIR / "rm2016_USA_crossstudy_additional8_ARG_mobility_overlap.tsv.gz"
OUT_STRICT = OUTDIR / "rm2016_USA_crossstudy_additional8_ARG_mobility_overlap_strict_perfect.tsv"
OUT_SUMMARY = OUTDIR / "rm2016_USA_crossstudy_additional8_ARG_mobility_sample_summary.tsv"
OUT_SUMMARY_STRICT = OUTDIR / "rm2016_USA_crossstudy_additional8_ARG_mobility_sample_summary_strict_perfect.tsv"
OUT_DRUG = OUTDIR / "rm2016_USA_crossstudy_additional8_ARG_mobility_drug_class_summary.tsv"
OUT_CKPT = ROOT / "checkpoints/RM2016_USA_CROSSSTUDY_ADDITIONAL8_ARG_MOBILITY_OVERLAP_COMPLETE.txt"

SAMPLES = [
    "ERR1713397",
    "ERR1713398",
    "ERR1713400",
    "ERR1713405",
    "ERR1713406",
    "ERR2729809",
    "SRR5007116",
    "SRR5007352",
]


def find_rgi_txt(sample: str) -> Path:
    candidates = list((RGI_DIR / sample).glob(f"{sample}_rgi_main_contigs.txt"))
    candidates += list((RGI_DIR / sample).glob("*_rgi_main_contigs.txt"))
    candidates = [p for p in candidates if p.is_file() and p.stat().st_size > 0]
    if not candidates:
        raise FileNotFoundError(f"No RGI-main contig txt found for {sample}")
    return candidates[0]


def find_genomad_summary(sample: str, kind: str) -> Path | None:
    # kind = plasmid or virus
    sample_dir = GENOMAD_DIR / "results" / sample
    candidates = []
    if sample_dir.exists():
        candidates += list(sample_dir.glob(f"final.contigs_summary/final.contigs_{kind}_summary.tsv"))
        candidates += list(sample_dir.glob(f"{sample}.contigs_summary/{sample}.contigs_{kind}_summary.tsv"))
        candidates += list(sample_dir.glob(f"*.contigs_summary/*.contigs_{kind}_summary.tsv"))

    # fallback: search whole stage if layout differs
    candidates += list(GENOMAD_DIR.glob(f"**/{sample}/**/*_{kind}_summary.tsv"))
    candidates += list(GENOMAD_DIR.glob(f"**/*{sample}*_{kind}_summary.tsv"))

    candidates = [p for p in candidates if p.is_file() and p.stat().st_size > 0]
    if candidates:
        return sorted(set(candidates))[0]
    return None


def read_mobile_tables(sample: str):
    plasmid_file = find_genomad_summary(sample, "plasmid")
    virus_file = find_genomad_summary(sample, "virus")

    plasmids = {}
    viruses = {}

    if plasmid_file is not None:
        p = pd.read_csv(plasmid_file, sep="\t")
        if "seq_name" in p.columns:
            for _, row in p.iterrows():
                seq = str(row["seq_name"])
                plasmids[seq] = row["plasmid_score"] if "plasmid_score" in p.columns else ""

    if virus_file is not None:
        v = pd.read_csv(virus_file, sep="\t")
        if "seq_name" in v.columns:
            for _, row in v.iterrows():
                seq = str(row["seq_name"])
                viruses[seq] = row["virus_score"] if "virus_score" in v.columns else ""

    return plasmids, viruses, plasmid_file, virus_file


rows = []
missing = []

for sample in SAMPLES:
    rgi_file = find_rgi_txt(sample)
    rgi = pd.read_csv(rgi_file, sep="\t")

    plasmids, viruses, plasmid_file, virus_file = read_mobile_tables(sample)

    if plasmid_file is None:
        missing.append((sample, "missing plasmid summary"))
    if virus_file is None:
        missing.append((sample, "missing virus summary"))

    for _, row in rgi.iterrows():
        contig = str(row.get("Contig", ""))
        in_plasmid = contig in plasmids
        in_virus = contig in viruses

        if in_plasmid and in_virus:
            context = "plasmid_like;virus_like"
        elif in_plasmid:
            context = "plasmid_like"
        elif in_virus:
            context = "virus_like"
        else:
            context = "chromosome_like"

        rows.append({
            "set": "RM2016_USA_crossstudy_additional8",
            "sample_id": sample,
            "contig_id": contig,
            "ORF_ID": row.get("ORF_ID", ""),
            "Best_Hit_ARO": row.get("Best_Hit_ARO", ""),
            "ARO": row.get("ARO", ""),
            "Cut_Off": row.get("Cut_Off", ""),
            "Best_Identities": row.get("Best_Identities", ""),
            "Percentage_Length_of_Reference_Sequence": row.get("Percentage Length of Reference Sequence", ""),
            "Drug_Class": row.get("Drug Class", ""),
            "Resistance_Mechanism": row.get("Resistance Mechanism", ""),
            "AMR_Gene_Family": row.get("AMR Gene Family", ""),
            "Antibiotic": row.get("Antibiotic", ""),
            "matched_to_genomad": bool(in_plasmid or in_virus),
            "genomad_context": context,
            "chromosome_score": "" if (in_plasmid or in_virus) else "1",
            "plasmid_score": plasmids.get(contig, ""),
            "virus_score": viruses.get(contig, ""),
        })

overlap = pd.DataFrame(rows)

overlap.to_csv(OUT_ALL, sep="\t", index=False)
with gzip.open(OUT_ALL_GZ, "wt") as gz:
    overlap.to_csv(gz, sep="\t", index=False)

strict = overlap[overlap["Cut_Off"].isin(["Strict", "Perfect"])].copy()
strict.to_csv(OUT_STRICT, sep="\t", index=False)


def summarize(df: pd.DataFrame, strict_only: bool = False) -> pd.DataFrame:
    out = []
    for sample, g in df.groupby("sample_id", sort=False):
        total = len(g)
        plasmid = g["genomad_context"].str.contains("plasmid_like", na=False).sum()
        virus = g["genomad_context"].str.contains("virus_like", na=False).sum()
        mobile = (
            g["genomad_context"].str.contains("plasmid_like", na=False)
            | g["genomad_context"].str.contains("virus_like", na=False)
        ).sum()
        chrom = (g["genomad_context"] == "chromosome_like").sum()
        matched = g["matched_to_genomad"].sum()
        strict_n = (g["Cut_Off"] == "Strict").sum()
        perfect_n = (g["Cut_Off"] == "Perfect").sum()
        loose_n = (g["Cut_Off"] == "Loose").sum()

        rec = {
            "set": "RM2016_USA_crossstudy_additional8",
            "sample_id": sample,
        }

        if strict_only:
            rec["strict_perfect_ARG_hits"] = total
        else:
            rec["total_ARG_hits"] = total
            rec["matched_to_genomad"] = matched

        rec.update({
            "chromosome_like_ARG_hits": chrom,
            "plasmid_like_ARG_hits": plasmid,
            "virus_like_ARG_hits": virus,
            "mobile_context_ARG_hits": mobile,
            "unclassified_ARG_hits": 0,
            "percent_ARG_hits_plasmid_like": (100 * plasmid / total) if total else 0,
            "percent_ARG_hits_mobile_context": (100 * mobile / total) if total else 0,
        })

        if strict_only:
            rec["strict_ARG_hits"] = strict_n
            rec["perfect_ARG_hits"] = perfect_n
        else:
            rec["loose_ARG_hits"] = loose_n
            rec["strict_ARG_hits"] = strict_n
            rec["perfect_ARG_hits"] = perfect_n

        out.append(rec)

    return pd.DataFrame(out)


summary = summarize(overlap, strict_only=False)
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

summary_strict = summarize(strict, strict_only=True)
summary_strict.to_csv(OUT_SUMMARY_STRICT, sep="\t", index=False)

drug_rows = []
for sample, g in overlap.groupby("sample_id", sort=False):
    for _, row in g.iterrows():
        classes = str(row.get("Drug_Class", "")).split(";")
        for cls in classes:
            cls = cls.strip()
            if not cls or cls.lower() == "nan":
                continue
            drug_rows.append({
                "set": "RM2016_USA_crossstudy_additional8",
                "sample_id": sample,
                "Drug_Class": cls,
                "ARG_hits": 1,
                "mobile_context": bool(row["matched_to_genomad"]),
                "Cut_Off": row.get("Cut_Off", ""),
            })

drug_df = pd.DataFrame(drug_rows)
if not drug_df.empty:
    drug_summary = (
        drug_df.groupby(["set", "sample_id", "Drug_Class"], as_index=False)
        .agg(
            ARG_hits=("ARG_hits", "sum"),
            mobile_context_ARG_hits=("mobile_context", "sum"),
        )
    )
else:
    drug_summary = pd.DataFrame(columns=["set", "sample_id", "Drug_Class", "ARG_hits", "mobile_context_ARG_hits"])

drug_summary.to_csv(OUT_DRUG, sep="\t", index=False)

with open(OUT_CKPT, "w") as f:
    f.write("RM2016 clean21 plus USA Chicago 11 ARG x mobility overlap complete.\n")
    f.write(f"Samples: {len(SAMPLES)}\n")
    f.write(f"Full overlap rows: {len(overlap)}\n")
    f.write(f"Strict+Perfect rows: {len(strict)}\n")
    f.write(f"Output directory: {OUTDIR}\n")
    if missing:
        f.write("Warnings:\n")
        for sample, msg in missing:
            f.write(f"- {sample}: {msg}\n")

print("DONE")
print(f"Samples: {len(SAMPLES)}")
print(f"Full overlap rows: {len(overlap)}")
print(f"Strict+Perfect rows: {len(strict)}")
print(f"Output directory: {OUTDIR}")
if missing:
    print("Warnings:")
    for sample, msg in missing:
        print(sample, msg)
