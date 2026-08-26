#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import csv
import re

ROOT = Path("rr-analysis")
OUTDIR = ROOT / "projectwide_ml_features"

ARO_IN = OUTDIR / "ResistanceRadar_projectwide_ARO_features_long.tsv"

GENE_LONG = OUTDIR / "ResistanceRadar_101_endpoint_matched_gene_features_long.tsv"
GENE_WIDE = OUTDIR / "ResistanceRadar_101_endpoint_matched_gene_features.tsv"
PANEL_OUT = OUTDIR / "ResistanceRadar_endpoint_gene_panel_preregistered.tsv"

QRDR_LONG = OUTDIR / "ResistanceRadar_101_Ecoli_QRDR_mutations_long.tsv"
QRDR_WIDE = OUTDIR / "ResistanceRadar_101_Ecoli_QRDR_mutations.tsv"
QRDR_AUDIT = OUTDIR / "ResistanceRadar_101_Ecoli_QRDR_source_audit.tsv"

FINAL_OUT = OUTDIR / "ResistanceRadar_101_endpoint_matched_AMR_features_FINAL.tsv"
DICT_OUT = OUTDIR / "ResistanceRadar_endpoint_feature_dictionary_FINAL.tsv"


# ==============================================================
# 1. LOAD PROJECT-WIDE CARD/RGI ARO FEATURES
# ==============================================================

df = pd.read_csv(ARO_IN, sep="\t", low_memory=False)

required = {
    "sample",
    "Best_Hit_ARO",
    "ARG_loci",
    "mobile_ARG_loci",
    "plasmid_ARG_loci",
}

missing = required - set(df.columns)

if missing:
    raise RuntimeError(f"Missing required ARO columns: {sorted(missing)}")

samples = sorted(df["sample"].astype(str).unique())

if len(samples) != 101:
    raise RuntimeError(f"Expected 101 samples, found {len(samples)}")


# ==============================================================
# 2. PREDEFINED ENDPOINT-MATCHED GENE PANEL
# ==============================================================

OXA48_LIKE = {
    "OXA-48",
    "OXA-162",
    "OXA-181",
    "OXA-204",
    "OXA-232",
    "OXA-244",
    "OXA-245",
    "OXA-370",
    "OXA-436",
    "OXA-438",
    "OXA-484",
    "OXA-519",
}

panel = [
    ["Fluoroquinolone", "qnr",
     "Qnr family", "existing_CARD_RGI"],

    ["Fluoroquinolone", "oqxA",
     "oqxA", "existing_CARD_RGI"],

    ["Fluoroquinolone", "oqxB",
     "oqxB", "existing_CARD_RGI"],

    ["Fluoroquinolone", "qepA",
     "qepA", "existing_CARD_RGI"],

    ["Fluoroquinolone", "aac6_Ib_cr",
     "AAC(6')-Ib-cr family", "existing_CARD_RGI"],

    ["Fluoroquinolone", "Ecoli_gyrA_QRDR",
     "E. coli gyrA resistance-associated QRDR mutation",
     "CARD_RGI_MAIN_StrictPerfect_protein_variant"],

    ["Fluoroquinolone", "Ecoli_parC_QRDR",
     "E. coli parC resistance-associated QRDR mutation",
     "CARD_RGI_MAIN_StrictPerfect_protein_variant"],

    ["Carbapenem", "KPC",
     "KPC carbapenemase family", "existing_CARD_RGI"],

    ["Carbapenem", "NDM",
     "NDM carbapenemase family", "existing_CARD_RGI"],

    ["Carbapenem", "VIM",
     "VIM carbapenemase family", "existing_CARD_RGI"],

    ["Carbapenem", "IMP",
     "IMP carbapenemase family", "existing_CARD_RGI"],

    ["Carbapenem", "OXA48_like",
     "OXA-48-like carbapenemase family", "existing_CARD_RGI"],
]

pd.DataFrame(
    panel,
    columns=["endpoint", "feature", "definition", "status"]
).to_csv(PANEL_OUT, sep="\t", index=False)


# ==============================================================
# 3. CLASSIFY PMQR + CARBAPENEMASE ARO NAMES
# ==============================================================

oxa48_upper = {x.upper() for x in OXA48_LIKE}


def classify_aro(name):
    n = str(name).strip()

    if re.fullmatch(r"Qnr[A-Za-z0-9_-]+", n, flags=re.I):
        return "Fluoroquinolone", "qnr"

    if re.fullmatch(r"oqxA", n, flags=re.I):
        return "Fluoroquinolone", "oqxA"

    if re.fullmatch(r"oqxB", n, flags=re.I):
        return "Fluoroquinolone", "oqxB"

    if re.fullmatch(r"qepA[A-Za-z0-9_-]*", n, flags=re.I):
        return "Fluoroquinolone", "qepA"

    if re.search(r"AAC\(6'\)-Ib-cr", n, flags=re.I):
        return "Fluoroquinolone", "aac6_Ib_cr"

    if re.fullmatch(r"KPC-\d+", n, flags=re.I):
        return "Carbapenem", "KPC"

    if re.fullmatch(r"NDM-\d+", n, flags=re.I):
        return "Carbapenem", "NDM"

    if re.fullmatch(r"VIM-\d+", n, flags=re.I):
        return "Carbapenem", "VIM"

    if re.fullmatch(r"IMP-\d+", n, flags=re.I):
        return "Carbapenem", "IMP"

    if n.upper() in oxa48_upper:
        return "Carbapenem", "OXA48_like"

    return None, None


classified = df["Best_Hit_ARO"].apply(classify_aro)

df["endpoint"] = [x[0] for x in classified]
df["panel_feature"] = [x[1] for x in classified]

hits = df[df["panel_feature"].notna()].copy()

agg = (
    hits.groupby(
        ["sample", "endpoint", "panel_feature"],
        as_index=False
    )
    .agg(
        ARG_loci=("ARG_loci", "sum"),
        mobile_ARG_loci=("mobile_ARG_loci", "sum"),
        plasmid_ARG_loci=("plasmid_ARG_loci", "sum"),
        detected_alleles=(
            "Best_Hit_ARO",
            lambda x: ";".join(sorted(set(map(str, x))))
        ),
    )
)

agg["present"] = (agg["ARG_loci"] > 0).astype(int)

agg["mobile_fraction"] = np.where(
    agg["ARG_loci"] > 0,
    agg["mobile_ARG_loci"] / agg["ARG_loci"],
    0.0,
)

agg["plasmid_fraction"] = np.where(
    agg["ARG_loci"] > 0,
    agg["plasmid_ARG_loci"] / agg["ARG_loci"],
    0.0,
)

existing_features = [
    ("Fluoroquinolone", "qnr"),
    ("Fluoroquinolone", "oqxA"),
    ("Fluoroquinolone", "oqxB"),
    ("Fluoroquinolone", "qepA"),
    ("Fluoroquinolone", "aac6_Ib_cr"),
    ("Carbapenem", "KPC"),
    ("Carbapenem", "NDM"),
    ("Carbapenem", "VIM"),
    ("Carbapenem", "IMP"),
    ("Carbapenem", "OXA48_like"),
]

grid = pd.MultiIndex.from_product(
    [samples, range(len(existing_features))],
    names=["sample", "feature_index"],
).to_frame(index=False)

grid["endpoint"] = grid["feature_index"].map(
    lambda i: existing_features[i][0]
)

grid["panel_feature"] = grid["feature_index"].map(
    lambda i: existing_features[i][1]
)

grid = grid.drop(columns="feature_index")

gene_long = grid.merge(
    agg,
    on=["sample", "endpoint", "panel_feature"],
    how="left",
)

numeric_cols = [
    "ARG_loci",
    "mobile_ARG_loci",
    "plasmid_ARG_loci",
    "present",
    "mobile_fraction",
    "plasmid_fraction",
]

for c in numeric_cols:
    gene_long[c] = gene_long[c].fillna(0)

gene_long["detected_alleles"] = (
    gene_long["detected_alleles"].fillna("")
)

gene_long.to_csv(GENE_LONG, sep="\t", index=False)

parts = []

for metric in [
    "ARG_loci",
    "present",
    "mobile_ARG_loci",
    "plasmid_ARG_loci",
    "mobile_fraction",
    "plasmid_fraction",
]:
    tmp = gene_long.pivot(
        index="sample",
        columns="panel_feature",
        values=metric,
    )

    tmp.columns = [
        f"{c}__{metric}" for c in tmp.columns
    ]

    parts.append(tmp)

gene_wide = pd.concat(parts, axis=1).reset_index()

fq_features = [
    "qnr",
    "oqxA",
    "oqxB",
    "qepA",
    "aac6_Ib_cr",
]

carb_features = [
    "KPC",
    "NDM",
    "VIM",
    "IMP",
    "OXA48_like",
]

gene_wide["FQ_PMQR_total_ARG_loci"] = sum(
    gene_wide[f"{f}__ARG_loci"] for f in fq_features
)

gene_wide["FQ_PMQR_feature_count"] = sum(
    gene_wide[f"{f}__present"] for f in fq_features
)

gene_wide["Carbapenemase_total_ARG_loci"] = sum(
    gene_wide[f"{f}__ARG_loci"] for f in carb_features
)

gene_wide["Carbapenemase_feature_count"] = sum(
    gene_wide[f"{f}__present"] for f in carb_features
)

gene_wide["FQ_QRDR_status"] = (
    "PENDING_gyrA_parC_mutation_calling"
)

gene_wide.to_csv(GENE_WIDE, sep="\t", index=False)


# ==============================================================
# 4. DISCOVER EXISTING RGI-MAIN FILES
# ==============================================================

candidate_files = []

for pattern in [
    "*_rgi_main_compact.tsv",
    "*_rgi_main.txt",
    "*_rgi_main_contigs.txt",
]:
    candidate_files.extend(ROOT.rglob(pattern))

candidate_files = sorted(set(candidate_files))

sample_files = {s: [] for s in samples}

for f in candidate_files:
    name = f.name
    pathstr = str(f)

    for s in samples:
        if s in name or f"/{s}/" in pathstr:
            sample_files[s].append(f)


# ==============================================================
# 5. EXTRACT STRICT/PERFECT E. COLI gyrA/parC MUTATIONS
# ==============================================================

wanted_names = {
    "Escherichia coli gyrA conferring resistance to fluoroquinolones":
        "gyrA",

    "Escherichia coli parC conferring resistance to fluoroquinolones":
        "parC",
}

records = []
audit = []

for sample in samples:

    files = sample_files[sample]

    readable = 0
    before = len(records)

    for f in files:

        try:
            with open(
                f,
                "r",
                encoding="utf-8",
                errors="replace",
                newline="",
            ) as handle:

                reader = csv.reader(handle, delimiter="\t")

                try:
                    header = next(reader)
                except StopIteration:
                    continue

                header = [str(x).strip() for x in header]

                required_rgi = {
                    "Best_Hit_ARO",
                    "Model_type",
                    "SNPs_in_Best_Hit_ARO",
                    "Cut_Off",
                }

                if not required_rgi.issubset(set(header)):
                    continue

                readable += 1
                idx = {c: i for i, c in enumerate(header)}

                for row in reader:

                    if len(row) < len(header):
                        continue

                    best = row[idx["Best_Hit_ARO"]].strip()

                    if best not in wanted_names:
                        continue

                    cutoff = row[idx["Cut_Off"]].strip()

                    if cutoff.lower() not in {"strict", "perfect"}:
                        continue

                    model = row[idx["Model_type"]].strip()

                    if "protein variant" not in model.lower():
                        continue

                    mutation = (
                        row[idx["SNPs_in_Best_Hit_ARO"]]
                        .strip()
                    )

                    if not mutation or mutation.lower() in {
                        "n/a", "na", "none", "nan"
                    }:
                        continue

                    rec = {
                        "sample": sample,
                        "gene": wanted_names[best],
                        "mutation": mutation,
                        "cutoff": cutoff,
                        "Best_Hit_ARO": best,
                        "Model_type": model,
                        "source_file": str(f),
                    }

                    for optional in [
                        "ORF_ID",
                        "Contig",
                        "Best_Identities",
                        "ARO",
                    ]:
                        if optional in idx:
                            rec[optional] = (
                                row[idx[optional]].strip()
                            )

                    records.append(rec)

        except Exception as exc:
            print("WARNING:", sample, f, exc)

    audit.append({
        "sample": sample,
        "n_candidate_files": len(files),
        "n_readable_files": readable,
        "n_strict_perfect_QRDR_calls":
            len(records) - before,
        "status":
            "OK" if readable > 0
            else "NO_READABLE_RGI_MAIN_TABLE",
    })

qrdr_hits = pd.DataFrame(records)

if len(qrdr_hits):

    dedup_cols = [
        c for c in [
            "sample",
            "gene",
            "mutation",
            "ORF_ID",
            "Contig",
            "ARO",
        ]
        if c in qrdr_hits.columns
    ]

    qrdr_hits = qrdr_hits.sort_values(
        ["sample", "gene", "mutation", "source_file"]
    )

    qrdr_hits = qrdr_hits.drop_duplicates(
        subset=dedup_cols,
        keep="first",
    )

else:
    qrdr_hits = pd.DataFrame(columns=[
        "sample",
        "gene",
        "mutation",
        "cutoff",
        "Best_Hit_ARO",
        "Model_type",
        "source_file",
    ])

qrdr_hits.to_csv(QRDR_LONG, sep="\t", index=False)

audit_df = pd.DataFrame(audit)
audit_df.to_csv(QRDR_AUDIT, sep="\t", index=False)


# ==============================================================
# 6. BUILD 101-SAMPLE QRDR MATRIX
# ==============================================================

qrdr_wide = pd.DataFrame({"sample": samples})

for gene in ["gyrA", "parC"]:

    sub = qrdr_hits[qrdr_hits["gene"] == gene].copy()

    positives = set(sub["sample"].astype(str))

    mutation_strings = (
        sub.groupby("sample")["mutation"]
        .apply(lambda x: ";".join(sorted(set(map(str, x)))))
        .to_dict()
    )

    qrdr_wide[f"Ecoli_{gene}_QRDR_present"] = (
        qrdr_wide["sample"].isin(positives).astype(int)
    )

    qrdr_wide[f"Ecoli_{gene}_QRDR_mutations"] = (
        qrdr_wide["sample"]
        .map(mutation_strings)
        .fillna("")
    )

qrdr_wide["Ecoli_QRDR_any_present"] = (
    (
        qrdr_wide["Ecoli_gyrA_QRDR_present"]
        + qrdr_wide["Ecoli_parC_QRDR_present"]
    ) > 0
).astype(int)

qrdr_wide["Ecoli_QRDR_gene_count"] = (
    qrdr_wide["Ecoli_gyrA_QRDR_present"]
    + qrdr_wide["Ecoli_parC_QRDR_present"]
)

if len(qrdr_hits):

    for gene in ["gyrA", "parC"]:

        sub = qrdr_hits[qrdr_hits["gene"] == gene]

        for mutation in sorted(
            sub["mutation"].dropna().unique()
        ):

            safe = re.sub(
                r"[^A-Za-z0-9]+",
                "_",
                str(mutation),
            ).strip("_")

            positives = set(
                sub.loc[
                    sub["mutation"] == mutation,
                    "sample",
                ].astype(str)
            )

            qrdr_wide[f"Ecoli_{gene}_{safe}"] = (
                qrdr_wide["sample"]
                .isin(positives)
                .astype(int)
            )

qrdr_wide.to_csv(QRDR_WIDE, sep="\t", index=False)


# ==============================================================
# 7. FINAL ENDPOINT-MATCHED MATRIX
# ==============================================================

genes = pd.read_csv(GENE_WIDE, sep="\t", low_memory=False)
qrdr = pd.read_csv(QRDR_WIDE, sep="\t", low_memory=False)

if "FQ_QRDR_status" in genes.columns:
    genes = genes.drop(columns=["FQ_QRDR_status"])

final = genes.merge(
    qrdr,
    on="sample",
    how="left",
    validate="one_to_one",
)

final["FQ_PMQR_any_present"] = (
    final[[f"{f}__present" for f in fq_features]]
    .max(axis=1)
    .astype(int)
)

final["FQ_PMQR_mechanism_count"] = (
    final[[f"{f}__present" for f in fq_features]]
    .sum(axis=1)
    .astype(int)
)

final["FQ_QRDR_any_present"] = (
    final["Ecoli_QRDR_any_present"]
    .fillna(0)
    .astype(int)
)

final["FQ_any_matched_resistance_feature"] = (
    (
        final["FQ_PMQR_any_present"]
        + final["FQ_QRDR_any_present"]
    ) > 0
).astype(int)

final["FQ_matched_feature_count"] = (
    final["FQ_PMQR_mechanism_count"]
    + final["Ecoli_gyrA_QRDR_present"].fillna(0).astype(int)
    + final["Ecoli_parC_QRDR_present"].fillna(0).astype(int)
)

final["Carbapenemase_any_present"] = (
    final[[f"{f}__present" for f in carb_features]]
    .max(axis=1)
    .astype(int)
)

final["Carbapenemase_family_count"] = (
    final[[f"{f}__present" for f in carb_features]]
    .sum(axis=1)
    .astype(int)
)

final["FQ_QRDR_source"] = (
    "CARD_RGI_MAIN_StrictPerfect_protein_variant_model"
)

final["FQ_QRDR_measurement_type"] = (
    "assembled_contig_mutation_presence"
)

final["FQ_QRDR_allele_fraction_available"] = 0

final = final.sort_values("sample").reset_index(drop=True)

final.to_csv(FINAL_OUT, sep="\t", index=False)


# ==============================================================
# 8. FEATURE DICTIONARY
# ==============================================================

dictionary = [
    ["qnr__present", "Fluoroquinolone", "PMQR",
     "Qnr family detected by CARD/RGI", "binary"],

    ["oqxA__present", "Fluoroquinolone", "PMQR_efflux",
     "oqxA detected by CARD/RGI", "binary"],

    ["oqxB__present", "Fluoroquinolone", "PMQR_efflux",
     "oqxB detected by CARD/RGI", "binary"],

    ["qepA__present", "Fluoroquinolone", "PMQR_efflux",
     "qepA detected by CARD/RGI", "binary"],

    ["aac6_Ib_cr__present", "Fluoroquinolone",
     "PMQR_drug_modification",
     "AAC(6')-Ib-cr family detected by CARD/RGI",
     "binary"],

    ["Ecoli_gyrA_QRDR_present", "Fluoroquinolone",
     "QRDR_target_mutation",
     "Strict/Perfect CARD protein-variant call in E. coli gyrA",
     "binary"],

    ["Ecoli_parC_QRDR_present", "Fluoroquinolone",
     "QRDR_target_mutation",
     "Strict/Perfect CARD protein-variant call in E. coli parC",
     "binary"],

    ["Ecoli_gyrA_QRDR_mutations", "Fluoroquinolone",
     "QRDR_target_mutation",
     "Exact CARD-reported gyrA amino-acid substitution(s)",
     "categorical_string"],

    ["Ecoli_parC_QRDR_mutations", "Fluoroquinolone",
     "QRDR_target_mutation",
     "Exact CARD-reported parC amino-acid substitution(s)",
     "categorical_string"],

    ["KPC__present", "Carbapenem", "carbapenemase",
     "KPC family detected", "binary"],

    ["NDM__present", "Carbapenem", "carbapenemase",
     "NDM family detected", "binary"],

    ["VIM__present", "Carbapenem", "carbapenemase",
     "VIM family detected", "binary"],

    ["IMP__present", "Carbapenem", "carbapenemase",
     "IMP family detected", "binary"],

    ["OXA48_like__present", "Carbapenem", "carbapenemase",
     "OXA-48-like family detected", "binary"],
]

pd.DataFrame(
    dictionary,
    columns=[
        "feature",
        "endpoint",
        "mechanism",
        "definition",
        "data_type",
    ],
).to_csv(DICT_OUT, sep="\t", index=False)


# ==============================================================
# 9. FINAL QC
# ==============================================================

print("==============================================================")
print("ENDPOINT-MATCHED AMR FEATURE BUILD")
print("==============================================================")
print("Samples:", len(final))
print("Unique samples:", final["sample"].nunique())
print("Duplicate samples:", final["sample"].duplicated().sum())
print()

print("RGI-MAIN readable:",
      int((audit_df["n_readable_files"] > 0).sum()),
      "/ 101")

print()
print("FLUOROQUINOLONE")
print("qnr:", int(final["qnr__present"].sum()))
print("oqxA:", int(final["oqxA__present"].sum()))
print("oqxB:", int(final["oqxB__present"].sum()))
print("qepA:", int(final["qepA__present"].sum()))
print("aac6-Ib-cr:", int(final["aac6_Ib_cr__present"].sum()))
print("PMQR any:", int(final["FQ_PMQR_any_present"].sum()))
print("QRDR any:", int(final["FQ_QRDR_any_present"].sum()))
print("gyrA:", int(final["Ecoli_gyrA_QRDR_present"].sum()))
print("parC:", int(final["Ecoli_parC_QRDR_present"].sum()))
print("Any matched FQ:",
      int(final["FQ_any_matched_resistance_feature"].sum()))

print()
print("Observed QRDR mutations:")
if len(qrdr_hits):
    print(
        qrdr_hits.groupby(["gene", "mutation"])["sample"]
        .nunique()
        .reset_index(name="samples")
        .sort_values(
            ["gene", "samples", "mutation"],
            ascending=[True, False, True],
        )
        .to_string(index=False)
    )

print()
print("CARBAPENEM")
print("Any carbapenemase:",
      int(final["Carbapenemase_any_present"].sum()))

for f in carb_features:
    print(f"{f:12s}:",
          int(final[f"{f}__present"].sum()))

final_pass = (
    len(final) == 101
    and final["sample"].nunique() == 101
    and final["sample"].duplicated().sum() == 0
    and int((audit_df["n_readable_files"] > 0).sum()) == 101
)

print()
print("FINAL PASS:", final_pass)

if not final_pass:
    raise RuntimeError("FINAL QC FAILED")
