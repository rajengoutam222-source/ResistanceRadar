#!/usr/bin/env python3
#@rgoutam
from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("rr-analysis/hendriksen_sensitivity")
TRUTH_XLSX = BASE / "truth_set" / "41467_2019_8853_MOESM7_ESM.xlsx"
MAP_FILE = BASE / "metadata" / "Hendriksen_candidate_8_ENA_mapping_curated.tsv"
RGI_DRUG_FILE = BASE / "tables" / "Stage3_4_RGI_sample_by_drug_class.tsv"

TABLE_DIR = BASE / "tables"
FIG_DIR = BASE / "figures"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

def norm(x):
    x = str(x).strip().lower()
    x = re.sub(r"[^a-z0-9]+", "_", x)
    return re.sub(r"_+", "_", x).strip("_")

def truth_gene_to_class(gene):
    g = str(gene).lower()
    if any(k in g for k in ["bla", "oxa", "ctx", "tem", "shv", "ampc", "cmy", "dha", "kpc", "ndm", "vim", "imp", "ges", "mec"]):
        return "beta_lactam"
    if any(k in g for k in ["aac", "aad", "aph", "ant", "arm", "rmt", "str"]):
        return "aminoglycoside"
    if "tet" in g:
        return "tetracycline"
    if any(k in g for k in ["erm", "mph", "msr", "lnu", "lsa", "vga", "vat", "vgb", "mef"]):
        return "macrolide_lincosamide_streptogramin"
    if "sul" in g:
        return "sulfonamide"
    if "dfr" in g:
        return "trimethoprim"
    if any(k in g for k in ["qnr", "oqx", "qep", "aac(6')-ib-cr"]):
        return "quinolone"
    if any(k in g for k in ["cat", "cml", "flo", "fex"]):
        return "phenicol"
    if "van" in g:
        return "glycopeptide"
    if "mcr" in g:
        return "colistin"
    if "arr" in g:
        return "rifamycin"
    if "fos" in g:
        return "fosfomycin"
    return "other_or_unmapped"

def card_class_to_broad(x):
    c = norm(x)
    if any(k in c for k in ["beta_lactam", "cephalosporin", "carbapenem", "penam", "monobactam"]):
        return "beta_lactam"
    if "aminoglycoside" in c:
        return "aminoglycoside"
    if "tetracycline" in c:
        return "tetracycline"
    if any(k in c for k in ["macrolide", "lincosamide", "streptogramin"]):
        return "macrolide_lincosamide_streptogramin"
    if "sulfonamide" in c:
        return "sulfonamide"
    if "trimethoprim" in c:
        return "trimethoprim"
    if "quinolone" in c or "fluoroquinolone" in c:
        return "quinolone"
    if "phenicol" in c or "chloramphenicol" in c:
        return "phenicol"
    if "glycopeptide" in c:
        return "glycopeptide"
    if "polymyxin" in c or "colistin" in c:
        return "colistin"
    if "rifamycin" in c or "rifampin" in c:
        return "rifamycin"
    if "fosfomycin" in c:
        return "fosfomycin"
    return "other_or_unmapped"

# sample mapping
mapping = pd.read_csv(MAP_FILE, sep="\t")
hend_col = [c for c in mapping.columns if "hendriksen" in c.lower() or "sample" in c.lower()][0]
run_col = [c for c in mapping.columns if "err" in c.lower() or "ena" in c.lower() or "run" in c.lower()][0]
mapping = mapping.rename(columns={hend_col: "hendriksen_sample", run_col: "sample_id"})
mapping["hendriksen_sample"] = mapping["hendriksen_sample"].astype(str).str.strip()
mapping["sample_id"] = mapping["sample_id"].astype(str).str.strip()

# truth set
xl = pd.ExcelFile(TRUTH_XLSX)
sheet = "ResFind.Gene.count"
truth = pd.read_excel(TRUTH_XLSX, sheet_name=sheet)
truth = truth.rename(columns={truth.columns[0]: "truth_gene"})
truth["truth_gene"] = truth["truth_gene"].astype(str).str.strip()

truth_rows = []
for _, m in mapping.iterrows():
    hend = m["hendriksen_sample"]
    run = m["sample_id"]
    if hend not in truth.columns:
        print(f"WARNING: {hend} not found in truth columns")
        continue
    sub = truth[["truth_gene", hend]].copy()
    sub[hend] = pd.to_numeric(sub[hend], errors="coerce").fillna(0)
    sub = sub[sub[hend] > 0].copy()
    sub["sample_id"] = run
    sub["hendriksen_sample"] = hend
    sub["truth_class"] = sub["truth_gene"].apply(truth_gene_to_class)
    truth_rows.append(sub)

truth_long = pd.concat(truth_rows, ignore_index=True)
truth_long.to_csv(TABLE_DIR / "Stage3_4B_truth_genes_long.tsv", sep="\t", index=False)

truth_classes = (
    truth_long[truth_long["truth_class"] != "other_or_unmapped"]
    .groupby("sample_id")["truth_class"]
    .apply(lambda x: set(x))
    .to_dict()
)

# RGI classes
rgi = pd.read_csv(RGI_DRUG_FILE, sep="\t")
drug_col = [c for c in rgi.columns if "drug" in c.lower() and "class" in c.lower()][0]
rgi["rgi_class"] = rgi[drug_col].apply(card_class_to_broad)

rgi_classes = (
    rgi[rgi["rgi_class"] != "other_or_unmapped"]
    .groupby("sample_id")["rgi_class"]
    .apply(lambda x: set(x))
    .to_dict()
)

# sensitivity
rows = []
for _, m in mapping.iterrows():
    sample = m["sample_id"]
    hend = m["hendriksen_sample"]
    tset = truth_classes.get(sample, set())
    rset = rgi_classes.get(sample, set())
    recovered = tset & rset
    missed = tset - rset
    extra = rset - tset
    sens = len(recovered) / len(tset) if len(tset) else None

    rows.append({
        "sample_id": sample,
        "hendriksen_sample": hend,
        "n_truth_classes": len(tset),
        "n_recovered_truth_classes": len(recovered),
        "sensitivity": sens,
        "sensitivity_percent": sens * 100 if sens is not None else None,
        "truth_classes": ";".join(sorted(tset)),
        "recovered_truth_classes": ";".join(sorted(recovered)),
        "missed_truth_classes": ";".join(sorted(missed)),
        "extra_rgi_classes_not_in_truth": ";".join(sorted(extra)),
    })

sens_df = pd.DataFrame(rows)
sens_df.to_csv(TABLE_DIR / "Stage3_4B_Hendriksen_class_level_sensitivity_by_sample.tsv", sep="\t", index=False)

overall_truth = set()
overall_recovered = set()
for _, row in sens_df.iterrows():
    overall_truth |= set([x for x in row["truth_classes"].split(";") if x])
    overall_recovered |= set([x for x in row["recovered_truth_classes"].split(";") if x])

overall = pd.DataFrame([{
    "n_samples": len(sens_df),
    "overall_truth_class_union": len(overall_truth),
    "overall_recovered_truth_class_union": len(overall_recovered),
    "overall_sensitivity": len(overall_recovered) / len(overall_truth) if overall_truth else None,
    "overall_sensitivity_percent": 100 * len(overall_recovered) / len(overall_truth) if overall_truth else None,
    "truth_classes_union": ";".join(sorted(overall_truth)),
    "recovered_classes_union": ";".join(sorted(overall_recovered)),
}])
overall.to_csv(TABLE_DIR / "Stage3_4B_Hendriksen_class_level_sensitivity_overall.tsv", sep="\t", index=False)

# plot 1
plot_df = sens_df.sort_values("sensitivity_percent")
plt.figure(figsize=(9, 5))
plt.barh(plot_df["sample_id"], plot_df["sensitivity_percent"])
plt.xlabel("Class-level sensitivity (%)")
plt.ylabel("ENA sample")
plt.title("Stage 3.4B Hendriksen validation: CARD/RGI class-level sensitivity")
plt.xlim(0, 105)
for i, (_, row) in enumerate(plot_df.iterrows()):
    plt.text(min(row["sensitivity_percent"] + 1, 101), i, f"{row['n_recovered_truth_classes']}/{row['n_truth_classes']}", va="center", fontsize=8)
plt.tight_layout()
plt.savefig(FIG_DIR / "Stage3_4B_Hendriksen_class_level_sensitivity_by_sample.png", dpi=300)
plt.savefig(FIG_DIR / "Stage3_4B_Hendriksen_class_level_sensitivity_by_sample.pdf")
plt.close()

# plot 2
plot_df = sens_df.sort_values("sample_id").copy()
plot_df["n_missed_truth_classes"] = plot_df["n_truth_classes"] - plot_df["n_recovered_truth_classes"]
x = range(len(plot_df))
plt.figure(figsize=(9, 5))
plt.bar(x, plot_df["n_recovered_truth_classes"], label="Recovered")
plt.bar(x, plot_df["n_missed_truth_classes"], bottom=plot_df["n_recovered_truth_classes"], label="Missed")
plt.xticks(x, plot_df["sample_id"], rotation=45, ha="right")
plt.ylabel("Truth drug classes")
plt.xlabel("ENA sample")
plt.title("Stage 3.4B Hendriksen validation: recovered vs missed truth classes")
plt.legend()
plt.tight_layout()
plt.savefig(FIG_DIR / "Stage3_4B_Hendriksen_recovered_vs_missed_classes.png", dpi=300)
plt.savefig(FIG_DIR / "Stage3_4B_Hendriksen_recovered_vs_missed_classes.pdf")
plt.close()

print("\nOverall:")
print(overall.to_string(index=False))

print("\nPer sample:")
print(sens_df[["sample_id","hendriksen_sample","n_recovered_truth_classes","n_truth_classes","sensitivity_percent"]].to_string(index=False))

print("\nWrote:")
print(TABLE_DIR / "Stage3_4B_Hendriksen_class_level_sensitivity_by_sample.tsv")
print(TABLE_DIR / "Stage3_4B_Hendriksen_class_level_sensitivity_overall.tsv")
print(FIG_DIR / "Stage3_4B_Hendriksen_class_level_sensitivity_by_sample.png")
print(FIG_DIR / "Stage3_4B_Hendriksen_recovered_vs_missed_classes.png")
