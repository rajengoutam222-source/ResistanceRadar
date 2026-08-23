
from pathlib import Path
import pandas as pd
import numpy as np
import json, hashlib, re
from datetime import datetime

ROOT = Path("/rsstu/users/s/sleblan/MismatchRepair/Project/pipeline/rr-pipeline")
MANIFEST = ROOT / "rr-analysis/ResistanceRadar_projectwide_ARG_mobility_manifest_CORRECTED.tsv"
OUT = ROOT / "rr-analysis/projectwide_ml_features"
OUT.mkdir(parents=True, exist_ok=True)

CLASSES = [
    "Aminoglycoside","Beta-lactam","Trimethoprim","Quinolone","Tetracycline",
    "MLS","Phenicol","Fosfomycin","Rifamycin","Sulfonamide"
]

CARD_MAP = {
    "aminoglycoside antibiotic":"Aminoglycoside",
    "carbapenem":"Beta-lactam","cephalosporin":"Beta-lactam",
    "cephamycin":"Beta-lactam","monobactam":"Beta-lactam",
    "penam":"Beta-lactam","penem":"Beta-lactam",
    "diaminopyrimidine antibiotic":"Trimethoprim",
    "fluoroquinolone antibiotic":"Quinolone",
    "tetracycline antibiotic":"Tetracycline","glycylcycline":"Tetracycline",
    "macrolide antibiotic":"MLS","lincosamide antibiotic":"MLS",
    "streptogramin antibiotic":"MLS","streptogramin a antibiotic":"MLS",
    "streptogramin b antibiotic":"MLS",
    "phenicol antibiotic":"Phenicol",
    "phosphonic acid antibiotic":"Fosfomycin",
    "rifamycin antibiotic":"Rifamycin",
    "sulfonamide antibiotic":"Sulfonamide",
}

def shannon(vals):
    a = np.asarray([x for x in vals if x > 0], dtype=float)
    if len(a) == 0:
        return 0.0
    p = a / a.sum()
    return float(-(p * np.log(p)).sum())

def bools(s):
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y","t"})

def clean_name(x):
    return re.sub(r"[^A-Za-z0-9]+","_",str(x)).strip("_")

manifest = pd.read_csv(MANIFEST, sep="\t")
all_samples = sorted(manifest["sample"].astype(str).unique())
N = len(all_samples)

coverage = (
    manifest[["mode","sample"]].drop_duplicates()
    .groupby("mode")["sample"].nunique()
    .to_dict()
)

strict_n = int(coverage.get("strict_perfect",0))
primary_n = int(coverage.get("primary",0))

if strict_n == N:
    MODE = "strict_perfect"
    COMPLETE = True
elif primary_n == N:
    MODE = "primary"
    COMPLETE = True
elif strict_n >= primary_n:
    MODE = "strict_perfect"
    COMPLETE = False
else:
    MODE = "primary"
    COMPLETE = False

print("===== MODE COVERAGE =====")
print("Global unique samples:", N)
print("Strict/perfect:", strict_n)
print("Primary:", primary_n)
print("Selected mode:", MODE)
print("Selected mode covers all:", COMPLETE)

# Save coverage before any heavy work
pd.DataFrame([
    {"global_unique_samples":N,"strict_perfect_samples":strict_n,
     "primary_samples":primary_n,"selected_mode":MODE,
     "selected_mode_complete":COMPLETE}
]).to_csv(OUT/"ResistanceRadar_projectwide_mode_coverage.tsv", sep="\t", index=False)

# Pick one source per sample, preferring sample-specific dataset files
m = manifest[manifest["mode"].eq(MODE)].copy()
memberships = (
    manifest.groupby("sample")["cohort"]
    .agg(lambda s: ";".join(sorted(set(map(str,s)))))
    .to_dict()
)

def rank_source(r):
    p = Path(r["source_file"])
    sample = str(r["sample"])
    score = 0
    if p.parent.name == sample: score += 100
    if sample in p.name: score += 50
    if "dataset_inputs" in str(p): score += 30
    if "kathmandu_kbase" in str(p) or "annual_wastewater_2017_2025" in str(p): score += 20
    if "hendriksen_sensitivity" in str(p): score -= 20
    return score

m["_rank"] = m.apply(rank_source, axis=1)
m = m.sort_values(["sample","_rank","source_file"], ascending=[True,False,True])
m = m.drop_duplicates("sample", keep="first")
m["cohort_memberships"] = m["sample"].map(memberships)

m[["sample","cohort_memberships","mode","source_file"]].to_csv(
    OUT/"ResistanceRadar_projectwide_source_selection.tsv", sep="\t", index=False
)

missing_mode = sorted(set(all_samples)-set(m["sample"]))
with open(OUT/"ResistanceRadar_projectwide_missing_selected_mode_samples.txt","w") as fh:
    for s in missing_mode:
        fh.write(s+"\n")

sample_rows = []
class_rows = []
mechanism_rows = []
aro_rows = []
problem_rows = []

for idx, r in m.reset_index(drop=True).iterrows():
    sample = str(r["sample"])
    f = Path(r["source_file"])
    try:
        d = pd.read_csv(f, sep="\t", low_memory=False)

        scol = next((c for c in ["sample_id","sample","Sample","SAMPLE","Sample_ID"] if c in d.columns), None)
        if scol is not None:
            d = d[d[scol].astype(str).str.strip().eq(sample)].copy()
        if d.empty:
            raise RuntimeError("no rows after sample filtering")

        for c in ["contig_id","ORF_ID","Drug_Class","Resistance_Mechanism",
                  "Best_Hit_ARO","ARO","matched_to_genomad","genomad_context"]:
            if c not in d.columns:
                d[c] = pd.NA

        if d["ORF_ID"].notna().any():
            d["locus_id"] = d["contig_id"].fillna("").astype(str)+"|"+d["ORF_ID"].fillna("").astype(str)
        else:
            d["locus_id"] = (
                d["contig_id"].fillna("").astype(str)+"|"+
                d["Best_Hit_ARO"].fillna("").astype(str)+"|"+
                d["ARO"].fillna("").astype(str)+"|"+d.index.astype(str)
            )

        d["_mobile"] = bools(d["matched_to_genomad"])
        ctx = d["genomad_context"].fillna("").astype(str).str.lower()
        d["_plasmid"] = d["_mobile"] & ctx.str.contains("plasmid_like", regex=False)
        d["_virus"] = d["_mobile"] & ctx.str.contains("virus_like", regex=False)
        d["_chromosome"] = ctx.str.contains("chromosome_like", regex=False)

        # One record per ARG locus for global sample features
        locus = (
            d.groupby("locus_id",as_index=False)
            .agg(mobile=("_mobile","max"),plasmid=("_plasmid","max"),
                 virus=("_virus","max"),chromosome=("_chromosome","max"))
        )
        total = len(locus)
        nmobile = int(locus["mobile"].sum())
        nplasmid = int(locus["plasmid"].sum())
        nvirus = int(locus["virus"].sum())
        nchrom = int(locus["chromosome"].sum())

        # Expand CARD drug classes and map to locked 10-class ontology
        c = d[["locus_id","Drug_Class","_mobile","_plasmid","_virus","_chromosome"]].copy()
        c["raw_class"] = c["Drug_Class"].fillna("").astype(str).str.split(";")
        c = c.explode("raw_class")
        c["raw_class"] = c["raw_class"].astype(str).str.strip().str.lower()
        c["common_class"] = c["raw_class"].map(CARD_MAP)
        c = c.dropna(subset=["common_class"])
        c = (
            c.groupby(["locus_id","common_class"],as_index=False)
            .agg(mobile=("_mobile","max"),plasmid=("_plasmid","max"),
                 virus=("_virus","max"),chromosome=("_chromosome","max"))
        )

        counts = {}
        for cls in CLASSES:
            q = c[c["common_class"].eq(cls)]
            n = len(q)
            mob = int(q["mobile"].sum())
            pla = int(q["plasmid"].sum())
            vir = int(q["virus"].sum())
            chrn = int(q["chromosome"].sum())
            counts[cls] = n
            class_rows.append({
                "sample":sample,"common_class":cls,
                "ARG_loci":n,"ARG_present":int(n>0),
                "mobile_ARG_loci":mob,"mobile_ARG_fraction":mob/n if n else 0.0,
                "plasmid_ARG_loci":pla,"plasmid_ARG_fraction":pla/n if n else 0.0,
                "virus_ARG_loci":vir,"virus_ARG_fraction":vir/n if n else 0.0,
                "chromosome_ARG_loci":chrn,"chromosome_ARG_fraction":chrn/n if n else 0.0,
            })

        denom10 = sum(counts.values())
        for row in class_rows[-len(CLASSES):]:
            row["ARG_locus_fraction"] = row["ARG_loci"]/denom10 if denom10 else 0.0

        ten_loci = int(c["locus_id"].nunique())
        ten_mobile = int(c.loc[c["mobile"],"locus_id"].nunique())
        ten_plasmid = int(c.loc[c["plasmid"],"locus_id"].nunique())

        sample_rows.append({
            "sample":sample,
            "cohort_memberships":r["cohort_memberships"],
            "source_mode":MODE,
            "total_ARG_loci":total,
            "mobile_ARG_loci":nmobile,
            "mobile_ARG_fraction":nmobile/total if total else 0.0,
            "plasmid_like_ARG_loci":nplasmid,
            "plasmid_like_ARG_fraction":nplasmid/total if total else 0.0,
            "virus_like_ARG_loci":nvirus,
            "virus_like_ARG_fraction":nvirus/total if total else 0.0,
            "chromosome_like_ARG_loci":nchrom,
            "chromosome_like_ARG_fraction":nchrom/total if total else 0.0,
            "canonical10_unique_ARG_loci":ten_loci,
            "canonical10_mobile_ARG_loci":ten_mobile,
            "canonical10_mobile_ARG_fraction":ten_mobile/ten_loci if ten_loci else 0.0,
            "canonical10_plasmid_ARG_loci":ten_plasmid,
            "canonical10_plasmid_ARG_fraction":ten_plasmid/ten_loci if ten_loci else 0.0,
            "ARG_10class_Shannon":shannon(counts.values()),
            "ARG_10class_richness":sum(v>0 for v in counts.values()),
        })

        # Resistance-mechanism long table
        me = d[["locus_id","Resistance_Mechanism","_mobile","_plasmid"]].copy()
        me["mechanism"] = me["Resistance_Mechanism"].fillna("").astype(str).str.split(";")
        me = me.explode("mechanism")
        me["mechanism"] = me["mechanism"].astype(str).str.strip()
        me = me[me["mechanism"].ne("")]
        me = me.groupby(["locus_id","mechanism"],as_index=False).agg(
            mobile=("_mobile","max"),plasmid=("_plasmid","max"))
        for mech,q in me.groupby("mechanism"):
            n=len(q); mob=int(q["mobile"].sum()); pla=int(q["plasmid"].sum())
            mechanism_rows.append({
                "sample":sample,"Resistance_Mechanism":mech,"ARG_loci":n,
                "mobile_ARG_loci":mob,"mobile_ARG_fraction":mob/n if n else 0.0,
                "plasmid_ARG_loci":pla,"plasmid_ARG_fraction":pla/n if n else 0.0
            })

        # ARO-level long table
        ar = d[["locus_id","Best_Hit_ARO","ARO","_mobile","_plasmid"]].copy()
        ar["Best_Hit_ARO"] = ar["Best_Hit_ARO"].fillna("").astype(str).str.strip()
        ar["ARO"] = ar["ARO"].fillna("").astype(str).str.strip()
        ar = ar[(ar["Best_Hit_ARO"].ne("")) | (ar["ARO"].ne(""))]
        ar = ar.groupby(["locus_id","Best_Hit_ARO","ARO"],as_index=False).agg(
            mobile=("_mobile","max"),plasmid=("_plasmid","max"))
        for (name,aro),q in ar.groupby(["Best_Hit_ARO","ARO"],dropna=False):
            n=len(q); mob=int(q["mobile"].sum()); pla=int(q["plasmid"].sum())
            aro_rows.append({
                "sample":sample,"Best_Hit_ARO":name,"ARO":aro,"ARG_loci":n,
                "mobile_ARG_loci":mob,"mobile_ARG_fraction":mob/n if n else 0.0,
                "plasmid_ARG_loci":pla,"plasmid_ARG_fraction":pla/n if n else 0.0
            })

    except Exception as e:
        problem_rows.append({"sample":sample,"source_file":str(f),"problem":str(e)})

    if (idx+1)%10==0 or idx+1==len(m):
        print(f"Processed {idx+1}/{len(m)} source samples")

sf = pd.DataFrame(sample_rows).sort_values("sample")
cf = pd.DataFrame(class_rows).sort_values(["sample","common_class"])
mf = pd.DataFrame(mechanism_rows).sort_values(["sample","Resistance_Mechanism"]) if mechanism_rows else pd.DataFrame()
af = pd.DataFrame(aro_rows).sort_values(["sample","Best_Hit_ARO","ARO"]) if aro_rows else pd.DataFrame()
pf = pd.DataFrame(problem_rows)

sf.to_csv(OUT/"ResistanceRadar_projectwide_sample_features.tsv",sep="\t",index=False)
cf.to_csv(OUT/"ResistanceRadar_projectwide_10class_features_long.tsv",sep="\t",index=False)
mf.to_csv(OUT/"ResistanceRadar_projectwide_mechanism_features_long.tsv",sep="\t",index=False)
af.to_csv(OUT/"ResistanceRadar_projectwide_ARO_features_long.tsv",sep="\t",index=False)
pf.to_csv(OUT/"ResistanceRadar_projectwide_problems.tsv",sep="\t",index=False)

# Wide 10-class ML matrix
wide = sf.set_index("sample")
for metric in ["ARG_loci","ARG_present","ARG_locus_fraction",
               "mobile_ARG_loci","mobile_ARG_fraction",
               "plasmid_ARG_loci","plasmid_ARG_fraction"]:
    w = cf.pivot(index="sample",columns="common_class",values=metric).fillna(0)
    w.columns = [f"{metric}__{clean_name(c)}" for c in w.columns]
    wide = wide.join(w,how="left")
wide.reset_index().to_csv(
    OUT/"ResistanceRadar_projectwide_10class_features_wide.tsv",sep="\t",index=False
)

pd.DataFrame(sorted(CARD_MAP.items()),columns=["card_rgi_main_term","common_class"]).to_csv(
    OUT/"ResistanceRadar_CARD_to_10class_mapping.tsv",sep="\t",index=False
)

success = sf["sample"].nunique()
qc = pd.DataFrame([{
    "global_unique_samples":N,
    "selected_mode":MODE,
    "selected_mode_manifest_samples":m["sample"].nunique(),
    "successfully_processed_samples":success,
    "processing_failures":len(pf),
    "missing_selected_mode_samples":len(missing_mode),
    "analysis_ready_all_samples":bool(success==N and len(pf)==0 and len(missing_mode)==0)
}])
qc.to_csv(OUT/"ResistanceRadar_projectwide_QC.tsv",sep="\t",index=False)

status = {
    "created":datetime.now().isoformat(),
    "global_unique_samples":N,
    "strict_perfect_samples":strict_n,
    "primary_samples":primary_n,
    "selected_mode":MODE,
    "successfully_processed_samples":int(success),
    "analysis_ready_all_samples":bool(success==N and len(pf)==0 and len(missing_mode)==0),
    "important_note":"No primary/strict_perfect mixing. Features are RGI-MAIN locus based; genome-equivalent normalization is not applied project-wide here."
}
with open(OUT/"ResistanceRadar_projectwide_feature_status.json","w") as fh:
    json.dump(status,fh,indent=2)

# Checksums and sizes so later Git push can include only small validated artifacts
arts=[]
for f in sorted(OUT.glob("*")):
    if f.is_file() and f.name!="ResistanceRadar_projectwide_artifact_manifest.tsv":
        h=hashlib.sha256()
        with open(f,"rb") as fh:
            for chunk in iter(lambda:fh.read(1024*1024),b""):
                h.update(chunk)
        arts.append({"file":str(f.relative_to(ROOT)),"bytes":f.stat().st_size,"sha256":h.hexdigest()})
pd.DataFrame(arts).to_csv(
    OUT/"ResistanceRadar_projectwide_artifact_manifest.tsv",sep="\t",index=False
)

print("\n===== FINAL =====")
print(qc.to_string(index=False))
print("\nSaved:",OUT)
print("Artifacts:",len(arts)+1)
