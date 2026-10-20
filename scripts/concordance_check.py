import pandas as pd
geno = pd.read_csv("tbprofiler.txt", sep="\t")
geno = geno.rename(columns={geno.columns[0]: "isolate_id"})
pheno = pd.read_csv("data/filtered_metadata.csv", dtype={"isolate_id": str})
drug_map = {
    "rifampicin": "RIF", "isoniazid": "INH", "pyrazinamide": "PZA",
    "ethambutol": "EMB", "streptomycin": "STR", "moxifloxacin": "MXF",
    "levofloxacin": "LEV", "ofloxacin": "OFX", "kanamycin": "KAN",
    "amikacin": "AMK", "capreomycin": "CAP", "ethionamide": "ETH",
}
merged = geno.merge(pheno, on="isolate_id", how="inner")
print(f"Matched isolates with both genotype and phenotype: {len(merged)}")
rows = []
for gcol, pcol in drug_map.items():
    if gcol not in merged.columns or pcol not in merged.columns:
        continue
    sub = merged[[gcol, pcol]].dropna()
    sub = sub[sub[pcol].isin(["R", "S"])]
    if sub.empty:
        continue
    geno_call = (sub[gcol] != "-").map({True: "R", False: "S"})
    pheno_call = sub[pcol]
    tp = ((geno_call == "R") & (pheno_call == "R")).sum()
    tn = ((geno_call == "S") & (pheno_call == "S")).sum()
    fp = ((geno_call == "R") & (pheno_call == "S")).sum()
    fn = ((geno_call == "S") & (pheno_call == "R")).sum()
    n = tp + tn + fp + fn
    rows.append({
        "drug": pcol, "n": n, "TP": tp, "TN": tn, "FP": fp, "FN": fn,
        "sensitivity_%": round(100*tp/(tp+fn), 2) if (tp+fn) else None,
        "specificity_%": round(100*tn/(tn+fp), 2) if (tn+fp) else None,
        "PPV_%": round(100*tp/(tp+fp), 2) if (tp+fp) else None,
        "NPV_%": round(100*tn/(tn+fn), 2) if (tn+fn) else None,
        "accuracy_%": round(100*(tp+tn)/n, 2) if n else None,
    })
out = pd.DataFrame(rows)
print(out.to_string(index=False))
out.to_csv("genotype_phenotype_concordance.csv", index=False)
