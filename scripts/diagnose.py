import pandas as pd
geno = pd.read_csv("tbprofiler.txt", sep="\t")
geno = geno.rename(columns={geno.columns[0]: "isolate_id"})
pheno = pd.read_csv("data/filtered_metadata.csv", dtype=str)

print("=== metadata columns ===")
print(list(pheno.columns))

print("\n=== RIF value counts ===")
print(pheno["RIF"].value_counts(dropna=False))

print("\n=== STR value counts ===")
print(pheno["STR"].value_counts(dropna=False))

print("\n=== INH value counts (control - this one worked) ===")
print(pheno["INH"].value_counts(dropna=False))
