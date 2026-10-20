import pandas as pd
pheno = pd.read_csv("data/filtered_metadata.csv", dtype=str)

print("=== RIF.1 value counts ===")
print(pheno["RIF.1"].value_counts(dropna=False))

print("\n=== RIF vs RIF.1 agreement ===")
both = pheno[["RIF","RIF.1"]].dropna()
agree = (both["RIF"] == both["RIF.1"]).sum()
print(f"Agree: {agree}/{len(both)} ({100*agree/len(both):.1f}%)")
print("\nSample disagreements:")
print(both[both["RIF"] != both["RIF.1"]].head(10))
