#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path

META    = "data/filtered_metadata.csv"
VCF_DIR = "05_filtered_vcfs"
OUT_DIR = "data/prevalence_results"

DRUGS           = ['RIF','INH','PZA','EMB','STR','MXF','LEV','OFX','KAN','AMK','CAP','ETH']
FLUOROQUINOLONES = ['MXF','LEV','OFX']
INJECTABLES      = ['KAN','AMK','CAP']
COUNTRY_COL      = "isolation_country"

os.makedirs(OUT_DIR, exist_ok=True)
meta = pd.read_csv(META)

vcf_ids = set(
    p.name.replace('_filtered.vcf.gz', '').replace('_filtered.vcf', '')
    for p in Path(VCF_DIR).glob('*_filtered.vcf*')
    if not p.name.endswith('.idx')
)

df = meta[meta['isolate_id'].isin(vcf_ids)].copy()
print(f"Cohort: {len(df)} samples with completed VCFs (out of {len(meta)} total in metadata)\n")

rows = []
for d in DRUGS:
    r = (df[d] == 'R').sum()
    s = (df[d] == 'S').sum()
    n = r + s
    rows.append({
        "drug": d, "R": int(r), "S": int(s), "tested": int(n),
        "resistance_pct": round(100 * r / n, 1) if n else None,
        "coverage_pct": round(100 * n / len(df), 1) if len(df) else None,
    })
prevalence_df = pd.DataFrame(rows)
prevalence_df.to_csv(f"{OUT_DIR}/drug_prevalence.csv", index=False)
print("=== Per-drug resistance prevalence ===")
print(prevalence_df.to_string(index=False))

mdr_mask = (df['RIF'] == 'R') & (df['INH'] == 'R')
mdr_count = int(mdr_mask.sum())
if len(df):
    print(f"\nMDR-TB (RIF+INH resistant): {mdr_count} / {len(df)} ({100*mdr_count/len(df):.1f}%)")

fq_r  = df[FLUOROQUINOLONES].eq('R').any(axis=1)
inj_r = df[INJECTABLES].eq('R').any(axis=1)
xdr_mask = mdr_mask & fq_r & inj_r
xdr_count = int(xdr_mask.sum())
if len(df):
    print(f"XDR-TB (MDR + fluoroquinolone-R + injectable-R): {xdr_count} / {len(df)} ({100*xdr_count/len(df):.1f}%)")

print(f"\n=== Top 10 countries by sample count ===")
if COUNTRY_COL in df.columns:
    print(df[COUNTRY_COL].value_counts().head(10).to_string())

print(f"\n=== MDR-TB rate by country (min 20 samples) ===")
if COUNTRY_COL in df.columns:
    counts = df[COUNTRY_COL].value_counts()
    for c in counts[counts >= 20].index:
        sub = df[df[COUNTRY_COL] == c]
        sub_mdr = int(((sub['RIF'] == 'R') & (sub['INH'] == 'R')).sum())
        print(f"  {c:<20} n={len(sub):>5}  MDR={sub_mdr:>4} ({100*sub_mdr/len(sub):.1f}%)")

if 'strain_type' in df.columns:
    print(f"\nUnique strain types represented: {df['strain_type'].nunique()}")

df.to_csv(f"{OUT_DIR}/cohort_full.csv", index=False)
print(f"\nSaved: {OUT_DIR}/drug_prevalence.csv, {OUT_DIR}/cohort_full.csv")
