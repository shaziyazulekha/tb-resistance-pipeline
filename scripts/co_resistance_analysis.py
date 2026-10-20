import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CSV  = "data/prevalence_results/cohort_full.csv"
PREV_CSV = "data/prevalence_results/drug_prevalence.csv"
OUT_DIR = "data/prevalence_results"
DRUGS = ['RIF','INH','PZA','EMB','STR','MXF','LEV','OFX','KAN','AMK','CAP','ETH']

df = pd.read_csv(IN_CSV)
print(f"Loaded {len(df)} samples\n")

rif_r = df[df['RIF'] == 'R']
print(f"=== Among RIF-resistant samples (n={len(rif_r)}), co-resistance to other drugs ===")
co_rows = []
for d in DRUGS:
    if d == 'RIF':
        continue
    sub = rif_r[rif_r[d].isin(['R', 'S'])]
    if len(sub) == 0:
        continue
    co_r = int((sub[d] == 'R').sum())
    pct = round(100 * co_r / len(sub), 1)
    co_rows.append({"drug": d, "n_tested": len(sub), "co_resistant": co_r, "co_resistance_pct": pct})
    print(f"  {d}: {co_r}/{len(sub)} ({pct:.1f}%)")

co_df = pd.DataFrame(co_rows).sort_values("co_resistance_pct", ascending=False)
co_df.to_csv(f"{OUT_DIR}/co_resistance_with_RIF.csv", index=False)
print(f"\nSaved: {OUT_DIR}/co_resistance_with_RIF.csv")

prev = pd.read_csv(PREV_CSV)
plt.figure(figsize=(10, 5))
plt.bar(prev['drug'], prev['resistance_pct'], color='steelblue')
plt.ylabel("Resistance (%)")
plt.title("Per-drug resistance prevalence")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/drug_prevalence_chart.png", dpi=150)
plt.close()
print(f"Saved chart: {OUT_DIR}/drug_prevalence_chart.png")

if 'isolation_country' in df.columns:
    top_countries = df['isolation_country'].value_counts().head(10)
    plt.figure(figsize=(10, 5))
    top_countries.plot(kind='bar', color='seagreen')
    plt.ylabel("Sample count")
    plt.title("Top 10 countries by sample count")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/country_distribution_chart.png", dpi=150)
    plt.close()
    print(f"Saved chart: {OUT_DIR}/country_distribution_chart.png")

mdr_mask = (df['RIF'] == 'R') & (df['INH'] == 'R')
fq_r  = df[['MXF', 'LEV', 'OFX']].eq('R').any(axis=1)
inj_r = df[['KAN', 'AMK', 'CAP']].eq('R').any(axis=1)
xdr_mask = mdr_mask & fq_r & inj_r

plt.figure(figsize=(6, 5))
plt.bar(['Total cohort', 'MDR-TB', 'XDR-TB'],
        [len(df), int(mdr_mask.sum()), int(xdr_mask.sum())],
        color=['steelblue', 'orange', 'crimson'])
plt.ylabel("Sample count")
plt.title("MDR-TB and XDR-TB burden in cohort")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/mdr_xdr_summary_chart.png", dpi=150)
plt.close()
print(f"Saved chart: {OUT_DIR}/mdr_xdr_summary_chart.png")

print(f"\nDone. All tables and charts saved in {OUT_DIR}/")
