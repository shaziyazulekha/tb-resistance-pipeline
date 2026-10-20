import pandas as pd
from pathlib import Path

meta = pd.read_csv('data/filtered_metadata.csv')
drugs = ['RIF','INH','PZA','EMB','STR','MXF','LEV','OFX','KAN','AMK','CAP','ETH']

vcf_ids = set(p.name.replace('_filtered.vcf.gz','').replace('_filtered.vcf','')
              for p in Path('05_filtered_vcfs').glob('*_filtered.vcf*'))
print(f"Completed VCFs right now: {len(vcf_ids)}")

sub = meta[meta['isolate_id'].isin(vcf_ids)]
print(f"Matched to metadata: {len(sub)}\n")

print(f"{'Drug':<6}{'FullPop_R':>10}{'FullPop_S':>10}{'FullPop_Cov%':>13}   |  {'VCF_R':>7}{'VCF_S':>7}{'VCF_Cov%':>10}")
for d in drugs:
    fr=(meta[d]=='R').sum(); fs=(meta[d]=='S').sum(); fc=100*(fr+fs)/len(meta)
    vr=(sub[d]=='R').sum(); vs=(sub[d]=='S').sum(); vc=100*(vr+vs)/len(sub) if len(sub) else 0
    print(f"{d:<6}{fr:>10}{fs:>10}{fc:>12.1f}%   |  {vr:>7}{vs:>7}{vc:>9.1f}%")
