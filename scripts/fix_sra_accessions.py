import pandas as pd
import requests
import time

df = pd.read_csv('data/filtered_metadata.csv', dtype=str)
print(f"Total samples: {len(df)}")

# Find ERS accessions
ers_mask = df['sra_accession'].str.startswith('ERS', na=False)
print(f"ERS accessions to convert: {ers_mask.sum()}")
print(f"ERR/SRR already correct  : {(~ers_mask).sum()}")

converted = 0
failed = 0

for idx, row in df[ers_mask].iterrows():
    ers = row['sra_accession']
    try:
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ers}&result=read_run&fields=run_accession&format=json"
        r = requests.get(url, timeout=30)
        data = r.json()
        if data:
            err = data[0]['run_accession']
            df.at[idx, 'sra_accession'] = err
            converted += 1
            print(f"  {ers} → {err}", end="\r")
        else:
            failed += 1
    except Exception as e:
        failed += 1
    time.sleep(0.2)

df.to_csv('data/filtered_metadata.csv', index=False)
print(f"\n✅ Done! Converted: {converted}  Failed: {failed}")
print(f"Saved → data/filtered_metadata.csv")
