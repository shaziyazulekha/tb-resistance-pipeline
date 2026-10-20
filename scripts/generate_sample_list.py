import csv, os
METADATA = "data/filtered_metadata.csv"
OUTPUT   = "data/slurm_sample_list.tsv"
os.makedirs("data", exist_ok=True)
rows = []
with open(METADATA, newline="") as f:
    for row in csv.DictReader(f):
        iso = str(row.get("isolate_id","")).strip()
        sra = str(row.get("sra_accession","")).strip()
        for ch in ["[","]","'",' ','"']: sra = sra.replace(ch,"")
        sra = sra.split(",")[0].strip()
        if iso and sra and sra not in ("","nan","None"):
            rows.append((iso,sra))
with open(OUTPUT,"w") as f:
    for iso,sra in rows:
        f.write(f"{iso}\t{sra}\n")
print(f"Written: {len(rows)} samples → {OUTPUT}")
print(f"Use in sbatch: --array=0-{len(rows)-1}%100")
