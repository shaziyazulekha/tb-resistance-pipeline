#!/bin/bash
#SBATCH --job-name=tbprofiler
#SBATCH --output=tbprofiler_%j.log
#SBATCH --reservation=nitk_res
#SBATCH --partition=standard
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

source ~/.bashrc
micromamba activate tbprofiler

mkdir -p tbprofiler_results

for vcf in 05_filtered_vcfs/*_filtered.vcf.gz; do
    size=$(stat -c%s "$vcf")
    if [ "$size" -lt 100 ]; then
        echo "SKIPPING corrupted: $vcf"
        continue
    fi
    sample=$(basename "$vcf" _filtered.vcf.gz)
    if [ -f "tbprofiler_results/results/${sample}.results.json" ]; then
        echo "SKIPPING already done: $sample"
        continue
    fi
    tb-profiler vcf_profile "$vcf" --dir tbprofiler_results
done
