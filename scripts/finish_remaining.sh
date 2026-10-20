#!/bin/bash
#SBATCH --job-name=tbprof_final
#SBATCH --output=logs_final/tbprof_final_%j.log
#SBATCH --reservation=nitk_res
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

source ~/.bashrc
micromamba activate tbprofiler

mkdir -p 05_renamed_vcfs tbprofiler_results

while read sample; do
    vcf="05_filtered_vcfs/${sample}_filtered.vcf.gz"
    renamed="05_renamed_vcfs/${sample}_filtered.vcf.gz"
    bcftools annotate --rename-chrs chr_rename.txt -Oz -o "$renamed" "$vcf"
    bcftools index -t "$renamed"
    tb-profiler vcf_profile "$renamed" --dir tbprofiler_results
done < missing_samples.txt
