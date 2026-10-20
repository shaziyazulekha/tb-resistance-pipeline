#!/bin/bash
#SBATCH --job-name=tbprof_arr
#SBATCH --output=logs_array/tbprof_%A_%a.log
#SBATCH --reservation=nitk_res
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

source ~/.bashrc
micromamba activate tbprofiler

mkdir -p 05_renamed_vcfs tbprofiler_results logs_array

LINE=$((SLURM_ARRAY_TASK_ID + OFFSET))
VCF=$(sed -n "${LINE}p" "$LISTFILE")
[ -z "$VCF" ] && exit 0

sample=$(basename "$VCF" _filtered.vcf.gz)
[ -f "tbprofiler_results/results/${sample}.results.json" ] && exit 0

renamed="05_renamed_vcfs/${sample}_filtered.vcf.gz"
bcftools annotate --rename-chrs chr_rename.txt -Oz -o "$renamed" "$VCF"
bcftools index -t "$renamed"
tb-profiler vcf_profile "$renamed" --dir tbprofiler_results 
