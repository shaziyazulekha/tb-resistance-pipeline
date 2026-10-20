#!/bin/bash
#SBATCH --job-name=tb_vcf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --partition=standard
#SBATCH --reservation=nitk_res
#SBATCH --output=$PROJECT_ROOT/tb_pipeline/logs/tb_%A_%a.out
#SBATCH --error=$PROJECT_ROOT/tb_pipeline/logs/tb_%A_%a.err
#SBATCH --array=0-999%100

bash $PROJECT_ROOT/tb_pipeline/submit_tb_array_final.sh
