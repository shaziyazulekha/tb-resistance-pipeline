#!/bin/bash
export PATH="$HOME/miniconda2/envs/tb_pipeline/bin:$PATH"
SCRATCH=$PROJECT_ROOT/tb_pipeline/tmp
set -euo pipefail
mkdir -p $SCRATCH
#SBATCH --array=1000-1999%100
#SBATCH --job-name=tb_vcf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --partition=standard
#SBATCH --output=$PROJECT_ROOT/tb_pipeline/logs/tb_%A_%a.out
#SBATCH --error=$PROJECT_ROOT/tb_pipeline/logs/tb_%A_%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<redacted>@cdac.in

ENV=$PROJECT_ROOT/tb_pipeline/sratoolkit.3.2.0-centos_linux64/bin
export PATH="$ENV:$PATH"
BWA=$PROJECT_ROOT/miniconda2/envs/tb_pipeline/bin/bwa
SAMTOOLS=$PROJECT_ROOT/miniconda2/envs/tb_pipeline/bin/samtools
GATK=$PROJECT_ROOT/miniconda2/envs/tb_pipeline/bin/gatk
TRIMMOMATIC=$PROJECT_ROOT/miniconda2/envs/tb_pipeline/bin/trimmomatic
FASTQC=$PROJECT_ROOT/miniconda2/envs/tb_pipeline/bin/fastqc
PREFETCH=$PROJECT_ROOT/tb_pipeline/sratoolkit.3.2.0-centos_linux64/bin/prefetch
FASTERQDUMP=$ENV/fasterq-dump

WORK=$PROJECT_ROOT/tb_pipeline
REF=$WORK/reference/H37Rv.fasta
ADAPTERS=$WORK/trimmomatic_adapters/All_Illumina_Adapters.fa
SAMPLE_LIST=$PROJECT_ROOT/tb_pipeline/data/missing_sample_list.tsv
VCF_DIR=$WORK/05_filtered_vcfs
QC_DIR=$WORK/01_raw_QC
PROGRESS=$WORK/data/slurm_progress.csv
THREADS=16
JAVA_MEM="-Xmx20g"

ACTUAL_IDX=$(( ${OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
SAMPLE=$(awk -v i=$((ACTUAL_IDX+1)) 'NR==i{print $1}' $SAMPLE_LIST)
SRA=$(awk -v i=$((ACTUAL_IDX+1)) 'NR==i{print $2}' $SAMPLE_LIST)

if [ -z "$SAMPLE" ] || [ -z "$SRA" ]; then
    echo "ERROR: No sample at index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

FINAL_VCF=$VCF_DIR/${SAMPLE}_filtered.vcf
LOG="[$(date '+%H:%M:%S')] [$SAMPLE]"

echo "$LOG START | SRA: $SRA | Node: ${SLURMD_NODENAME:-unknown}"

if [ -f "$FINAL_VCF" ] && [ -s "$FINAL_VCF" ]; then
    echo "$LOG Already done, skipping."
    exit 0
fi

SCRATCH=$PROJECT_ROOT/tb_pipeline/tmp/${SAMPLE}_${SLURM_JOB_ID}
mkdir -p $SCRATCH
cleanup() { rm -rf $SCRATCH 2>/dev/null; }
trap cleanup EXIT

echo "$LOG [1/7] Downloading $SRA..."
for attempt in 1 2 3; do
    $PREFETCH $SRA -O $SCRATCH --max-size 50G >> logs/${SAMPLE}_prefetch.log 2>&1
    if [ -d "$SCRATCH/$SRA" ]; then
        $FASTERQDUMP $SCRATCH/$SRA -O $SCRATCH \
            --split-files --force --threads $THREADS >> logs/${SAMPLE}_fasterqdump.log 2>&1
        if ls $SCRATCH/*.fastq 2>/dev/null | grep -q .; then
            for fq in $SCRATCH/*.fastq; do gzip -1 "$fq" & done
            wait
            rm -rf $SCRATCH/$SRA
            echo "$LOG Download OK attempt $attempt"
            break
        fi
    fi
    echo "$LOG Attempt $attempt failed, retrying..."
    sleep $((attempt * 10))
done

FQ1=$(ls $SCRATCH/*_1.fastq.gz 2>/dev/null | head -1)
FQ2=$(ls $SCRATCH/*_2.fastq.gz 2>/dev/null | head -1)
[ -z "$FQ1" ] && FQ1=$(ls $SCRATCH/*.fastq.gz 2>/dev/null | head -1)

if [ -z "$FQ1" ]; then
    echo "$LOG FAILED: no FASTQ"
    echo "${SAMPLE},${SRA},failed,,$(date -Iseconds),download_failed" >> $PROGRESS
    exit 1
fi

[ -n "$FQ1" ] && [ -n "$FQ2" ] && READ_TYPE=paired || READ_TYPE=single
echo "$LOG Read type: $READ_TYPE"

echo "$LOG [2/7] FastQC..."
mkdir -p $QC_DIR
if [ "$READ_TYPE" = "paired" ]; then
    $FASTQC "$FQ1" "$FQ2" -o $QC_DIR -t $THREADS -q 2>/dev/null || true
else
    $FASTQC "$FQ1" -o $QC_DIR -t $THREADS -q 2>/dev/null || true
fi

echo "$LOG [3/7] Trimmomatic..."
TRIM1=$SCRATCH/${SAMPLE}_1_paired.fastq.gz
TRIM2=$SCRATCH/${SAMPLE}_2_paired.fastq.gz
TRIM_SE=$SCRATCH/${SAMPLE}_se.fastq.gz

if [ "$READ_TYPE" = "paired" ]; then
    _JAVA_OPTIONS="-Xmx8g" $TRIMMOMATIC PE -threads $THREADS "$FQ1" "$FQ2" \
        "$TRIM1" $SCRATCH/${SAMPLE}_1_up.gz \
        "$TRIM2" $SCRATCH/${SAMPLE}_2_up.gz \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
        -quiet 2>/dev/null
    rm -f $SCRATCH/${SAMPLE}_1_up.gz $SCRATCH/${SAMPLE}_2_up.gz
else
    _JAVA_OPTIONS="-Xmx8g" $TRIMMOMATIC SE -threads $THREADS "$FQ1" "$TRIM_SE" \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
        -quiet 2>/dev/null
fi
rm -f "$FQ1" "$FQ2"
echo "$LOG FASTQs deleted after trimming."

echo "$LOG [4/7] BWA MEM..."
RG="@RG\tID:${SAMPLE}\tLB:lib1\tPL:ILLUMINA\tSM:${SAMPLE}\tPU:unit1"
SAM=$SCRATCH/${SAMPLE}.sam
BAM=$SCRATCH/${SAMPLE}_sorted.bam

if [ "$READ_TYPE" = "paired" ]; then
    $BWA mem -t $THREADS -R "$RG" $REF "$TRIM1" "$TRIM2" > $SAM 2>/dev/null
else
    $BWA mem -t $THREADS -R "$RG" $REF "$TRIM_SE" > $SAM 2>/dev/null
fi
rm -f "$TRIM1" "$TRIM2" "$TRIM_SE"

$SAMTOOLS sort -@ $THREADS -o $BAM $SAM 2>/dev/null
$SAMTOOLS index $BAM 2>/dev/null
rm -f $SAM
echo "$LOG SAM and trimmed FASTQs deleted."

echo "$LOG [5/7] MarkDuplicates..."
DEDUP=$SCRATCH/${SAMPLE}_dedup.bam
METRICS=$SCRATCH/${SAMPLE}_metrics.txt
$GATK --java-options "$JAVA_MEM" MarkDuplicates \
    -I $BAM -O $DEDUP -M $METRICS \
    --REMOVE_DUPLICATES true \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY LENIENT \
    >>logs/${SAMPLE}_gatk_debug.log 2>&1
rm -f $BAM ${BAM}.bai
echo "$LOG Sorted BAM deleted."

echo "$LOG [6/7] HaplotypeCaller..."
RAW_VCF=$SCRATCH/${SAMPLE}_raw.vcf
$GATK --java-options "$JAVA_MEM" HaplotypeCaller \
    -R $REF \
    -I $DEDUP \
    -O $RAW_VCF \
    -ploidy 1 \
    --native-pair-hmm-threads $THREADS \
    --minimum-mapping-quality 20 \
    --min-base-quality-score 20 \
    >> logs/${SAMPLE}_haplotypecaller.log 2>&1
rm -f $DEDUP ${DEDUP%.bam}.bai $METRICS
echo "$LOG Dedup BAM deleted."

echo "$LOG [7/7] VariantFiltration..."
mkdir -p $VCF_DIR
$GATK VariantFiltration \
    -R $REF \
    -V $RAW_VCF \
    -O $FINAL_VCF \
    --filter-name "LowQD"   --filter-expression "QD < 2.0"  \
    --filter-name "HighFS"  --filter-expression "FS > 60.0" \
    --filter-name "LowMQ"   --filter-expression "MQ < 40.0" \
    --filter-name "HighSOR" --filter-expression "SOR > 3.0" \
    >> logs/${SAMPLE}_variantfiltration.log 2>&1
rm -f $RAW_VCF ${RAW_VCF}.idx

if [ ! -s "$FINAL_VCF" ]; then
    echo "$LOG FAILED: empty VCF"
    echo "${SAMPLE},${SRA},failed,,$(date -Iseconds),empty_vcf" >> $PROGRESS
    exit 1
fi

VARIANTS=$(grep -vc "^#" $FINAL_VCF 2>/dev/null || echo 0)
echo "$LOG SUCCESS: $VARIANTS variants → $FINAL_VCF"
(flock 200; echo "${SAMPLE},${SRA},done,${FINAL_VCF},$(date -Iseconds),${VARIANTS}") \
    200>>${PROGRESS}.lock >> $PROGRESS
echo "$LOG DONE"
