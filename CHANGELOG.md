# Pipeline debugging log

- prefetch (conda build) failed silently; resolved by using the official
  sratoolkit 3.2.0 binary distribution directly instead of the conda package
- fasterq-dump --split-files desynced paired-end reads; switched to --split-3
- Trimmomatic OOM on paired FASTQ trimming; fixed with _JAVA_OPTIONS="-Xmx8g"
- GATK4 CLI syntax error on java options; fixed to --java-options "$JAVA_MEM"
- Node-local scratch is not shared across array tasks; switched temp
  directory to a path under the shared home filesystem
- run_batch.sh and submit_tb_array_final.sh both applied an OFFSET value,
  double-counting it and silently processing wrong/duplicate samples for
  OFFSET > 0; fixed by applying OFFSET in exactly one place
- GATK step errors were redirected to /dev/null, hiding real failures;
  corrected to log to per-sample files
- GATK crashed with a JNI UnsupportedClassVersionError from invoking system
  Java 8 instead of the conda environment's Java 17; fixed by exporting the
  conda env's bin/ onto PATH at the top of the script
- Sample list parsing with cut -f1 broke on a space-delimited file; fixed
  to awk '{print $1}'
- Home quota hard-capped with no grace period; quota exhaustion from
  undeleted intermediate FASTQ/BAM files caused silent job failures under
  high concurrency
- A subset of filtered VCFs used chromosome naming incompatible with
  TB-Profiler; resolved by renaming contigs (bcftools annotate
  --rename-chrs) into 05_renamed_vcfs/ before profiling (finish_remaining.sh)
- Download steps (prefetch, fasterq-dump) also redirected errors to
  /dev/null, in addition to the GATK steps already fixed; this meant
  download failures left no trace in SLURM logs either. All four
  redirects now write to per-sample log files under logs/.
