#!/bin/bash
LISTFILE=$1
TOTAL=$(wc -l < "$LISTFILE")
CHUNK=1000
OFFSET=0
while [ $OFFSET -lt $TOTAL ]; do
    REMAIN=$((TOTAL - OFFSET))
    SIZE=$CHUNK
    [ $REMAIN -lt $CHUNK ] && SIZE=$REMAIN
    sbatch --array=1-${SIZE}%300 --export=ALL,LISTFILE=$LISTFILE,OFFSET=$OFFSET run_tbprofiler_array.sh
    OFFSET=$((OFFSET + CHUNK))
done
