#!/bin/bash

# Ensure PRANK and IQ-TREE2 are installed
if ! command -v prank &> /dev/null || ! command -v iqtree2 &> /dev/null; then
    echo "Error: PRANK or IQ-TREE2 is not installed. Please install them in WSL."
    exit 1
fi

# Limit the number of concurrent jobs
MAX_JOBS=4
CURRENT_JOBS=0

# Process each text file containing FASTA sequences
for txt_file in *.txt; do
    [[ -e "$txt_file" ]] || { echo "No text files found!"; exit 1; }

    # Extract the base name without extension
    base_name="${txt_file%.txt}"

    echo "Processing: $txt_file -> Quality Check"

    # Run Quality Checker
    passed_fasta="${base_name}_passed.fasta"
    trimmed_fasta="${base_name}_trimmed.fasta"
    python3 SeqQualityChecker.py "$txt_file"

    # Run PRANK alignment in the background
    prank -d="$passed_fasta" -codon -o="${base_name}_prank" &

    # Track the number of running jobs
    ((CURRENT_JOBS++))

    # Wait for jobs if the limit is reached
    if [[ "$CURRENT_JOBS" -ge "$MAX_JOBS" ]]; then
        wait
        CURRENT_JOBS=0
    fi
done

# Wait for PRANK jobs to finish
wait

# Process PRANK outputs with IQ-TREE2
for prank_output in *_prank.best.fas; do
    [[ -e "$prank_output" ]] || continue

    echo "Running IQ-TREE2 on: $prank_output"

    # Run IQ-TREE2 in the background
    iqtree2 -s "$prank_output" -st CODON -T AUTO -B 10000 -alrt 10000 -bnni &

    # Track running jobs
    ((CURRENT_JOBS++))

    # Wait for jobs if the limit is reached
    if [[ "$CURRENT_JOBS" -ge "$MAX_JOBS" ]]; then
        wait
        CURRENT_JOBS=0
    fi
done

# Wait for remaining jobs
wait

echo "All files processed!"
