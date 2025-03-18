#!/bin/bash

# Ensure PRANK and IQ-TREE2 are installed
if ! command -v prank &> /dev/null || ! command -v iqtree2 &> /dev/null; then
    echo "Error: PRANK or IQ-TREE2 is not installed. Please install them in WSL."
    exit 1
fi

# Process each FASTA file in the directory
for fasta_file in *.fasta; do
    # Skip if no FASTA files found
    [[ -e "$fasta_file" ]] || { echo "No FASTA files found!"; exit 1; }

    # Extract the base name without extension
    base_name="${fasta_file%.fasta}"

    echo "Processing: $fasta_file -> PRANK -> IQ-TREE2"

    # Run PRANK alignment
    prank -d="$fasta_file" -codon -o="${base_name}_prank"

    # Find the main output file from PRANK (usually ends with .best.fas)
    prank_output="${base_name}_prank.best.fas"

    # Check if PRANK output was generated
    if [[ ! -f "$prank_output" ]]; then
        echo "Error: PRANK failed for $fasta_file. Skipping IQ-TREE2."
        continue
    fi

    echo "Running IQ-TREE2 on: $prank_output"

    # Run IQ-TREE2 on the PRANK output
    iqtree2 -s "$prank_output" -st CODON -T AUTO -B 10000 -alrt 10000 -bnni

    echo "Completed: $fasta_file"
done

echo "All files processed!"
