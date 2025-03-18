#!/bin/bash

# Define the Python script name
PYTHON_SCRIPT="SeqQC.py"

# Check if the Python script exists
if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "Error: $PYTHON_SCRIPT not found!"
    exit 1
fi

# Create an output directory if it doesn't exist
mkdir -p QC_Passed

# Process all .txt and .fasta files
for file in *.txt *.fasta; do
    # Skip if no matching files are found
    [[ -e "$file" ]] || continue

    # Extract filename without extension
    BASENAME=$(basename "$file" | sed 's/\.[^.]*$//')

    # Output file names
    OUTPUT_PASSED="QC_Passed/${BASENAME}_passed.fasta"
    OUTPUT_TRIMMED="QC_Passed/${BASENAME}_trimmed.fasta"

    echo "Processing: $file -> $OUTPUT_PASSED, $OUTPUT_TRIMMED"

    # Run the Python script
    python3 "$PYTHON_SCRIPT" "$file" "$OUTPUT_PASSED" "$OUTPUT_TRIMMED"

    # Check if Python script executed successfully
    if [[ $? -eq 0 ]]; then
        echo "Completed: $file"
    else
        echo "Error processing: $file"
    fi
done

echo "All files processed!"
