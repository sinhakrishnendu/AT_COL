import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd

def is_valid_sequence(seq):
    """Checks if the sequence meets quality criteria."""
    seq = seq.upper()  # Convert to uppercase
    if not seq.startswith("ATG") or len(seq) < 300 or len(seq) % 3 != 0:
        return False
    if any(base not in "ATGC" for base in seq):  # Non-standard nucleotides
        return False
    if "*" in seq[:-3]:  # Internal stop codon check (except at the end)
        return False
    if seq[-3:] not in ["TAA", "TAG", "TGA"]:  # Valid stop codon check
        return False
    return True

def remove_length_outliers(sequences):
    """Removes length outliers using the IQR method."""
    lengths = np.array([len(seq) for seq in sequences.values()])
    q1, q3 = np.percentile(lengths, [25, 75])
    iqr = q3 - q1
    lower_bound, upper_bound = q1 - 1.5 * iqr, q3 + 1.5 * iqr

    return {id_: seq for id_, seq in sequences.items() if lower_bound <= len(seq) <= upper_bound}

def process_fasta(input_file, output_passed, output_trimmed):
    """Processes the FASTA file to filter sequences based on quality criteria."""
    
    unique_sequences = {}  # Store unique sequences by ID

    for record in SeqIO.parse(input_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq).upper()  # Convert sequence to uppercase

        if seq_id not in unique_sequences:  # Ensure uniqueness by ID
            unique_sequences[seq_id] = seq
    
    print(f"Total unique sequences after duplicate removal: {len(unique_sequences)}")

    # Filter sequences based on quality checks
    valid_sequences = {id_: seq for id_, seq in unique_sequences.items() if is_valid_sequence(seq)}
    print(f"Sequences passing quality check: {len(valid_sequences)}")

    # Remove outliers
    filtered_sequences = remove_length_outliers(valid_sequences)
    print(f"Sequences after outlier removal: {len(filtered_sequences)}")

    # Save sequences that passed QC
    with open(output_passed, "w") as passed_fasta:
        for id_, seq in filtered_sequences.items():
            passed_fasta.write(f">{id_}\n{seq}\n")

    # Save sequences with stop codon removed
    with open(output_trimmed, "w") as trimmed_fasta:
        for id_, seq in filtered_sequences.items():
            trimmed_fasta.write(f">{id_}\n{seq[:-3]}\n")  # Remove last 3 bases (stop codon)

    # Print summary of discarded sequences
    discarded = len(unique_sequences) - len(filtered_sequences)
    print(f"Total sequences discarded: {discarded}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and process FASTA sequences.")
    parser.add_argument("input", help="Input FASTA file (.txt)")
    parser.add_argument("output_passed", help="Output file for sequences that passed QC")
    parser.add_argument("output_trimmed", help="Output file for sequences without stop codon")

    args = parser.parse_args()
    process_fasta(args.input, args.output_passed, args.output_trimmed)
