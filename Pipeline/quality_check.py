import sys
import re
import numpy as np
from Bio import SeqIO

def is_valid_sequence(seq):
    """Checks if a sequence is valid for analysis (case insensitive)."""
    seq = seq.upper()  # Convert everything to uppercase
    return (seq.startswith("ATG") and
            seq[-3:] in ["TAA", "TAG", "TGA"] and
            len(seq) % 3 == 0 and
            re.fullmatch("[ATGC]+", seq) and  # Ensures only ATGC are present
            "*" not in seq[:-3])  # Ensures no premature stop codon

def remove_duplicates(records):
    """Removes duplicate sequences by key name."""
    seen = set()
    unique_records = []
    for record in records:
        if record.id not in seen:
            seen.add(record.id)
            unique_records.append(record)
    return unique_records

def filter_length_outliers(records):
    """Removes sequences that are extreme length outliers using IQR method."""
    lengths = np.array([len(record.seq) for record in records])
    
    # Compute Q1, Q3, and IQR
    Q1 = np.percentile(lengths, 25)
    Q3 = np.percentile(lengths, 75)
    IQR = Q3 - Q1

    # Define outlier bounds
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    # Filter out sequences outside IQR range
    filtered_records = [rec for rec in records if lower_bound <= len(rec.seq) <= upper_bound]
    
    print(f"Filtered length outliers: {len(records) - len(filtered_records)} sequences removed.")
    return filtered_records

# Read input and output file names from command-line arguments
input_file, output_passed, output_trimmed = sys.argv[1:4]

# Load sequences
records = list(SeqIO.parse(input_file, "fasta"))

# Step 1: Convert all sequences to uppercase
for rec in records:
    rec.seq = rec.seq.upper()

# Step 2: Validate sequences
valid_records = [rec for rec in records if is_valid_sequence(str(rec.seq))]

# Step 3: Remove duplicate sequences
valid_records = remove_duplicates(valid_records)

# Step 4: Remove length outliers using IQR
valid_records = filter_length_outliers(valid_records)

# Save passed sequences
SeqIO.write(valid_records, output_passed, "fasta")

# Create trimmed version (remove stop codon)
for rec in valid_records:
    rec.seq = rec.seq[:-3]  # Trim stop codon
SeqIO.write(valid_records, output_trimmed, "fasta")

print(f"Final sequences: {len(valid_records)}")
