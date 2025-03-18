from Bio import SeqIO
import numpy as np
import sys
import os

def quality_check_fasta(input_file):
    sequences = {}  # Dictionary to store unique sequences
    discarded_counts = {
        "duplicate": 0, "invalid_start": 0, "invalid_stop": 0,
        "premature_stop": 0, "not_divisible_by_3": 0,
        "short_length": 0, "invalid_chars": 0, "outliers": 0
    }
    
    valid_seqs = []

    for record in SeqIO.parse(input_file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq).upper()

        # Remove duplicates (keep only the first occurrence)
        if seq_id in sequences:
            discarded_counts["duplicate"] += 1
            continue
        sequences[seq_id] = sequence

        # Check for standard criteria
        if not sequence.startswith("ATG"):
            discarded_counts["invalid_start"] += 1
            continue
        if not sequence[-3:] in {"TAA", "TAG", "TGA"}:
            discarded_counts["invalid_stop"] += 1
            continue
        if "*" in sequence[3:-3]:  # Check for premature stop codons
            discarded_counts["premature_stop"] += 1
            continue
        if len(sequence) % 3 != 0:
            discarded_counts["not_divisible_by_3"] += 1
            continue
        if len(sequence) < 300:
            discarded_counts["short_length"] += 1
            continue
        if any(char not in "ATGC" for char in sequence):
            discarded_counts["invalid_chars"] += 1
            continue

        valid_seqs.append(record)

    # Identify outliers in length
    lengths = np.array([len(record.seq) for record in valid_seqs])
    mean_len = np.mean(lengths)
    std_len = np.std(lengths)

    filtered_seqs = [
        record for record in valid_seqs if mean_len - 2 * std_len <= len(record.seq) <= mean_len + 2 * std_len
    ]
    discarded_counts["outliers"] = len(valid_seqs) - len(filtered_seqs)

    # Output filenames
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    passed_file = f"{base_name}_passed.fasta"
    trimmed_file = f"{base_name}_trimmed.fasta"

    # Save the quality-checked sequences
    SeqIO.write(filtered_seqs, passed_file, "fasta")

    # Save the trimmed sequences (without stop codon)
    for record in filtered_seqs:
        record.seq = record.seq[:-3]  # Remove stop codon
    SeqIO.write(filtered_seqs, trimmed_file, "fasta")

    # Print summary
    print(f"File: {input_file}")
    for key, count in discarded_counts.items():
        print(f" - {key.replace('_', ' ').capitalize()} discarded: {count}")
    print(f"Final sequences: {len(filtered_seqs)}")

    return passed_file, trimmed_file


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 SeqQualityChecker.py <input_fasta>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    quality_check_fasta(input_fasta)
