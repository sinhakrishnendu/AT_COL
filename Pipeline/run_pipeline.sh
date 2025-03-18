#!/bin/bash

# Run quality check
for file in *.txt; do
    base_name=$(basename "$file" .txt)
    python3 quality_check.py "$file" "${base_name}_passed.fasta" "${base_name}_trimmed.fasta"
done

# Run PRANK alignment (2 parallel executions)
find *_passed.fasta | xargs -P 2 -I {} prank -d={} -codon -o={}_output

# Run IQ-TREE2 (2 parallel executions)
find *_output.best.fas | xargs -P 2 -I {} iqtree2 -s {} -st CODON -T AUTO -B 10000 -alrt 10000 -bnni

# Run CODEML for selection analysis
for tree in *_output.best.fas.treefile; do
    alignment="${tree%.treefile}"
    python3 run_codeml.py "$alignment" "$tree"
done

echo "Pipeline execution complete!"
