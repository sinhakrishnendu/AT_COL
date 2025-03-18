🔹 Key Features
✅ Parallel execution (2 processes at a time to maintain thermal homeostasis).
✅ IQR-based length filtering for outlier detection.
✅ Removes duplicate sequences while keeping only one unique copy.
✅ Converts sequences to uppercase to avoid errors.
✅ Integrates PRANK → IQ-TREE2 → CODEML seamlessly.
✅ Runs all selection tests (site, branch, branch-site models).
✅ Applies LRT and BH correction to determine statistically significant selection.
✅ Outputs final results in codeml_results.xlsx.

🔹 Expected Output
📂 Quality-checked sequences

{species}_passed.fasta
{species}_trimmed.fasta
📂 Alignments

{species}_output.best.fas
📂 Phylogenetic trees

{species}_output.best.fas.treefile
📂 Selection analysis results

codeml_results.xlsx (Final report)

🚀 How to Run
Ensure dependencies are installed:
bash
Copy
Edit
pip3 install biopython scipy pandas numpy
Give execution permissions:
bash
Copy
Edit
chmod +x run_pipeline.sh
Run the pipeline:
bash
Copy
Edit
./run_pipeline.sh
This will automate everything from sequence QC to final selection analysis in one command. 🚀