Step	Action
Step 1	Runs SeqQC4RunPhyloPipelineFromRawSeq2Tree.py on each .txt file to remove duplicates & filter sequences.
Step 2	Runs PRANK on the filtered FASTA files (*_passed.fasta).
Step 3	Runs IQ-TREE2 on the PRANK alignment (*_prank.best.fas).
Step 4	Uses parallel execution (4 files at a time) to improve efficiency.
Step 5	Outputs filtered, aligned, and phylogenetic tree results.