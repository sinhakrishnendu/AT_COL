This pair of programme SeqQC.py and SeqQCPipeline.sh are paired. bash command ./SeqQCPipeline.sh will fetch the .py file and quality check codon sequences based on criteria like:
starts with start codon
stops with stop codon 
no intermediate stop codon
length divisible by 3
no non-nucleotide characters
making final output in all caps
check for length outlier and discards if found
check for duplicate sequence based on ids and discards if found