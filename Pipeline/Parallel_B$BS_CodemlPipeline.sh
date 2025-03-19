#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define the base control file
base_ctl_file="codeml.ctl"

# Check if the control file exists
if [ ! -f "$base_ctl_file" ]; then
    echo "Error: Base control file '$base_ctl_file' not found."
    exit 1
fi

# Extract seqfile from the base control file
seqfile=$(grep -E '^ *seqfile *= *' "$base_ctl_file" | sed -E 's/.*= *([^ ]*).*/\1/')
seqfile="${seqfile//\"/}"  # Trim possible quotes around the seqfile path

# Check if the sequence file exists
if [ ! -f "$seqfile" ]; then
    echo "Error: Sequence file '$seqfile' not found."
    exit 1
fi

# Function to run codeml for a given treefile
run_codeml() {
    treefile="$1"
    base_name="${treefile%.treefile}"
    
    # ----- Branch Model -----
    branch_folder="${base_name}_B"
    mkdir -p "$branch_folder"
    branch_ctl_file="$branch_folder/$base_name.B.ctl"
    cp "$base_ctl_file" "$branch_ctl_file"
    sed -i "s|^ *treefile *=.*|treefile = $treefile|" "$branch_ctl_file"
    sed -i "s|^ *model *=.*|model = 2|" "$branch_ctl_file"
    sed -i "s|^ *NSsites *=.*|NSsites = 0|" "$branch_ctl_file"
    cp "$seqfile" "$treefile" "$branch_folder/"
    (cd "$branch_folder" && codeml "$base_name.B.ctl")
    
    # ----- Branch-Site Model -----
    bs_folder="${base_name}_BS"
    mkdir -p "$bs_folder"
    bs_ctl_file="$bs_folder/$base_name.BS.ctl"
    cp "$base_ctl_file" "$bs_ctl_file"
    sed -i "s|^ *treefile *=.*|treefile = $treefile|" "$bs_ctl_file"
    sed -i "s|^ *model *=.*|model = 2|" "$bs_ctl_file"
    sed -i "s|^ *NSsites *=.*|NSsites = 2|" "$bs_ctl_file"
    cp "$seqfile" "$treefile" "$bs_folder/"
    (cd "$bs_folder" && codeml "$base_name.BS.ctl")
    
    # ----- Null Model -----
    null_folder="${base_name}_BS_NULL"
    mkdir -p "$null_folder"
    null_ctl_file="$null_folder/$base_name.BS_NULL.ctl"
    cp "$base_ctl_file" "$null_ctl_file"
    sed -i "s|^ *treefile *=.*|treefile = $treefile|" "$null_ctl_file"
    sed -i "s|^ *model *=.*|model = 2|" "$null_ctl_file"
    sed -i "s|^ *NSsites *=.*|NSsites = 2|" "$null_ctl_file"
    sed -i "s|^ *fix_omega *=.*|fix_omega = 1|" "$null_ctl_file"
    sed -i "s|^ *omega *=.*|omega = 1|" "$null_ctl_file"
    cp "$seqfile" "$treefile" "$null_folder/"
    (cd "$null_folder" && codeml "$base_name.BS_NULL.ctl")
    
}

# Export function for parallel execution
export -f run_codeml
export base_ctl_file seqfile

# Run codeml in parallel with a maximum of 4 jobs at a time
parallel -j 8 run_codeml ::: *.treefile

echo "All processing completed."