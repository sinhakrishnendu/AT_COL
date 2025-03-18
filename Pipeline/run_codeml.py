import os
import subprocess
import pandas as pd
import numpy as np
from Bio import Phylo

# Define models
SITE_MODELS = {
    "M0": "model = 0\nNSsites = 0",
    "M1a": "model = 0\nNSsites = 1",
    "M2a": "model = 0\nNSsites = 2",
    "M7": "model = 0\nNSsites = 7",
    "M8": "model = 0\nNSsites = 8"
}

BRANCH_MODELS = {
    "One-ratio": "model = 0",
    "Free-ratio": "model = 1",
    "Two-ratio": "model = 2"
}

BRANCH_SITE_MODEL = "model = 2\nNSsites = 2"

def create_codeml_control(alignment, tree, model_name, model_content, out_file):
    """Generates a CODEML control file dynamically."""
    ctl_content = f"""
    seqfile = {alignment}
    treefile = {tree}
    outfile = {out_file}
    noisy = 2
    verbose = 1
    seqtype = 1
    ndata = 1
    icode = 0
    cleandata = 0
    CodonFreq = 7
    {model_content}
    fix_omega = 0
    omega = 0.4
    RateAncestor = 2
    """
    ctl_file = f"{model_name}.ctl"
    with open(ctl_file, "w") as f:
        f.write(ctl_content.strip())
    return ctl_file

def run_codeml(ctl_file):
    """Runs CODEML with the specified control file."""
    subprocess.run(["codeml", ctl_file], check=True)

def extract_log_likelihood(codeml_output):
    """Extracts log likelihood values from CODEML output."""
    with open(codeml_output) as f:
        for line in f:
            if "lnL" in line and "np:" in line:
                return float(line.split()[4])  # Extracts log-likelihood value
    return None

def likelihood_ratio_test(l1, l2, df):
    """Computes the likelihood ratio test statistic and p-value."""
    from scipy.stats import chi2
    lr_stat = 2 * (l2 - l1)
    p_value = 1 - chi2.cdf(lr_stat, df)
    return lr_stat, p_value

def benjamini_hochberg(p_values):
    """Applies the Benjamini-Hochberg correction."""
    p_values = np.array(p_values)
    ranked_pvals = np.argsort(p_values)
    adjusted_pvals = np.zeros_like(p_values, dtype=float)
    m = len(p_values)

    for i, p_idx in enumerate(ranked_pvals):
        adjusted_pvals[p_idx] = min((p_values[p_idx] * m) / (i + 1), 1)

    return adjusted_pvals

def run_pipeline(alignment, tree):
    """Runs site, branch, and branch-site models."""
    results = []

    # 1. Run Site Models
    for model, content in SITE_MODELS.items():
        ctl_file = create_codeml_control(alignment, tree, model, content, f"{model}.out")
        run_codeml(ctl_file)
        log_likelihood = extract_log_likelihood(f"{model}.out")
        results.append((model, log_likelihood))

    # 2. Run Branch Models
    for model, content in BRANCH_MODELS.items():
        ctl_file = create_codeml_control(alignment, tree, model, content, f"{model}.out")
        run_codeml(ctl_file)
        log_likelihood = extract_log_likelihood(f"{model}.out")
        results.append((model, log_likelihood))

    # 3. Run Branch-Site Model for Each Branch
    tree_obj = Phylo.read(tree, "newick")
    for branch in tree_obj.get_terminals() + tree_obj.get_nonterminals():
        branch_name = branch.name if branch.name else "internal"
        tree_out = f"branch_{branch_name}.tree"
        Phylo.write(tree_obj, tree_out, "newick")

        ctl_file = create_codeml_control(alignment, tree_out, f"branch_{branch_name}", BRANCH_SITE_MODEL, f"branch_{branch_name}.out")
        run_codeml(ctl_file)
        log_likelihood = extract_log_likelihood(f"branch_{branch_name}.out")
        results.append((f"Branch-Site_{branch_name}", log_likelihood))

    # Perform LRT and BH correction
    df = pd.DataFrame(results, columns=["Model", "LogLikelihood"])
    df["LRT"], df["p-value"] = zip(*[likelihood_ratio_test(df.iloc[i, 1], df.iloc[i+1, 1], 2) for i in range(0, len(df)-1, 2)])
    df["BH-corrected p"] = benjamini_hochberg(df["p-value"].values)

    # Save to Excel
    df.to_excel("codeml_results.xlsx", index=False)
    print("Analysis complete. Results saved in codeml_results.xlsx")

# Run pipeline
if __name__ == "__main__":
    import sys
    alignment, tree = sys.argv[1:3]
    run_pipeline(alignment, tree)
