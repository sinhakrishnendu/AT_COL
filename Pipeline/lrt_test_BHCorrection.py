import os
import pandas as pd
import numpy as np
from scipy.stats import chi2

# Get the current working directory
directory = os.getcwd()
print(f"Checking directory: {directory}")

def benjamini_hochberg(p_values):
    """Apply Benjamini-Hochberg FDR correction to a list of p-values."""
    m = len(p_values)  # Number of tests
    sorted_indices = np.argsort(p_values)  # Sort indices by p-values
    sorted_p_values = np.array(p_values)[sorted_indices]
    adjusted_p_values = np.zeros(m)

    min_value = 1
    for i in range(m, 0, -1):
        rank = i
        min_value = min(min_value, (m / rank) * sorted_p_values[i - 1])
        adjusted_p_values[i - 1] = min_value

    # Rearrange adjusted p-values back to original order
    corrected_p_values = np.zeros(m)
    corrected_p_values[sorted_indices] = adjusted_p_values
    return corrected_p_values

# Find all CSV files in the directory
csv_files = [f for f in os.listdir(directory) if f.endswith(".csv")]
print(f"CSV files found: {csv_files}")

if not csv_files:
    print("No CSV files found. Exiting.")
    exit()

for filename in csv_files:
    file_path = os.path.join(directory, filename)
    df = pd.read_csv(file_path)

    # Check if required columns exist
    required_cols = {"Folder", "lnL"}
    if not required_cols.issubset(df.columns):
        print(f"Skipping {filename}: Missing required columns {required_cols}")
        continue

    # Extract gene names by removing suffixes (_B, _BS, _BS_NULL, _M0)
    df["Gene"] = df["Folder"].str.replace(r'(_B|_BS|_BS_NULL|_M0)$', '', regex=True)
    unique_genes = df["Gene"].unique()

    # Perform LRT for Branch-site Model (_BS vs _BS_NULL)
    branchsite_results = []
    for gene in unique_genes:
        bs_data = df[df["Folder"].str.contains(f"{gene}_BS", regex=False)]
        null_data = df[df["Folder"].str.contains(f"{gene}_BS_NULL", regex=False)]

        if not bs_data.empty and not null_data.empty:
            lnL_bs = bs_data.iloc[0]['lnL']
            lnL_null = null_data.iloc[0]['lnL']
            lrt_value = 2 * (lnL_bs - lnL_null)
            p_value = 1 - chi2.cdf(lrt_value, df=1)

            branchsite_results.append({
                "Gene": gene,
                "lnL_BS": lnL_bs,
                "lnL_BS_NULL": lnL_null,
                "LRT": lrt_value,
                "p_value": p_value
            })

    # Convert to DataFrame
    branchsite_df = pd.DataFrame(branchsite_results)

    # Apply Benjamini-Hochberg correction for Branch-site Model
    if not branchsite_df.empty:
        branchsite_df['FDR_Corrected_P'] = benjamini_hochberg(branchsite_df["p_value"].values)

    # Perform LRT for Branch Model (_B vs M0)
    branch_model_results = []
    for gene in unique_genes:
        branch_data = df[df['Folder'].str.contains(f'{gene}_B', regex=False)]
        m0_data = df[df['Folder'].str.contains(f'{gene}_M0', regex=False)]

        if not branch_data.empty and not m0_data.empty:
            lnL_branch = branch_data.iloc[0]['lnL']
            lnL_m0 = m0_data.iloc[0]['lnL']
            lrt_value_branch = 2 * (lnL_branch - lnL_m0)
            p_value_branch = chi2.sf(lrt_value_branch, df=1)

            branch_model_results.append({
                "Gene": gene,
                "lnL_B": lnL_branch,
                "lnL_M0": lnL_m0,
                "LRT": lrt_value_branch,
                "p_value": p_value_branch
            })

    # Convert to DataFrame
    branch_model_df = pd.DataFrame(branch_model_results)

    # Apply Benjamini-Hochberg correction for Branch Model
    if not branch_model_df.empty:
        branch_model_df["BH_FDR"] = benjamini_hochberg(branch_model_df["p_value"].values)

    # Save results to an Excel file
    output_file = f"LRT_results_{filename.replace('.csv', '.xlsx')}"
    with pd.ExcelWriter(output_file) as writer:
        if not branchsite_df.empty:
            branchsite_df.to_excel(writer, sheet_name="Branchsite_Model", index=False)
        if not branch_model_df.empty:
            branch_model_df.to_excel(writer, sheet_name="Branch_Model", index=False)

    print(f"Results saved to {output_file}")

print("LRT analysis with Benjamini-Hochberg correction completed.")
