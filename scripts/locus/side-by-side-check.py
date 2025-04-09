import pandas as pd

# Load files
matrix_path = "output/scaled_fragment_ratios_matrix.xlsx"
map_path = "output/gene_cgi_map.csv"

# Load CpG island names from matrix
df = pd.read_excel(matrix_path)
start_idx = df[df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
matrix_cpg = df.iloc[start_idx:, 0].astype(str).str.strip().unique()

# Load gene mapping
map_df = pd.read_csv(map_path)
map_df.columns = map_df.columns.str.strip()
map_df["cgi_id"] = map_df["cgi_id"].astype(str).str.strip()
map_cpg = map_df["cgi_id"].unique()

# Print comparison
print("Example CpG islands from methylation matrix:")
print(matrix_cpg[:5])

print("\nExample CpG islands from mapping file:")
print(map_cpg[:5])

# Check exact overlap
overlap = set(matrix_cpg) & set(map_cpg)
print(f"\nNumber of CpG islands in BOTH files: {len(overlap)}")
