import os
import glob
import pandas as pd

# Define the folder to search
output_dir = "output"

# Search for both .xlsx and .csv files that contain 'gene_annotation' in their name
candidate_files = glob.glob(os.path.join(output_dir, "*gene_annotation*.xlsx")) + \
                  glob.glob(os.path.join(output_dir, "*gene_annotation*.csv"))

# Debug info if no matches
if not candidate_files:
    available_files = os.listdir(output_dir)
    print("❌ No gene annotation file found.")
    print("📂 Files in 'output/' folder:")
    for f in available_files:
        print(" -", f)
    raise FileNotFoundError("No Excel or CSV file with 'gene_annotation' in filename was found in the 'output/' folder.")

# Use the first matching file
gene_annotation_file = candidate_files[0]
print(f"📄 Found gene annotation file: {gene_annotation_file}")

# Load the file
if gene_annotation_file.endswith('.xlsx'):
    gene_annot = pd.read_excel(gene_annotation_file)
elif gene_annotation_file.endswith('.csv'):
    gene_annot = pd.read_csv(gene_annotation_file)
else:
    raise ValueError("Unsupported file format. Only .xlsx or .csv are supported.")

# Get all 'Gene' columns
gene_cols = [col for col in gene_annot.columns if col.startswith("Gene")]
if not gene_cols:
    raise ValueError("No columns starting with 'Gene' found in the annotation file.")

# Melt into long format
gene_annot_long = gene_annot.melt(
    id_vars=["chr", "start genomic coordinate", "end genomic coordinate"],
    value_vars=gene_cols,
    var_name="gene_col",
    value_name="gene_name"
)
gene_annot_long = gene_annot_long.dropna(subset=["gene_name"])

# Build 'cgi_id' in chr:start-end format
gene_annot_long["cgi_id"] = gene_annot_long.apply(
    lambda row: f"{row['chr']}:{int(row['start genomic coordinate'])}-{int(row['end genomic coordinate'])}",
    axis=1
)

# Keep only what's needed
gene_annot_final = gene_annot_long[["cgi_id", "gene_name"]].drop_duplicates()

# Show and save
print(gene_annot_final.head())
gene_annot_final.to_csv(os.path.join(output_dir, "gene_cgi_map.csv"), index=False)
print("✅ Saved: gene_cgi_map.csv")
