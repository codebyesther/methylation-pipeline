import os
import glob
import pandas as pd

# Define the folder to search
output_dir = "output"

# Automatically find the Excel gene annotation file
gene_annotation_files = glob.glob(os.path.join(output_dir, "*gene_annotation*.xlsx"))
if not gene_annotation_files:
    raise FileNotFoundError("No Excel file with 'gene_annotation' found in the 'output/' folder.")
gene_annotation_file = gene_annotation_files[0]
print(f"📄 Found gene annotation file: {gene_annotation_file}")

# Load the gene annotation file
gene_annot = pd.read_excel(gene_annotation_file)

# Find all columns that start with 'Gene'
gene_cols = [col for col in gene_annot.columns if col.startswith("Gene")]
if not gene_cols:
    raise ValueError("No columns starting with 'Gene' found in the annotation file.")

# Melt the gene table into long format
gene_annot_long = gene_annot.melt(
    id_vars=["chr", "start genomic coordinate", "end genomic coordinate"],
    value_vars=gene_cols,
    var_name="gene_col",
    value_name="gene_name"
)

# Drop rows with missing gene names
gene_annot_long = gene_annot_long.dropna(subset=["gene_name"])

# Build 'cgi_id' as 'chr:start-end'
gene_annot_long["cgi_id"] = gene_annot_long.apply(
    lambda row: f"{row['chr']}:{int(row['start genomic coordinate'])}-{int(row['end genomic coordinate'])}",
    axis=1
)

# Final clean DataFrame
gene_annot_final = gene_annot_long[["cgi_id", "gene_name"]].drop_duplicates()

# Show results
print(gene_annot_final.head())

# Optionally: save to CSV
gene_annot_final.to_csv(os.path.join(output_dir, "gene_cgi_map.csv"), index=False)
print("✅ Saved: gene_cgi_map.csv")
