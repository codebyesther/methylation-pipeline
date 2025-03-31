import pandas as pd
import os

def find_excel_file(directory, keyword):
    for file_name in os.listdir(directory):
        if keyword in file_name and file_name.endswith(".xlsx"):
            return os.path.join(directory, file_name)
    raise FileNotFoundError(f"No Excel file with keyword '{keyword}' found in directory '{directory}'")

def generate_gene_annotation(input_excel, output_csv):
    # Load the Excel file without skipping rows
    df_raw = pd.read_excel(input_excel, sheet_name=0, header=None)

    # Extract CGI names starting from row 1, column 0
    cgi_names_clean = df_raw.iloc[1:, 0].dropna().astype(str)
    cgi_names_clean = cgi_names_clean[cgi_names_clean.str.startswith("CGI_")]

    # Split and process each CGI name
    split_cgi = cgi_names_clean.str.split("_").tolist()
    processed_rows = []
    for parts in split_cgi:
        if parts[0] == "CGI":
            parts = parts[1:]

        if parts and parts[-1].isdigit():
            probe_id = parts.pop()
        else:
            probe_id = ""

        chr_part = parts[0] if len(parts) > 0 else ""
        start = parts[1] if len(parts) > 1 else ""
        end = parts[2] if len(parts) > 2 else ""
        genes = parts[3:] if len(parts) > 3 else []

        row = [chr_part, start, end] + genes + [probe_id]
        processed_rows.append(row)

    # Determine max number of gene columns
    max_genes = max(len(row) - 4 for row in processed_rows)

    # Create column names
    column_names = ["chr", "start genomic coordinate", "end genomic coordinate"]
    gene_columns = [f"Gene{i+1}" for i in range(max_genes)]
    final_columns = column_names + gene_columns + ["CGI index or probe ID"]

    # Normalize row lengths
    normalized_rows = [
        row[:3] + row[3:-1] + [""] * (max_genes - len(row[3:-1])) + [row[-1]]
        for row in processed_rows
    ]

    # Build and save DataFrame
    final_df = pd.DataFrame(normalized_rows, columns=final_columns)
    final_df.to_csv(output_csv, index=False)

# Example usage
if __name__ == "__main__":
    input_excel = find_excel_file("../output", "matrix")
    output_csv = "../output/structured_gene_annotation.csv"
    generate_gene_annotation(input_excel, output_csv)
