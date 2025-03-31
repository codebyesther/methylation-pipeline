# gene_heatmap_lineplots_local.py

import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# Argument parser for local execution
parser = argparse.ArgumentParser(description='Generate gene-level methylation heatmaps and line plots.')
parser.add_argument('--output_dir', type=str, default='output', help='Directory to save plots')
args = parser.parse_args()

# Create output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

# Helper to find file containing a keyword
def find_file(directory, keyword):
    files = glob.glob(os.path.join(directory, f"*{keyword}*"))
    if not files:
        raise FileNotFoundError(f"No file containing '{keyword}' found in directory '{directory}'")
    return files[0]

# Define the directories
output_folder = 'output'
data_folder = 'data'

# Load inputs based on known filename patterns
try:
    cpg_matrix_file = find_file(output_folder, "matrix")
    patient_list_file = find_file(data_folder, "patient")
    gene_annotation_file = find_file(output_folder, "gene_annotation")
except FileNotFoundError as e:
    print(e)
    exit(1)

# Specify the encoding to handle decoding issues
try:
    cpg_matrix = pd.read_csv(cpg_matrix_file, sep="\t", index_col=0, encoding='ISO-8859-1')
except pd.errors.ParserError as e:
    print(f"Error parsing {cpg_matrix_file}: {e}")
    exit(1)

try:
    patient_df = pd.read_excel(patient_list_file)
except Exception as e:
    print(f"Error reading {patient_list_file}: {e}")
    exit(1)

try:
    gene_annot = pd.read_csv(gene_annotation_file, encoding='ISO-8859-1')
except Exception as e:
    print(f"Error reading {gene_annotation_file}: {e}")
    exit(1)

# Function to classify sample timepoints
def classify_timepoint(sample_name):
    if "Baseline" in sample_name:
        return "Baseline"
    elif "Off-tx" in sample_name:
        return "Post-Treatment"
    elif "INNOV" in sample_name:
        return "Healthy"
    else:
        return "On-Treatment"

# Annotate sample timepoints
sample_timepoints = {s: classify_timepoint(s) for s in cpg_matrix.columns}
timepoint_df = pd.DataFrame.from_dict(sample_timepoints, orient='index', columns=['Timepoint'])

# Filter for genes with multiple CpG islands
cpg_gene_counts = gene_annot['gene_name'].value_counts()
multicpg_genes = cpg_gene_counts[cpg_gene_counts > 1].index.tolist()

# Average methylation for each gene
gene_means = {}
for gene in multicpg ▋
