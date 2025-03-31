#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob

class Args:
    patients = ""
    methylation = ""
    outdir = "plots"
args = Args()

# Locate the files based on the specified criteria
data_folder = "data"
output_folder = "output"

patient_files = glob.glob(os.path.join(data_folder, '*[pP][aA][tT][iI][eE][nN][tT] [iI][dD]*.xlsx'))
matrix_files = glob.glob(os.path.join(output_folder, '*[mM][aA][tT][rR][iI][xX]*.xlsx'))

if patient_files:
    args.patients = patient_files[0]
else:
    raise FileNotFoundError("No patient ID file found in the data folder.")

if matrix_files:
    args.methylation = matrix_files[0]
else:
    raise FileNotFoundError("No matrix file found in the output folder.")

methylation_dfs = {}
patient_ids = []

os.makedirs(args.outdir, exist_ok=True)

# Read files
if "patient" in args.patients.lower():
    patient_df = pd.read_excel(args.patients)
    patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

df = pd.read_excel(args.methylation)
methylation_dfs[os.path.basename(args.methylation)] = df

# Metadata functions
def get_patient(sample):
    for pid in patient_ids:
        if pid in sample:
            return pid
    return None

def normalize_timepoint(sample):
    if "Baseline" in sample:
        return "Baseline"
    elif "Off-tx" in sample:
        return "Post-Treatment"
    elif "INNOV" in sample:
        return "Healthy"
    else:
        return "On-Treatment"

# Process the methylation data
for fname, df in methylation_dfs.items():
    print(f"\n=== Processing file: {fname} ===")
    base_fname = os.path.splitext(fname)[0]  # For cleaner filenames

    # Locate locus-level matrix
    start_idx = df[df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
    cpg_island_df = df.iloc[start_idx:].reset_index(drop=True)
    cpg_island_df.rename(columns={cpg_island_df.columns[0]: "CpG_Island"}, inplace=True)
    cpg_island_df.dropna(how="all", subset=cpg_island_df.columns[1:], inplace=True)

    samples = cpg_island_df ▋
