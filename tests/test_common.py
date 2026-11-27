import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from common import (GLOB20_TOTAL_ROW, GLOBMIN80_TOTAL_ROW, SCALED_TOTAL_ROW,  # noqa: E402
                    cgi_gene_pairs, cgi_genes, classify_timepoint, contains_id,
                    detailed_timepoint, find_file, gene_matrix, get_patient,
                    read_scaled_matrix)


def test_patient_ids_match_whole_tokens():
    ids = ["PT01", "PT010"]
    assert get_patient("PT010_C2D1", ids) == "PT010"
    assert get_patient("PT01_Baseline", ids) == "PT01"
    assert get_patient("LOI_PT01_C7D1-Off-tx", ids) == "PT01"
    assert get_patient("OTHER99_Baseline", ids) is None
    assert contains_id("INNOV_control_01", "INNOV")


@pytest.mark.parametrize("name,grouped,detailed", [
    ("PT001_Baseline", "Baseline", "Baseline"),
    ("PT001_C2D1", "On-Treatment", "C2"),
    ("PT001_C1D15", "On-Treatment", "C1"),
    ("PT001_C10D1", "On-Treatment", "C10"),
    ("PT001_C7D1-Off-tx", "Post-Treatment", "Off-Tx"),
    ("PT001_C7D1_off_tx", "Post-Treatment", "Off-Tx"),
    ("INNOV_control_01", "Healthy", "Healthy"),
])
def test_timepoints(name, grouped, detailed):
    assert classify_timepoint(name) == grouped
    assert detailed_timepoint(name) == detailed


def test_cgi_genes_exact():
    assert cgi_genes("CGI_chr3_100_900_GATA2_GATA2-AS1_2") == ["GATA2", "GATA2-AS1"]
    assert cgi_genes("CGI_chrX_100_900_1") == []
    pairs = cgi_gene_pairs(["CGI_chr3_1_2_GATA2-AS1_0", "CGI_chr3_3_4_GATA2_1"])
    assert pairs.groupby("gene_name")["cgi_id"].count().to_dict() == {"GATA2": 1, "GATA2-AS1": 1}


def test_gene_matrix_does_not_mix_similar_symbols():
    m = pd.DataFrame({"S1": [10, 5, 1]},
                     index=["CGI_chr3_1_2_GATA2_0", "CGI_chr3_3_4_GATA2-AS1_1", "Total row"])
    g = gene_matrix(m)
    assert g.loc["GATA2", "S1"] == 10 and g.loc["GATA2-AS1", "S1"] == 5


def test_find_file_ignores_case_and_lock_files(tmp_path):
    (tmp_path / "~$Patient ID list.xlsx").write_text("")
    (tmp_path / "Patient ID list.xlsx").write_text("")
    assert os.path.basename(find_file(str(tmp_path), "patient")) == "Patient ID list.xlsx"


def test_read_scaled_matrix_reads_labeled_row_and_drops_inf(tmp_path):
    df = pd.DataFrame({
        "Header": [SCALED_TOTAL_ROW, GLOB20_TOTAL_ROW, GLOBMIN80_TOTAL_ROW, "CGI_chr1_1_2_A_0", "CGI_chr1_3_4_B_1"],
        "S1": [500.0, 5, 1000, "INF", 20.0],
    })
    path = tmp_path / "scaled_fragment_ratios_matrix.xlsx"
    df.to_excel(path, index=False)
    cgi, scaled = read_scaled_matrix(str(path))
    assert scaled["S1"] == 500.0
    assert np.isnan(cgi.loc["CGI_chr1_1_2_A_0", "S1"])   # "INF" -> missing, not infinity
    assert cgi.loc["CGI_chr1_3_4_B_1", "S1"] == 20.0


def test_read_scaled_matrix_computes_ratio_for_older_files(tmp_path):
    df = pd.DataFrame({"Header": [GLOB20_TOTAL_ROW, GLOBMIN80_TOTAL_ROW, "CGI_chr1_1_2_A_0"],
                       "S1": [5, 1000, 1.0]})
    path = tmp_path / "scaled_fragment_ratios_matrix.xlsx"
    df.to_excel(path, index=False)
    _, scaled = read_scaled_matrix(str(path))
    assert scaled["S1"] == pytest.approx(500.0)
