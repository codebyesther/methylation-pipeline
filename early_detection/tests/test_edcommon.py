import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from edcommon import load_count_file, parse_cgi_name, subject_of, window_of  # noqa: E402
from step_2_merge_runs import _keyword_exact  # noqa: E402


@pytest.mark.parametrize("name,expected", [
    ("HL01_12m", "pre12"),
    ("256960_12-months", "pre12"),
    ("HL01_6m", "pre6"),
    ("256960_6-months", "pre6"),
    ("HL06_6m", "pre6"),
    ("HL01_near-Dx", "near"),
    ("256960_at-cancer", "near"),
    ("CTRL03_control", "control"),
    ("INNOV_control_01", "control"),
    ("HL16_baseline", None),
    ("HL06_unlabeled", None),
])
def test_window_of(name, expected):
    assert window_of(name) == expected


def test_subject_of_uses_id_list_first():
    assert subject_of("256960_12-months", ["256960"]) == "256960"
    assert subject_of("HL01_12m") == "HL01"


def test_parse_cgi_name_genes_and_index():
    p = parse_cgi_name("CGI_chr1_778604_779167_LOC100288069_0")
    assert (p["chr"], p["start"], p["end"], p["genes"], p["cgi_index"]) == (
        "chr1", 778604, 779167, ["LOC100288069"], "0")
    p = parse_cgi_name("CGI_chr3_100_900_GATA2_GATA2-AS1_2")
    assert p["genes"] == ["GATA2", "GATA2-AS1"]
    p = parse_cgi_name("CGI_chr7_100_900_1")
    assert p["genes"] == [] and p["cgi_index"] == "1"


def test_keyword_exact_keeps_glob20_and_globmin80_apart():
    assert _keyword_exact("RunA_Tot10_Glob20_cpgi.xlsx", "Glob20")
    assert not _keyword_exact("RunA_Tot10_Glob200_cpgi.xlsx", "Glob20")
    assert not _keyword_exact("RunA_Tot10_GlobMin80_cpgi.xlsx", "Glob20")


def test_load_count_file_reads_total_row_and_ignores_info(tmp_path):
    df = pd.DataFrame({
        "Header": ["Info_Pool", "Total CpG island fragments counts for this particular spreadsheet",
                   "CGI_chr1_1_2_A_0", "CGI_chr1_3_4_B_1"],
        "S1_12m": ["pool1", 10, 4, 6],
        "S2_control": ["pool1", 5, 5, 0],
    })
    path = tmp_path / "x_Glob20.xlsx"
    df.to_excel(path, index=False)
    cgi, totals = load_count_file(str(path))
    assert list(cgi.index) == ["CGI_chr1_1_2_A_0", "CGI_chr1_3_4_B_1"]
    assert cgi.loc["CGI_chr1_3_4_B_1", "S2_control"] == 0
    assert totals.to_dict() == {"S1_12m": 10, "S2_control": 5}
