import pytest
from pathlib import Path

from grid.utils.hi_inference import read_dip_cn_file, read_neighbors_file


# --- read_dip_cn_file ---

def test_read_dip_cn_file_basic(tmp_path):
    f = tmp_path / "dipcn.tsv"
    f.write_text("1001 2.5\n1002 3.0\n1003 1.8\n")
    IRRs, IDs, IDtoInd = read_dip_cn_file(str(f), [], [], {})
    assert IDs == [1001, 1002, 1003]
    assert IRRs == pytest.approx([2.5, 3.0, 1.8])
    assert IDtoInd[1001] == 0
    assert IDtoInd[1003] == 2

def test_read_dip_cn_file_skips_blank_lines(tmp_path):
    f = tmp_path / "dipcn.tsv"
    f.write_text("1001 2.5\n\n1002 3.0\n")
    IRRs, IDs, _ = read_dip_cn_file(str(f), [], [], {})
    assert len(IDs) == 2

def test_read_dip_cn_file_empty(tmp_path):
    f = tmp_path / "dipcn.tsv"
    f.write_text("")
    IRRs, IDs, IDtoInd = read_dip_cn_file(str(f), [], [], {})
    assert IDs == []
    assert IRRs == []


# --- read_neighbors_file ---

def _make_neighbors_file(tmp_path, lines):
    f = tmp_path / "nbrs.txt"
    f.write_text("ID\thap\tnbr_ind\tcM_len\tcM_edge\tid_nbr\thap_nbr\n" + "\n".join(lines) + "\n")
    return str(f)

def test_read_neighbors_file_basic(tmp_path):
    # IDs: 1001→0, 1002→1; hap_nbrs has 4 slots (2 per sample)
    IDtoInd = {1001: 0, 1002: 1}
    hap_nbrs = [[] for _ in range(4)]
    f = _make_neighbors_file(tmp_path, [
        "1001\t1\t0\t1.5\t0.1\t1002\t2",  # sample 0 hap1 gets neighbor (2*1+2-1=3)
    ])
    result = read_neighbors_file(f, IDtoInd, hap_nbrs, MAX_NBR=10)
    assert result[0] == [3]  # 2*j + hap_nbr - 1 = 2*1+2-1 = 3

def test_read_neighbors_file_skips_negative_ids(tmp_path):
    IDtoInd = {1001: 0, 1002: 1}
    hap_nbrs = [[] for _ in range(4)]
    f = _make_neighbors_file(tmp_path, [
        "-1\t1\t0\t1.0\t0.0\t1002\t1",
    ])
    result = read_neighbors_file(f, IDtoInd, hap_nbrs, MAX_NBR=10)
    assert all(lst == [] for lst in result)

def test_read_neighbors_file_respects_max_nbr(tmp_path):
    IDtoInd = {1001: 0, 1002: 1, 1003: 2}
    hap_nbrs = [[] for _ in range(6)]
    lines = [
        "1001\t1\t0\t1.0\t0.0\t1002\t1",
        "1001\t1\t0\t1.0\t0.0\t1003\t1",
    ]
    f = _make_neighbors_file(tmp_path, lines)
    result = read_neighbors_file(f, IDtoInd, hap_nbrs, MAX_NBR=1)
    assert len(result[0]) == 1  # capped at MAX_NBR=1

def test_read_neighbors_file_skips_unknown_ids(tmp_path):
    IDtoInd = {1001: 0}
    hap_nbrs = [[] for _ in range(2)]
    f = _make_neighbors_file(tmp_path, [
        "1001\t1\t0\t1.0\t0.0\t9999\t1",  # 9999 not in IDtoInd
    ])
    result = read_neighbors_file(f, IDtoInd, hap_nbrs, MAX_NBR=10)
    assert result[0] == []

def test_read_neighbors_file_empty_after_header(tmp_path):
    IDtoInd = {1001: 0}
    hap_nbrs = [[] for _ in range(2)]
    f = tmp_path / "nbrs.txt"
    f.write_text("ID\thap\tnbr_ind\tcM_len\tcM_edge\tid_nbr\thap_nbr\n")
    result = read_neighbors_file(str(f), IDtoInd, hap_nbrs, MAX_NBR=10)
    assert result == [[], []]
