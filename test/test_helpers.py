import threading
import pytest
from pathlib import Path

from grid.utils.helper_dir.write_result_to_file import write_result_to_file
from grid.utils.helper_dir.find_all_cram_files import find_cram_files
from grid.utils.helper_dir.setup_output_file import setup_output_file
from grid.utils.helper_dir.create_region import create_region_string
from grid.utils.helper_dir.display_results import (
    print_individual_success,
    print_individual_error,
)


# --- write_result_to_file ---

def test_write_result_to_file_writes_line(tmp_path):
    f = tmp_path / "out.tsv"
    f.write_text("")
    lock = threading.Lock()
    write_result_to_file(f, "S1", 42, lock)
    assert "S1\t42\n" in f.read_text()

def test_write_result_to_file_thread_safe(tmp_path):
    f = tmp_path / "out.tsv"
    f.write_text("")
    lock = threading.Lock()
    threads = [threading.Thread(target=write_result_to_file, args=(f, f"S{i}", i, lock)) for i in range(10)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    lines = f.read_text().strip().split("\n")
    assert len(lines) == 10

def test_write_result_to_file_error_string(tmp_path):
    f = tmp_path / "out.tsv"
    f.write_text("")
    lock = threading.Lock()
    write_result_to_file(f, "S1", "Error", lock)
    assert "S1\tError\n" in f.read_text()


# --- find_cram_files ---

def test_find_cram_files_finds_crams(tmp_path):
    (tmp_path / "a.cram").touch()
    (tmp_path / "b.cram").touch()
    (tmp_path / "c.bam").touch()
    result = find_cram_files(str(tmp_path))
    assert len(result) == 2
    assert all(r.endswith(".cram") for r in result)

def test_find_cram_files_empty_dir(tmp_path):
    assert find_cram_files(str(tmp_path)) == []


# --- helper_dir setup_output_file ---

def test_helper_setup_output_file(tmp_path):
    out = tmp_path / "sub" / "result.tsv"
    path = setup_output_file(str(out), "chr6", 1000, 2000)
    assert path.exists()
    assert path.read_text() == "Sample\tchr6:1000-2000\n"


# --- create_region_string (helper_dir version) ---

def test_helper_create_region_string_with_coords():
    result = create_region_string(None, "chr6", 1000, 2000)
    assert result == "chr6:1000-2000"

def test_helper_create_region_string_full_region():
    result = create_region_string("chr6:1000-2000", None, None, None)
    assert result == "chr6:1000-2000"


# --- display_results ---

def test_print_individual_success_no_console(capsys):
    print_individual_success("S1", "indexed")
    # rich prints to its own console; just check no exception is raised

def test_print_individual_error_no_console(capsys):
    print_individual_error("S1", "file not found")
    # just check no exception is raised

def test_print_individual_success_with_console():
    from unittest.mock import MagicMock
    mock_console = MagicMock()
    print_individual_success("S1", "done", progress_console=mock_console)
    mock_console.print.assert_called_once()
    assert "S1" in mock_console.print.call_args[0][0]

def test_print_individual_error_with_console():
    from unittest.mock import MagicMock
    mock_console = MagicMock()
    print_individual_error("S1", "something went wrong", progress_console=mock_console)
    mock_console.print.assert_called_once()
    assert "S1" in mock_console.print.call_args[0][0]
