"""
normalize_mosdepth tests.

Internal helper functions (norm_chrom, load_repeat_mask, normalize_matrix,
select_high_variance_regions, build_matrix_from_regions, filter_empty_samples,
map_mosdepth_files_to_samples, find_bed_gz_for_individual,
write_normalized_output) were consolidated into the top-level
normalize_mosdepth() orchestrator in the current implementation.

Integration-level tests are in test_pipeline.py.
"""

import pytest


def test_normalize_mosdepth_placeholder():
    """Placeholder — unit tests for internal helpers removed when the
    implementation was consolidated into normalize_mosdepth().
    See test_pipeline.py for integration coverage."""
    pass
