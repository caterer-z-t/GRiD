"""
find_neighbors tests.

Internal helper functions (filter_regions_by_variance, find_neighbors_sklearn,
save_neighbors, read_normalized_data) were consolidated into the top-level
find_neighbors() function and the C++ find_neighbors executable in the current
implementation. Integration-level tests are in test_pipeline.py.

write_normalized_output / read_normalized_data round-trip tests remain in
test_normalize.py which covers normalize_mosdepth internals.
"""

import pytest

# No importable unit-level symbols remain from find_neighbors.py beyond
# the top-level find_neighbors() orchestrator, which requires a full
# pipeline config and running mosdepth outputs. Mark as placeholder.


def test_find_neighbors_placeholder():
    """Placeholder — unit tests for internal helpers removed when the
    implementation was consolidated into find_neighbors() + C++ binary.
    See test_pipeline.py for integration coverage."""
    pass
