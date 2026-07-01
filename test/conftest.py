"""
Stub out heavy optional dependencies that aren't needed for pure-function tests.

pysam  — only used in create_index_for_file(); not needed by any tested function.
"""
import sys
from unittest.mock import MagicMock

# Stub pysam before any grid module is imported
if "pysam" not in sys.modules:
    sys.modules["pysam"] = MagicMock()
