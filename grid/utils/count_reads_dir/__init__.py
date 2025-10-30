# grid/utils/count_reads_dir/__init__.py
"""
Count reads utilities.

This package contains modular functions for counting reads
in directories and processing read count results.
"""

from .count_reads_in_region import count_reads_in_region
from .process_single_cram import process_single_cram


__all__ = [
    'count_reads_in_region',
    'process_single_cram'
]
