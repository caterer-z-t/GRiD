from .create_region import create_region_string
from .display_results import print_individual_success, print_individual_error
from .find_all_cram_files import find_cram_files
from .load_flags_from_yaml import load_flags
from .setup_output_file import setup_output_file
from .write_result_to_file import write_result_to_file

__all__ = [
    "create_region_string",
    "print_individual_success",
    "print_individual_error",
    "find_cram_files",
    "load_flags",
    "setup_output_file",
    "write_result_to_file",
]


from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("GRiD")
except PackageNotFoundError:
    __version__ = "unknown"

__author__ = "Zachary Caterer"
__email__ = "ztcaterer@colorado.edu"
