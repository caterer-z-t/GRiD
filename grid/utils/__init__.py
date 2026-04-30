from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("GRiD")
except PackageNotFoundError:
    __version__ = "unknown"

__author__ = "Zachary Caterer"
__email__ = "ztcaterer@colorado.edu"
