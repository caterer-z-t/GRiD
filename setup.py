# setup.py at the root of ~/epi/lpa
from setuptools import setup, find_packages

setup(
    name="GRiD",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "pysam",
    ],
    entry_points={
        "console_scripts": [
            "grid-count=GRiD.cli:main",  # optional if you want CLI
        ],
    },
)
