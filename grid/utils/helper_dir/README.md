# Helper Utilities

Shared utility functions used across multiple GRiD modules. These provide common functionality for file handling, configuration management, and output formatting.

## Modules

### `create_region.py`
**Purpose:** Format genomic region strings for samtools/pysam  
**Functions:**
- `create_region_string(region, chrom, start, end)` - Creates standardized region strings (e.g., "chr6:160000000-160100000")

**Usage:**
```python
from grid.utils.helper_dir.create_region import create_region_string

# From separate components
region = create_region_string(None, "chr6", 160000000, 160100000)
# Returns: "chr6:160000000-160100000"

# From existing region string
region = create_region_string("chr6:160000000-160100000", None, None, None)
# Returns: "chr6:160000000-160100000"
```

---

### `display_results.py`
**Purpose:** Pretty-print results to console using Rich library  
**Functions:**
- `display_results(results, title)` - Display tabular results with formatting

**Usage:**
```python
from grid.utils.helper_dir.display_results import display_results

results = [
    {"sample": "S001", "count": 1234},
    {"sample": "S002", "count": 2345}
]
display_results(results, "Read Counts")
```

---

### `find_all_cram_files.py`
**Purpose:** Discover CRAM files in directories  
**Functions:**
- `find_all_cram_files(directory)` - Recursively find all .cram files

**Usage:**
```python
from grid.utils.helper_dir.find_all_cram_files import find_all_cram_files

cram_files = find_all_cram_files("data/crams/")
# Returns: List of Path objects
```

---

### `load_flags_from_yaml.py`
**Purpose:** Load SAM flag filters from YAML configuration  
**Functions:**
- `load_flags_from_yaml(config_path)` - Parse SAM flags from config

**Usage:**
```python
from grid.utils.helper_dir.load_flags_from_yaml import load_flags_from_yaml

flags = load_flags_from_yaml("config.yaml")
# Returns: Dict with 'required_flags' and 'excluded_flags'
```

**Config Format:**
```yaml
sam_flags:
  required_flags:
    - PROPER_PAIR
    - READ_MAPPED
  excluded_flags:
    - DUPLICATE
    - SECONDARY
```

---

### `setup_output_file.py`
**Purpose:** Prepare output files and directories  
**Functions:**
- `setup_output_file(output_path, create_dirs=True)` - Ensure output path exists

**Usage:**
```python
from grid.utils.helper_dir.setup_output_file import setup_output_file

output_file = setup_output_file("results/analysis/output.tsv")
# Creates 'results/analysis/' directory if needed
```

---

### `write_result_to_file.py`
**Purpose:** Write analysis results to disk  
**Functions:**
- `write_result_to_file(results, output_file, header)` - Write TSV output
- `append_result_to_file(result, output_file)` - Append single result

**Usage:**
```python
from grid.utils.helper_dir.write_result_to_file import write_result_to_file

results = [
    {"sample": "S001", "reads": 1234, "coverage": 45.2},
    {"sample": "S002", "reads": 2345, "coverage": 52.1}
]
write_result_to_file(
    results, 
    "output.tsv", 
    header=["sample", "reads", "coverage"]
)
```

---

## Design Principles

These helpers follow DRY (Don't Repeat Yourself) principles:
- **Reusable:** Used across multiple pipeline modules
- **Focused:** Each module does one thing well
- **Flexible:** Accept various input formats when reasonable
- **Type-safe:** Use pathlib.Path for file operations

## Common Patterns

Most utilities follow these conventions:
- Accept `Path` or `str` for file paths
- Return standardized types (`Path`, `dict`, `list`)
- Raise descriptive exceptions on errors
- Support both programmatic and CLI usage