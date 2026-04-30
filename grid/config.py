# config.py
from pathlib import Path
from .utils.utils import log

# --- Top-level required fields ---
REQUIRED_TOP_LEVEL = {
    "samples_file": str,
    "directory_loc": str,
    "reference_genome": str,
    "output_dir": str,
    "threads": int,
    "file_type": str,
    "chrom": str,
    "start_bp": int,
    "end_bp": int,
    "output_file_type": str,
}

REQUIRED_FILES_TOP_LEVEL = ["samples_file", "reference_genome"]

STEP_SCHEMA = [
    # index
    {"path": ("index", "output_file_prefix"), "default": "'output'"},
    # count_reads
    {"path": ("count_reads", "min_mapq"), "gate": ("count_reads",), "default": "1"},
    {
        "path": ("count_reads", "output_file_prefix"),
        "gate": ("count_reads",),
        "default": "'output'",
    },
    {"path": ("count_reads", "flags"), "gate": ("count_reads",), "required": True},
    # mosdepth
    {"path": ("mosdepth", "output_file_prefix"), "gate": ("mosdepth",), "default": "'output'"},
    {"path": ("mosdepth", "bin_size"), "gate": ("mosdepth",), "default": "1000"},
    {"path": ("mosdepth", "mode"), "gate": ("mosdepth",), "default": "'fast'"},
    {
        "path": ("mosdepth", "work_dir"),
        "gate": ("mosdepth",),
        "default": "output_dir/mosdepth_workdir",
    },
    {"path": ("mosdepth", "remove_intermediate"), "gate": ("mosdepth",), "default": "True"},
    # mosdepth.normalize
    {
        "path": ("mosdepth", "normalize", "min_depth"),
        "gate": ("mosdepth", "normalize"),
        "default": "20",
    },
    {
        "path": ("mosdepth", "normalize", "max_depth"),
        "gate": ("mosdepth", "normalize"),
        "default": "100",
    },
    {
        "path": ("mosdepth", "normalize", "top_frac"),
        "gate": ("mosdepth", "normalize"),
        "default": "0.1",
    },
    {
        "path": ("mosdepth", "normalize", "output_file_prefix"),
        "gate": ("mosdepth", "normalize"),
        "default": "'output'",
    },
    {
        "path": ("mosdepth", "normalize", "repeat_mask_file"),
        "gate": ("mosdepth", "normalize"),
        "required": True,
        "is_file": True,
    },
    # mosdepth.neighbors
    {
        "path": ("mosdepth", "neighbors", "output_file_prefix"),
        "gate": ("mosdepth", "neighbors"),
        "default": "'output'",
    },
    {
        "path": ("mosdepth", "neighbors", "num_neighbors"),
        "gate": ("mosdepth", "neighbors"),
        "default": "5",
    },
    {
        "path": ("mosdepth", "neighbors", "zmax"),
        "gate": ("mosdepth", "neighbors"),
        "default": "2.0",
    },
    {
        "path": ("mosdepth", "neighbors", "sigma2_max"),
        "gate": ("mosdepth", "neighbors"),
        "default": "1000",
    },
    # compute_diploid_genotypes
    {
        "path": ("compute_diploid_genotypes", "output_file_prefix"),
        "gate": ("compute_diploid_genotypes",),
        "default": "'output'",
    },
    # compute_haploid_genotypes
    {
        "path": ("compute_haploid_genotypes", "phased_vcf"),
        "gate": ("compute_haploid_genotypes",),
        "required": True,
        "is_file": True,
    },
    {
        "path": ("compute_haploid_genotypes", "ibs_output"),
        "gate": ("compute_haploid_genotypes",),
        "required": True,
    },
    {
        "path": ("compute_haploid_genotypes", "output_file_prefix"),
        "gate": ("compute_haploid_genotypes",),
        "default": "'output'",
    },
    {
        "path": ("compute_haploid_genotypes", "min_neighbors"),
        "gate": ("compute_haploid_genotypes",),
        "default": "1",
    },
    {
        "path": ("compute_haploid_genotypes", "max_neighbors"),
        "gate": ("compute_haploid_genotypes",),
        "default": "10",
    },
    {
        "path": ("compute_haploid_genotypes", "n_iters"),
        "gate": ("compute_haploid_genotypes",),
        "default": "100",
    },
]


def _get_nested(config, *keys):
    """Navigate a nested dict by key path, returning None if any key is missing."""
    node = config
    for key in keys:
        if not isinstance(node, dict):
            return None
        node = node.get(key)
    return node


def _is_enabled(config, gate):
    """Return True if the config section at gate path has run=True."""
    section = _get_nested(config, *gate)
    return isinstance(section, dict) and section.get("run") is True


def validate_top_level(config, errors, warnings):
    for key, expected_type in REQUIRED_TOP_LEVEL.items():
        if key not in config:
            errors.append(f"Missing required field: '{key}'")
        elif not isinstance(config[key], expected_type):
            errors.append(f"'{key}' must be {expected_type.__name__}")

    for key in REQUIRED_FILES_TOP_LEVEL:
        val = config.get(key)
        if val and not Path(val).exists():
            errors.append(f"File not found: {key} = {val}")


def validate_steps(config, errors, warnings):
    for entry in STEP_SCHEMA:
        gate = entry.get("gate")
        if gate and not _is_enabled(config, gate):
            continue

        value = _get_nested(config, *entry["path"])
        field_name = ".".join(entry["path"])

        if value is None:
            if entry.get("required"):
                errors.append(f"{field_name} not set.")
            else:
                warnings.append(f"{field_name} not set. Defaulting to {entry['default']}.")
        elif entry.get("is_file") and not Path(value).exists():
            errors.append(f"File not found: {field_name} = {value}")


def error_check_config(config, console):
    errors = []
    warnings = []

    validate_top_level(config, errors, warnings)
    validate_steps(config, errors, warnings)

    if errors:
        for e in errors:
            log(console, e, style="danger")
        raise ValueError(f"{len(errors)} config error(s) found. Aborting.")

    if warnings:
        for w in warnings:
            log(console, w, style="warning")
        log(
            console,
            f"{len(warnings)} config warning(s) found. Please review. This may affect the results.",
            style="warning",
        )
