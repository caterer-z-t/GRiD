# In[0]: Imports
from pathlib import Path
import subprocess

from .utils import log


# In[1]: Main Function
def normalize_mosdepth(config, console=None):
    """
    Normalize mosdepth coverage by calling the precompiled C++ executable,
    which matches Mukamel's pipeline exactly.

    Args:
        config (dict): Configuration data loaded from the config file.
        console: Rich console object for logging.
    """
    try:
        output_dir = Path(config.get("output_dir", "."))
        mosdepth_cfg = config.get("mosdepth", {})
        normalize_cfg = mosdepth_cfg.get("normalize", {})

        # Combined bins file — default to mosdepth output_file_prefix in output_dir
        combined_bins_gz_cfg = normalize_cfg.get("combined_bins_gz", None)
        if combined_bins_gz_cfg is not None:
            combined_bins_gz = Path(combined_bins_gz_cfg)
        else:
            # Fall back to mosdepth output file
            mosdepth_prefix = mosdepth_cfg.get("output_file_prefix", "mosdepth_results")
            output_file_type = config.get("output_file_type", "tsv")
            combined_bins_gz = output_dir / f"{mosdepth_prefix}.{output_file_type}.gz"

        # Repeat mask BED file
        repeat_mask = Path(normalize_cfg["repeat_mask_file"])

        # Example regions bed
        mosdepth_dir = Path(mosdepth_cfg.get("work_dir", output_dir / "mosdepth"))
        example_bed_cfg = mosdepth_cfg.get("example_bed", None)
        if example_bed_cfg is not None:
            example_bed = Path(example_bed_cfg)
        else:
            example_bed = mosdepth_dir / "example.regions.bed.gz"

        # Output file
        output_file_prefix = normalize_cfg.get("output_file_prefix", "ID_scale_zdepths")
        output_gz = output_dir / f"{output_file_prefix}.txt.gz"

        # Number of samples
        samples_file = config["samples_file"]
        n_samples = sum(1 for line in open(samples_file) if line.strip())

        # optional min/max depth parameters
        min_depth = normalize_cfg.get("min_depth", 20.0)
        max_depth = normalize_cfg.get("max_depth", 100.0)

        # Path to the compiled executable — lives in grid/utils/
        exe = Path(__file__).parent / "normalize_mosdepth"

    except KeyError as e:
        log(console, f"Config error — missing required key: {e}", style="danger")
        return
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    # Sanity checks
    if not exe.exists():
        log(
            console,
            f"Executable not found: {exe}\n"
            "Please compile normalize_mosdepth.cpp and place the binary at "
            f"{exe}",
            style="danger",
        )
        return

    if not combined_bins_gz.exists():
        log(console, f"Combined bins file not found: {combined_bins_gz}", style="danger")
        return

    if not repeat_mask.exists():
        log(console, f"Repeat mask file not found: {repeat_mask}", style="danger")
        return

    if not example_bed.exists():
        log(console, f"Example regions bed not found: {example_bed}", style="danger")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # Build command — matches C++ argument order exactly:
    # arg1: combined_bins_gz
    # arg2: repeat_mask_bed
    # arg3: example_regions_bed_gz
    # arg4: n_samples
    # arg5: output_gz
    cmd = [
        str(exe),
        str(combined_bins_gz),
        str(repeat_mask),
        str(example_bed),
        str(n_samples),
        str(output_gz),
        str(min_depth),
        str(max_depth),
    ]

    log(console, "Running normalize_mosdepth executable...", style="info")
    log(console, f"  Combined bins : {combined_bins_gz}", style="info")
    log(console, f"  Repeat mask   : {repeat_mask}", style="info")
    log(console, f"  Example bed   : {example_bed}", style="info")
    log(console, f"  N samples     : {n_samples}", style="info")
    log(console, f"  Output        : {output_gz}", style="info")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        # Stream C++ stdout so user can see normalization progress
        for line in result.stdout.splitlines():
            log(console, line, style="info")

    except subprocess.CalledProcessError as e:
        log(
            console,
            f"normalize_mosdepth executable failed (exit {e.returncode}):",
            style="danger",
        )
        for line in (e.stdout or "").splitlines():
            log(console, line, style="danger")
        return
    except FileNotFoundError:
        log(
            console,
            f"Could not execute {exe} — check permissions (chmod +x).",
            style="danger",
        )
        return

    if output_gz.exists():
        log(console, f"Normalization complete. Output: {output_gz}", style="success")
    else:
        log(
            console,
            f"Executable finished but output file not found: {output_gz}",
            style="danger",
        )
