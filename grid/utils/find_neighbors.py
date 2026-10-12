# In[0]: Imports
from pathlib import Path
import subprocess

from .utils import log


# In[1]: Main Function
def find_neighbors(config, console=None):
    """
    Find nearest neighbors by calling the precompiled C++ executable,
    which matches Mukamel's pipeline exactly.

    C++ argument order:
      arg1: batch number (always 0 — single pass)
      arg2: total batches (always 1 — single pass)
      arg3: zMax
      arg4: ID_scale_zdepths input file
      arg5: output prefix (C++ appends .zMax<val>.txt.gz)

    Args:
        config (dict): Configuration data loaded from the config file.
        console: Rich console object for logging.
    """
    try:
        output_dir = Path(config.get("output_dir", "."))
        mosdepth_cfg = config.get("mosdepth", {})
        neighbors_cfg = mosdepth_cfg.get("neighbors", {})
        normalize_cfg = mosdepth_cfg.get("normalize", {})

        # Input: zdepths file produced by normalize_mosdepth
        input_prefix = normalize_cfg.get("output_file_prefix", "ID_scale_zdepths")
        input_file = output_dir / f"{input_prefix}.txt.gz"

        # Parameters
        zmax = float(neighbors_cfg.get("zmax", 2.0))

        # Output prefix — C++ will append .zMax<val>.txt.gz automatically
        output_prefix = neighbors_cfg.get("output_file_prefix", "ID_scale_neighbors")
        output_prefix_path = output_dir / output_prefix

        # Expected output file the C++ will produce
        output_file = output_dir / f"{output_prefix}.zMax{zmax:g}.txt.gz"

        # Path to compiled executable
        exe = Path(__file__).parent / "find_neighbors"

        # Optional parameters for the C++ executable
        N_OUTPUT = int(neighbors_cfg.get("num_neighbors", 500))
        sigma2max = float(neighbors_cfg.get("sigma2_max", 1000))
        fracR = float(neighbors_cfg.get("frac_r", 1.0))

    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    # Sanity checks
    if not exe.exists():
        log(
            console,
            f"Executable not found: {exe}\n"
            "Please compile find_neighbors.cpp and place the binary at "
            f"{exe}",
            style="danger",
        )
        return

    if not input_file.exists():
        log(console, f"Input zdepths file not found: {input_file}", style="danger")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # Build command — always batch 0 of 1 (single pass over all samples)
    cmd = [
        str(exe),
        "0",  # batch number
        "1",  # total batches
        str(zmax),  # zMax
        str(input_file),
        str(output_prefix_path),
        str(N_OUTPUT),
        str(sigma2max),
        str(fracR),
    ]

    log(console, "Running find_neighbors executable...", style="info")
    log(console, f"  Input zdepths : {input_file}", style="info")
    log(console, f"  zMax          : {zmax}", style="info")
    log(console, f"  Output prefix : {output_prefix_path}", style="info")
    log(console, f"  Expected out  : {output_file}", style="info")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        # Stream C++ stdout so user can see progress dots
        for line in result.stdout.splitlines():
            log(console, line, style="info")

    except subprocess.CalledProcessError as e:
        log(
            console,
            f"find_neighbors executable failed (exit {e.returncode}):",
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

    if output_file.exists():
        log(console, f"find_neighbors complete. Output: {output_file}", style="success")
    else:
        log(
            console,
            f"Executable finished but output file not found: {output_file}\n"
            f"C++ appends .zMax{zmax:g}.txt.gz to the output prefix — check the prefix is correct.",
            style="danger",
        )
