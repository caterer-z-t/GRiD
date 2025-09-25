# In[0]: Imports
import os
import subprocess
from pathlib import Path
from GRiD.utils.utils import format_region, arg_check, argument_parser

# In[1]: Functions
def run_mosdepth(cram_file: Path, out_prefix: Path, reference: Path, window: int = 1000, n_threads: int = None):
    """Run mosdepth on a single CRAM file."""
    n_threads = n_threads or (os.cpu_count() or 1)
    cmd = [
        "mosdepth",
        "-n",
        "--fast-mode",
        "--by", str(window),
        "-f", str(reference),
        "-t", str(n_threads),
        str(out_prefix),
        str(cram_file)
    ]
    subprocess.run(cmd, check=True)

# In[2]: Main
def main():
    parser = argument_parser()
    args = parser.parse_args()
    arg_check(args)

    if args.region is None:
        args.region = format_region(args)

    # Ensure mosdepth is installed/accessible
    try:
        subprocess.run(["mosdepth", "--version"], check=True, capture_output=True)
    except FileNotFoundError:
        raise EnvironmentError("mosdepth is not installed or not found in PATH.")
    
    # Run mosdepth
    run_mosdepth(Path(args.input), Path(args.output).with_suffix(""), Path(args.reference), n_threads=args.threads)


# In[3]: Run
if __name__ == "__main__":
    main()
