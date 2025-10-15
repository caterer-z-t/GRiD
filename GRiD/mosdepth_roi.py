# In[0]: Imports
import gzip
from pathlib import Path
from GRiD.utils.utils import format_region, arg_check, argument_parser

# In[1]: Functions

def _strip_chr(x: str) -> str:
    return x[3:] if x.startswith("chr") else x

"""
def extract_region_cov(regions_bed_gz: Path, region: str) -> int:
    '''
    Extract coverage for a given region from a mosdepth .regions.bed.gz file.

    Parameters
    ----------
    regions_bed_gz : Path
        Path to mosdepth .regions.bed.gz file.
    region : str
        Genomic region string ("chr:start-end").

    Returns
    -------
    int
        Coverage scaled by 100 and rounded (like int(100*$4+0.5) in bash).
    '''
    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))
    
    # extract number from chromosome if prefixed with "chr"
    if chrom.startswith("chr"):
        chrom = chrom[3:]

    with gzip.open(regions_bed_gz, "rt") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            c, s, e, cov = fields[:4]
            if c == chrom and int(s) == start and int(e) == end:
                return int(100 * float(cov) + 0.5)

    raise ValueError(f"Region {chrom}:{start}-{end} not found in {regions_bed_gz}")
"""
def extract_region_cov(
    regions_bed_gz: Path,
    region: str,
    verbose: bool = False,
    compute_fallback: bool = False,
    max_show: int = 50
) -> int:
    """
    Debug-friendly extractor that:
      - tries strict match first (like your zgrep+awk),
      - if not found, prints helpful lines to diagnose,
      - optionally computes a weighted average coverage across overlapping bins.

    Parameters
    ----------
    regions_bed_gz : Path
        Path to mosdepth .regions.bed.gz file.
    region : str
        Genomic region string ("chr:start-end" or "start-end" with/without "chr").
    verbose : bool
        Print debugging info (matching lines, sample head).
    compute_fallback : bool
        If True and exact match isn't found, compute weighted average over overlapping bins.
    max_show : int
        Max number of matching lines to print.
    """
    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))

    target_norm = _strip_chr(chrom)

    if verbose:
        print(f"[debug] looking for region: {chrom}:{start}-{end} (norm '{target_norm}')", file=sys.stderr)
        print(f"[debug] reading file: {regions_bed_gz}", file=sys.stderr)

    matches = []  # will store tuples (s,e,cov,line)
    strict_found = None

    with gzip.open(regions_bed_gz, "rt") as fh:
        for lineno, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            c, s_text, e_text, cov_text = fields[:4]
            try:
                s = int(s_text); e = int(e_text); cov = float(cov_text)
            except ValueError:
                # skip lines that aren't numeric where expected
                continue

            # normalize chromosome for comparison
            if _strip_chr(c) == target_norm:
                # collect matching-chromosome lines
                if len(matches) < max_show:
                    matches.append((s, e, cov, line.rstrip("\n")))
                else:
                    # keep counting but don't store more than max_show
                    matches.append(None)

                # check strict exact-match
                if s == start and e == end:
                    strict_found = cov
                    if verbose:
                        print(f"[debug] exact match at line {lineno}: {c}\t{s}\t{e}\t{cov}", file=sys.stderr)
                    return int(100 * strict_found + 0.5)

    # if we reach here, strict match not found
    if verbose:
        print(f"[debug] exact region not found.", file=sys.stderr)
        if matches:
            print(f"[debug] {len(matches)} line(s) found for chromosome '{chrom}' (showing up to {max_show}):", file=sys.stderr)
            for m in matches[:max_show]:
                if m is None:
                    break
                s, e, cov, line = m
                overlap = max(0, min(end, e) - max(start, s))
                print(f"[debug]   {line}\t(overlap={overlap})", file=sys.stderr)
        else:
            print(f"[debug] no lines for chromosome '{chrom}' found in the file.", file=sys.stderr)

        # show small head of the file to confirm format (first 20 lines)
        try:
            with gzip.open(regions_bed_gz, "rt") as fh2:
                print("[debug] file head (first 20 non-empty lines):", file=sys.stderr)
                cnt = 0
                for l in fh2:
                    if not l.strip():
                        continue
                    print(f"[debug]   {l.rstrip()}", file=sys.stderr)
                    cnt += 1
                    if cnt >= 20:
                        break
        except Exception as ex:
            print(f"[debug] error reading file head: {ex}", file=sys.stderr)

    # optional fallback: compute weighted average coverage over overlapping windows
    if compute_fallback and matches:
        total_cov_times_len = 0.0
        total_len = 0
        # Re-open and iterate all lines to be exact (we only stored a sample above)
        with gzip.open(regions_bed_gz, "rt") as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                c, s_text, e_text, cov_text = fields[:4]
                if _strip_chr(c) != target_norm:
                    continue
                try:
                    s = int(s_text); e = int(e_text); cov = float(cov_text)
                except ValueError:
                    continue
                # overlap check
                if e <= start or s >= end:
                    continue
                overlap = min(end, e) - max(start, s)
                if overlap > 0:
                    total_cov_times_len += cov * overlap
                    total_len += overlap

        if total_len > 0:
            avg_cov = total_cov_times_len / total_len
            cov_int = int(round(avg_cov * 100))
            if verbose:
                print(f"[debug] fallback weighted avg cov over overlapping windows: {avg_cov} -> {cov_int}", file=sys.stderr)
            return cov_int
        else:
            if verbose:
                print("[debug] fallback found no overlapping windows either.", file=sys.stderr)

    # if still not found, raise with detailed message
    raise ValueError(
        f"Region {region} not found in {regions_bed_gz}.\n"
        f"  - Checked chrom norm: '{target_norm}'.\n"
        f"  - If mosdepth was run in fixed windows, the exact start/end may not appear.\n"
        f"  - Try running mosdepth with a BED containing the target interval or use compute_fallback=True."
    )

# In[2]: Main
def main():
    parser = argument_parser()
    args = parser.parse_args()
    arg_check(args)

    if args.region is None:
        args.region = format_region(args)

    # regions.bed.gz is expected to share prefix with mosdepth output
    regions_bed_gz = Path(args.output).with_suffix(".regions.bed.gz")

    if not regions_bed_gz.exists():
        raise FileNotFoundError(f"{regions_bed_gz} does not exist. Run mosdepth first.")

    sample_id = Path(args.input).stem
    coverage = extract_region_cov(regions_bed_gz, args.region)

    with open(args.output, "a") as out:
        out.write(f"{sample_id}\t{coverage}\n")


# In[3]: Run
if __name__ == "__main__":
    main()