# grid/utils/normalize_mosdepth_dir/process_one_individual.py
# In[1]: Imports
import gzip
from pathlib import Path
from .get_regions_bed import regions_bed_gz

# In[2]: Worker function to process one individual
def _process_one_individual(args_tuple):
    """
    This runs in its own process. Returns (individual_id, [(start,end,depth),...])
    args_tuple contains:
    
    individual_id, global_dist_file, chromosome, start, end, min_depth, max_depth, excluded
    """
    (individual_id, global_dist_file, chromosome, start, end, min_depth, max_depth, excluded) = args_tuple
    bed_gz = regions_bed_gz(global_dist_file)
    results = []
    if not bed_gz.exists():
        return individual_id, results

    # Normalize chromosome name - handle both 'chr6' and '6' formats
    # The BED file likely has 'chr6' format, so ensure we're checking correctly
    chrom_to_match = chromosome if chromosome.startswith("chr") else f"chr{chromosome}"

    try:
        with gzip.open(bed_gz, "rt") as f:
            for line in f:
                # Check if line starts with the target chromosome
                if not line.startswith(chrom_to_match):
                    continue
                    
                fields = line.strip().split("\t")
                if len(fields) < 4:
                    continue
                
                chrom = fields[0]
                region_start = int(fields[1])
                region_end = int(fields[2])
                depth = float(fields[3])

                # Keep if depth>0 AND overlaps target window
                if not (depth > 0 and region_end >= start and region_start <= end):
                    continue

                # depth filter
                if depth < min_depth or depth > max_depth:
                    continue

                # repeat exclusion (kb bins)
                if excluded is not None:
                    region_kb = set(range(region_start // 1000, region_end // 1000 + 1))
                    if region_kb & excluded.get(chrom, set()):
                        continue

                results.append((region_start, region_end, depth))
    except Exception:
        # If file corrupt or other IO error, just return empty
        return individual_id, []

    return individual_id, results