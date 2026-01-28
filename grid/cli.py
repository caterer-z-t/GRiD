"""Command-line interface for LPA KIV-CNV pipeline."""
# In[0]: Imports
import sys
import click
from pathlib import Path
from rich.console import Console

from . import __version__

# In[1]: CLI Setup
console = Console()

def print_banner():
    """Print ASCII banner."""
    banner = """
    ╔═══════════════════════════════════════════════════════════╗
    ║                                                           ║
    ║      GRiD - Genomic Repeat inference from Depth           ║
    ║           LPA KIV-2 Copy Number Variant Pipeline          ║
    ║                      Version {version}                        ║
    ║                                                           ║
    ╚═══════════════════════════════════════════════════════════╝
    """.format(version=__version__)
    console.print(banner, style="bold red")


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-v', '--version', is_flag=True, is_eager=True, expose_value=False,
              help="Show the GRiD version", callback=lambda ctx, param, value: print_version_and_exit(value))

def cli():
    """
    GRiD - Genomic Repeat inference from Depth
    
    A modular pipeline for estimating copy number variants in the LPA gene's KIV-2 VNTR region.
    """
    pass

def print_version_and_exit(value):
    if value:
        console.print(f"GRiD version: {__version__}")
        raise SystemExit(0)

# In[2]: CLI Command for Ensuring CRAI Index
@cli.command()
@click.option("-c", "--cram", required=True, type=click.Path(exists=True), help="Input CRAM file")
@click.option("-r", "--reference", required=True, type=click.Path(exists=True), help="Reference genome FASTA")
def crai(cram, reference):
    """
    Ensure CRAI index exists for a CRAM file.
    
    Args:
        cram: Path to the CRAM file.
        reference: Path to the reference genome FASTA.

    Returns:
        Path to the CRAI index file.
    """
    from .utils.ensure_crai import ensure_crai

    print_banner()
    try:
        with console.status("[bold green]Checking/creating CRAI...", spinner="dots"):
            crai_path = ensure_crai(cram_path=cram, reference=reference)
        console.print(f"[green]✓ CRAI index ready:[/green] {crai_path}")
    except Exception as e:
        console.print(f"[red]✗ Failed to create CRAI: {str(e)}[/red]")
        sys.exit(1)

# In[3]: CLI Command for Batch Indexing CRAMs
@cli.command()
@click.option("-C", "--cram-dir", required=True, type=click.Path(exists=True), help="Directory with CRAM files")
@click.option("-r", "--reference", required=True, type=click.Path(exists=True), help="Reference genome FASTA")
@click.option("-t", "--threads", default=1, type=int, help="Number of threads for parallel processing")

def batch_crai(cram_dir, reference, threads):
    """
    Ensure CRAI indices exist for all CRAM files in a directory.
    
    Args:
        cram_dir: Directory containing CRAM files.
        reference: Path to the reference genome FASTA.
        threads: Number of threads for parallel processing.

    Returns:
        None
    """
    from .utils.batch_crai import batch_crai

    print_banner()
    try:
        with console.status("[bold green]Batch checking/creating CRAIs...", spinner="dots"):
            batch_crai(cram_dir=cram_dir, reference=reference, console=console, threads=threads)
        console.print(f"[green]✓ All CRAI indices are ready.[/green]")
    except Exception as e:
        console.print(f"[red]✗ Batch CRAI creation failed: {str(e)}[/red]")
        sys.exit(1)
        
# In[4]: CLI Command for Subsetting CRAM
@cli.command()
@click.option('-c', '--cram', required=True, type=click.Path(exists=True), help='Input CRAM file')
@click.option('-r', '--region', required=False, help="Genomic region (e.g., 'chr6:160000000-160100000')")
@click.option('-C', '--chrom', required=False, help='Chromosome (e.g., chr6)')
@click.option('-s', '--start', required=False, type=int, help='Start position (e.g., 160000000)')
@click.option('-e', '--end', required=False, type=int, help='End position (e.g., 160100000)')
@click.option('-o', '--output', required=True, type=click.Path(), help='Output subset CRAM file')
@click.option('-R', '--reference', required=True, type=click.Path(exists=True), help='Reference genome FASTA')
def subset(cram, region, output, reference):
    """
    Subset a CRAM file to a specific genomic region.
    
    Args:
        cram: Input CRAM file.
        region: Genomic region in 'chr:start-end' format.
        chrom: Chromosome (if region not provided).
        start: Start position (if region not provided).
        end: End position (if region not provided).
        output: Output CRAM file path.
        reference: Reference genome FASTA.

    Returns:
        Path to the subset CRAM file.
    """
    from .utils.subset_cram import subset_cram
    from .helper.create_region import create_region_string
    
    region = create_region_string(region, chrom, start, end)

    print_banner()
    try:
        with console.status(f"[bold green]Subsetting CRAM for region {region}...", spinner="dots"):
            subset_path = subset_cram(cram_path=cram, region=region, output_path=output, reference=reference)
        console.print(f"[green]✓ Subset CRAM created:[/green] {subset_path}")
    except Exception as e:
        console.print(f"[red]✗ Failed to subset CRAM: {str(e)}[/red]")
        sys.exit(1)

# In[]: CLI Command for Batch Subsetting CRAMs
@cli.command()
@click.option("-C", "--cram-dir", required=True, type=click.Path(exists=True), help="Directory with CRAM files")
@click.option("-r", "--region", required=False, help="Genomic region (e.g., 'chr6:160000000-160100000')")
@click.option("-c", "--chrom", required=False, help="Chromosome (e.g., chr6)")
@click.option("-s", "--start", required=False, type=int, help="Start position (e.g., 160000000)")
@click.option("-e", "--end", required=False, type=int, help="End position (e.g., 160100000)")
@click.option("-o", "--output-dir", required=True, type=click.Path(), help="Output directory for subset CRAM files")
@click.option("-R", "--reference", required=True, type=click.Path(exists=True), help="Reference genome FASTA")
def batch_subset(cram_dir, region, chrom, start, end, output_dir, reference):
    """
    Subset all CRAM files in a directory to a specific genomic region.
    
    Args:
        cram_dir: Directory containing CRAM files.
        region: Genomic region in 'chr:start-end' format.
        chrom: Chromosome (if region not provided).
        start: Start position (if region not provided).
        end: End position (if region not provided).
        output_dir: Directory to save subset CRAM files.
        reference: Reference genome FASTA.
    Returns:
        None
    """
    from .utils.batch_subset_cram import batch_subset_cram
    from .helper.create_region import create_region_string
    
    region = create_region_string(region, chrom, start, end)

    print_banner()
    try:
        with console.status(f"[bold green]Batch subsetting CRAMs for region {region}...", spinner="dots"):
            batch_subset_cram(
                cram_dir=cram_dir,
                region=region,
                output_dir=output_dir,
                reference=reference,
                console=console
            )
        console.print(f"[green]✓ All subset CRAM files created in {output_dir}.[/green]")
    except Exception as e:
        console.print(f"[red]✗ Batch CRAM subsetting failed: {str(e)}[/red]")
        sys.exit(1)

# In[5]: CLI Command for Google Cloud Copy
@cli.command()
@click.option(
    "-b", "--bucket-path",
    required=True,
    type=str,
    help="Google Cloud Storage bucket path (e.g., gs://your-bucket/path/)")
@click.option("-C", "--cram-dir", default="~/cram", type=click.Path(), help="Local directory to copy CRAM files into")
def gcloud_copy(bucket_path, cram_dir):
    """Copy CRAM files from a Google Cloud bucket to a local directory."""
    console.print(f"[bold red]This command is deprecated. We currently expect the user to have access to the CRAM files before continuing, please ensure you have said files.[/bold red]")
    exit(1)
    from .utils.google_cloud_copy import copy_crams_from_bucket
    print_banner()
    copy_crams_from_bucket(bucket_path, cram_dir)

# In[6]: CLI Command for Counting Reads in Directory
@cli.command()
@click.option("-C", "--cram-dir", required=True, type=click.Path(exists=True), help="Directory with CRAM files")
@click.option("-o", "--output-file", required=True, type=click.Path(), help="Output TSV file")
@click.option("-r", "--ref-fasta", required=True, type=click.Path(exists=True), help="Reference genome FASTA")
@click.option("-c", "--chrom", required=True, help="Chromosome (e.g., chr6)")
@click.option("-s", "--start", required=True, type=int, help="Start position")
@click.option("-e", "--end", required=True, type=int, help="End position")
@click.option("--config", required=True, type=click.Path(exists=True), help="Path to YAML config file")
@click.option("-t", "--threads", default=1, type=int, help="Number of threads for parallel processing (default: 1)")
def count_reads(cram_dir, output_file, ref_fasta, chrom, start, end, config, threads):
    """Count properly paired reads in a specified LPA VNTR region for all CRAMs."""
    console.print(f"[bold green]Counting reads in CRAMs...[/bold green]")
    print_banner()

    from .utils.count_reads import count_reads
    
    try:
        count_reads(cram_dir, output_file, ref_fasta, chrom, start, end, config, threads)
    except Exception as e:
        console.print(f"[red]✗ Read counting failed: {str(e)}[/red]")
        sys.exit(1)

# In[7]: CLI Command for Mosdepth Coverage
@cli.command()
@click.option("-C", "--cram-dir", required=True, type=click.Path(exists=True), help="Directory with CRAM files")
@click.option("-o", "--output-file", required=True, type=click.Path(), help="Output TSV file")
@click.option("-r", "--ref-fasta", required=True, type=click.Path(exists=True), help="Reference genome FASTA")
@click.option("-c", "--chrom", required=True, help="Chromosome (e.g., chr6)")
@click.option("-s", "--start", required=True, type=int, help="Start position")
@click.option("-e", "--end", required=True, type=int, help="End position")
@click.option('-n', '--region-name', required=True, type=str, help="Name of the region for mosdepth output (default: LPA_VNTR)")
@click.option("-w", "--work-dir", required=False, default="~/mosdepth_work", type=click.Path(), help="Working directory for intermediate files")
@click.option("--by", required=False, default=1000, type=int, help="Bin size for mosdepth --by (Default: 1000)")
@click.option("--fast/--no-fast", required=False, default=True, help="Enable or disable mosdepth --fast-mode")
@click.option("-t", "--threads", required=False, default=1, type=int, help="Number of threads for parallel processing (default: 1)")
def mosdepth(
    cram_dir, output_file, ref_fasta, chrom, start, end, region_name, work_dir, by, fast, threads
):
    """Run mosdepth on all CRAMs in a directory and extract coverage for the LPA VNTR region."""
    print_banner()
    
    from .utils.run_mosdepth import run_mosdepth 
    
    try:
        run_mosdepth(
            cram_dir=cram_dir,
            output_file=output_file,
            ref_fasta=ref_fasta,
            chrom=chrom,
            start=start,
            end=end,
            region_name=region_name,
            work_dir=work_dir,
            by=by,
            fast_mode=fast,
            threads=threads
        )
    except Exception as e:
        console.print(f"[red]✗ Mosdepth coverage failed: {str(e)}[/red]")
        sys.exit(1)

# In[8]: CLI Command for Normalize Mosdepth Coverage
@cli.command()
@click.option("-M", "--mosdepth-dir", required=True, type=click.Path(exists=True), help="Directory with mosdepth files")
@click.option("-o", "--output-file", required=True, type=click.Path(), help="Normalized output file (.gz)")
@click.option("-r", "--repeat-mask", required=True, type=click.Path(exists=True), help="Repeat mask BED file")
@click.option("-c", "--chrom", required=True, help="Chromosome (e.g., chr6)")
@click.option("-s", "--start", required=True, type=int, help="Start coordinate")
@click.option("-e", "--end", required=True, type=int, help="End coordinate")
@click.option("--min-depth", required=False, default=20, type=int, help="Minimum depth to keep (default: 20)")
@click.option("--max-depth", required=False, default=100, type=int, help="Maximum depth to keep (default: 100)")
@click.option("--top-frac", required=False, default=0.1, type=float, help="Top fraction of high-variance regions to keep (default: 0.1)")
@click.option("-t", "--threads", required=False, default=1, type=int, help="Number of threads for parallel processing (default: 1)")
def normalize_mosdepth(
    mosdepth_dir, output_file, repeat_mask, chrom, start, end, 
    min_depth, max_depth, top_frac, threads
):
    """Normalize mosdepth coverage across all samples in a directory."""
    print_banner()
    
    from .utils.run_normalized_mosdepth import run_normalize_mosdepth
    
    try:
        run_normalize_mosdepth(
            mosdepth_dir=mosdepth_dir,
            output_file=output_file,
            repeat_mask=repeat_mask,
            chrom=chrom,
            start=start,
            end=end,
            min_depth=min_depth,
            max_depth=max_depth,
            top_frac=top_frac,
            threads=threads
        )
    except Exception as e:
        console.print(f"[red]✗ Mosdepth normalization failed: {str(e)}[/red]")
        sys.exit(1)

# In[]: updated normalize mosdepth
@cli.command()
@click.option("-M", "--mosdepth-dir", required=True, type=click.Path(exists=True), help="Directory with mosdepth files")
@click.option("-r", "--repeat-mask", required=True, type=click.Path(exists=True), help="Repeat mask BED file")
@click.option("-c", "--chrom", required=True, help="Chromosome (e.g., chr6)")
@click.option("-s", "--start", required=True, type=int, help="Start coordinate")
@click.option("-e", "--end", required=True, type=int, help="End coordinate")
@click.option("-o", "--output-file", required=True, type=click.Path(), help="Normalized output file (.gz)")
@click.option("--min-depth", required=False, default=20, type=int, help="Minimum depth to keep (default: 20)")
@click.option("--max-depth", required=False, default=100, type=int, help="Maximum depth to keep (default: 100)")
@click.option("--top-frac", required=False, default=0.1, type=float, help="Top fraction of high-variance regions to keep (default: 0.1)")
@click.option("-t", "--threads", required=False, default=1, type=int, help="Number of threads for parallel processing (default: 1)")
def normalize_mosdepth_update(
    mosdepth_dir, repeat_mask, chrom, start, end, output_file,
    min_depth, max_depth, top_frac, threads
):
    """Normalize mosdepth coverage across all samples in a directory."""
    print_banner()
    
    from .utils.normalize_mosdepth_dir.normalize_mosdepth import run_normalize_mosdepth
    
    try:
        run_normalize_mosdepth(
            mosdepth_dir=mosdepth_dir,
            output_file=output_file,
            repeat_mask=repeat_mask,
            chrom=chrom,
            start=start,
            end=end,
            min_depth=min_depth,
            max_depth=max_depth,
            top_frac=top_frac,
            threads=threads
        )
    except Exception as e:
        console.print(f"[red]✗ Mosdepth normalization failed: {str(e)}[/red]")
        sys.exit(1)
        
# In[9]: Find Neighbors Command
@cli.command()
@click.option("-i", "--input-file", required=True, type=click.Path(exists=True),
              help="Input normalized z-depth file (.gz) from normalize_mosdepth")
@click.option("-o", "--output-file", required=True, type=click.Path(),
              help="Output file (will be gzipped and converted to .zMax{zmax}.txt.gz)")
@click.option("-z", "--zmax", default=2.0, show_default=True, type=float,
              help="Maximum z-score clipping threshold")
@click.option("-n", "--n-neighbors", default=500, show_default=True, type=int,
              help="Number of nearest neighbors to find")
@click.option("--sigma2-max", default=1000.0, show_default=True, type=float,
              help="Maximum variance threshold for filtering regions")
def find_neighbors(input_file, output_file, zmax, n_neighbors, sigma2_max):
    """Find nearest neighbors among individuals using normalized z-depths."""
    from .utils.find_neighbors import find_neighbors

    print_banner()
    try:
        find_neighbors(
            input_file=input_file,
            output_file=output_file,
            zmax=zmax,
            n_neighbors=n_neighbors,
            sigma2_max=sigma2_max
        )
    except Exception as e:
        console.print(f"[red]✗ Neighbor finding failed:[/red] {str(e)}")
        sys.exit(1)

# In[10]: Extract Reference Command
@cli.command()
@click.option("-r", "--reference-fa", required=True, type=click.Path(exists=True),
              help="Path to reference genome FASTA (e.g., hs37d5.fa)")
@click.option("-b", "--bed-file", required=True, type=click.Path(exists=True),
              help="BED file defining regions to extract")
@click.option("-o", "--output-dir", required=True, type=click.Path(),
              help="Output directory for FASTA file")
@click.option("-f", "--output-prefix", default="ref_lpa", show_default=True,
              help="Prefix for output FASTA file")
def extract_reference(reference_fa, bed_file, output_dir, output_prefix):
    """Extract FASTA sequences from a reference genome based on BED regions."""
    from .utils.extract_reference import extract_reference as extract_reference_fn

    print_banner()
    try:
        extract_reference_fn(reference_fa, bed_file, output_dir, output_prefix)
    except Exception as e:
        console.print(f"[red]✗ Reference extraction failed:[/red] {str(e)}")
        sys.exit(1)

# In[11]: CLI Command for Realignment
@cli.command()
@click.option("-C", "--cram-dir", required=True, type=click.Path(exists=True), help="Directory with CRAM files")
@click.option("-o", "--output-file", required=True, type=click.Path(), help="Output TSV file")
@click.option("-r", "--ref-fasta", required=True, type=click.Path(exists=True), help="Reference genome FASTA")
@click.option("-f", "--lpa-ref-fasta", required=True, type=click.Path(exists=True), help="LPA KIV-2 reference FASTA for realignment")
@click.option("-P", "--positions-file", required=True, type=click.Path(exists=True), help="Hardcoded positions file for LPA KIV-2")
@click.option("-g", "--genome-build", default="hg38", type=click.Choice(['hg19', 'hg37', 'hg38']), show_default=True, help="Genome build")
@click.option("-c", "--chrom", required=True, help="Chromosome (e.g., chr6)")
@click.option("-s", "--start", required=True, type=int, help="Start position (0-based)")
@click.option("-e", "--end", required=True, type=int, help="End position (0-based)")
@click.option("-t", "--threads", default=1, type=int, help="Number of threads for parallel processing (default: 1)")
def lpa_realign(cram_dir, output_file, ref_fasta, lpa_ref_fasta, positions_file, genome_build, chrom, start, end, threads):
    """Realign reads in LPA KIV-2 region for all CRAMs in a directory."""
    from .utils.align_lpa import run_lpa_realignments
    print_banner()

    try:
        run_lpa_realignments(
        cram_dir=cram_dir,
        reference_fa=ref_fasta,
        lpa_ref_fasta=lpa_ref_fasta,
        positions_file=positions_file,
        genome_build=genome_build,
        chrom=chrom,
        start=start,
        end=end,
        output_file=output_file,
        threads=threads
        )

    except Exception as e:
        console.print(f"[red]✗ LPA realignment failed:[/red] {str(e)}")
        sys.exit(1)

# In[12]: CLI Command for Computing Diploid Copy Numbers
@cli.command()
@click.option(
    "-c", "--count-file", required=True, type=click.Path(exists=True),
    help="Realignment output file with read counts"
)
@click.option(
    "-n", "--neighbor-file", required=True, type=click.Path(exists=True),
    help="Find neighbors output file (can be gzipped)"
)
@click.option(
    "-o", "--output-prefix", required=True, type=click.Path(),
    help="Output prefix for diploid CN files"
)
@click.option(
    "-N", "--n-neighbors", default=200, show_default=True, type=int,
    help="Number of top neighbors to use"
)
def compute_dipcn(count_file, neighbor_file, output_prefix, n_neighbors):
    """
    Compute diploid copy numbers for LPA KIV-2 repeats.

    Uses neighbor-based normalization to compute dipCN for each exon type:
    1B_KIV3, 1B_notKIV3, 1B, and 1A.
    """
    from .utils.compute_dipcn import compute_dipcn_pipeline
    print_banner()

    # Run the main Python computation
    try:
        compute_dipcn_pipeline(
            count_file=Path(count_file).expanduser(),
            neighbor_file=Path(neighbor_file).expanduser(),
            output_prefix=output_prefix,
            n_neighbors=n_neighbors,
            console=console
        )
    except Exception as e:
        console.print(f"[red]✗ Diploid CN computation failed: {e}[/red]")
        sys.exit(1)

    console.print(f"[green]Diploid copy number computation completed successfully![/green]")

# In[13]: Estimate KIV Copy Numbers Command
@cli.command()
@click.option(
    "-a", "--exon1a", required=True, type=click.Path(exists=True),
    help="Exon1A diploid CN file (*.exon1A.dipCN.txt)"
)
@click.option(
    "-b", "--exon1b", required=True, type=click.Path(exists=True),
    help="Exon1B diploid CN file (*.exon1B.dipCN.txt)"
)
@click.option(
    "-o", "--output", required=True, type=click.Path(),
    help="Output file path for KIV2 copy number estimates"
)
@click.option(
    "-f", "--format", type=click.Choice(['tsv', 'txt', 'csv']), default='tsv', show_default=True,
    help="Output file format"
)
def estimate_kiv(exon1a, exon1b, output, format):
    """
    Estimate LPA KIV2 copy numbers from exon1A and exon1B diploid copy numbers.

    Formula:
        diploid_estimate = 34.9 × exon1A + 5.2 × exon1B - 1  
        haploid_estimate = diploid_estimate / 2
    """
    from .utils.estimate_kiv import compute_kiv_estimates
    print_banner()
    
    try:
        # Compute estimates - pass string paths
        _ = compute_kiv_estimates(
            output=output,
            exon1a_file=exon1a,
            exon1b_file=exon1b,
            console=console
        )
    
    except Exception as e:
        console.print(f"[red]✗ KIV2 estimation failed: {e}[/red]")
        sys.exit(1)

# In[ ]: Main Entry Point
def main():
    """Entry point for CLI."""
    try:
        cli()
    except KeyboardInterrupt:
        console.print("\n[yellow]Pipeline interrupted by user[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"\n[bold red]Unexpected error:[/bold red] {str(e)}")
        sys.exit(1)

# In[ ]: Run Main
if __name__ == '__main__':
    main()