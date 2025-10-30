"""Extract reference FASTA sequences from a BED file."""
import sys
import time
from pathlib import Path
from rich.console import Console
from pyfaidx import Fasta

console = Console()

def extract_reference(reference_fa: str, bed_file: str, output_dir: str, output_prefix: str = "ref_lpa"):
    """
    Extracts sequences from a reference FASTA for regions listed in a BED file.

    Args:
        reference_fa (str): Path to reference genome FASTA.
        bed_file (str): Path to BED file with regions to extract.
        output_dir (str): Directory to write output FASTA file.
        output_prefix (str): Output FASTA file prefix (default: ref_lpa).
    """
    start_time = time.time()
    ref_path = Path(reference_fa)
    bed_path = Path(bed_file)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    output_fasta = out_dir / f"{output_prefix}.fasta"

    console.print(f"[bold blue]Starting reference extraction[/bold blue]")
    console.print(f"Reference: {ref_path}")
    console.print(f"BED file:  {bed_path}")
    console.print(f"Output:    {output_fasta}")

    try:
        ref = Fasta(str(ref_path))
        with open(bed_path) as bed, open(output_fasta, "w") as out_fa:
            for line in bed:
                if line.startswith("#") or not line.strip():
                    continue
                chrom, start, end, *rest = line.strip().split("\t")
                seq_id = f"{chrom}:{start}-{end}"
                start, end = int(start), int(end)
                seq = ref[chrom][start:end].seq
                out_fa.write(f">{seq_id}\n{seq}\n")

        elapsed = time.time() - start_time
        console.print(f"[green]✓ Extraction complete[/green] ({elapsed:.2f}s)")
        console.print(f"[bold green]Output:[/bold green] {output_fasta}")

        # Show preview of output
        with open(output_fasta) as f:
            preview = "".join(f.readlines()[:10])
        console.print("\n[yellow]Preview:[/yellow]\n" + preview)

        return str(output_fasta)

    except Exception as e:
        console.print(f"[red]✗ Reference extraction failed:[/red] {e}")
        sys.exit(1)
