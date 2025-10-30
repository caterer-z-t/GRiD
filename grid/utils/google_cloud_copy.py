# grid/utils/google_cloud_copy.py
# In[1]: Imports
import os
import subprocess
from pathlib import Path
from rich.console import Console

# In[2]: Setup Console
console = Console()

# In[3]: Function to Copy CRAMs from Google Cloud Bucket
def copy_crams_from_bucket(bucket_path: str, local_dir: str):
    local_dir = Path(local_dir).expanduser().resolve()
    local_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"[bold green]Scanning bucket:[/bold green] {bucket_path}")

    # Get the current environment
    env = os.environ.copy()

    # Ensure gcloud is in PATH
    try:
        gcloud_path = subprocess.run(
            ["which", "gcloud"],
            capture_output=True,
            text=True,
            check=True,
            shell=False
        ).stdout.strip()
        if gcloud_path:
            gcloud_dir = str(Path(gcloud_path).parent)
            env["PATH"] = gcloud_dir + os.pathsep + env.get("PATH", "")
    except subprocess.CalledProcessError:
        console.print(f"[red]gcloud command not found. Please install Google Cloud SDK.[/red]")
        return

    # List all files in bucket
    try:
        result = subprocess.run(
            ["gcloud", "storage", "ls", "-r", bucket_path],
            capture_output=True,
            text=True,
            check=True,
            env=env,
            shell=False
        )
        all_files = result.stdout.splitlines()
    except subprocess.CalledProcessError as e:
        console.print(f"[red]Failed to list bucket contents:[/red] {e.stderr}")
        return
    except FileNotFoundError:
        console.print(f"[yellow]Current PATH: {env.get('PATH', 'Not set')}[/yellow]")
        return

    # Filter CRAM files
    cram_files = [f for f in all_files if f.endswith(".cram")]

    for cram_file in cram_files:
        file_name = Path(cram_file).name
        dest_path = local_dir / file_name

        if dest_path.exists():
            console.print(f"[yellow]Skipping {file_name} (already exists)[/yellow]")
            continue

        console.print(f"[blue]Copying {cram_file} → {local_dir}[/blue]")
        try:
            subprocess.run(
                ["gcloud", "storage", "cp", cram_file, str(dest_path)],
                check=True,
                env=env  # <-- env is now defined at function scope
            )
        except subprocess.CalledProcessError as e:
            console.print(f"[red]Failed to copy {cram_file}: {e}[/red]")

    # Remove empty files
    for f in local_dir.glob("*"):
        if f.stat().st_size == 0:
            console.print(f"[red]Removing empty file: {f.name}[/red]")
            f.unlink()

    console.print(f"[green]CRAM copy complete! Files in:[/green] {local_dir}")