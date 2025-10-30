# grid/utils/helper_dir/display_results.py
# In[1]: Imports
from rich.console import Console

# In[2]: Console Setup
console = Console()

# In[3]: Print Individual Success Message
def print_individual_success(
    basename: str,
    message: str,
    progress_console=None
) -> None:
    """
    Print success message for an individual item.
    
    Args:
        basename: Item basename/name
        message: Success message (e.g., "done: reads=1234" or "indexed")
        progress_console: Optional Progress console for proper rendering
    """
    output = f"[green]✓ {basename} {message}[/green]"
    
    if progress_console:
        progress_console.print(output)
    else:
        console.print(output)

# In[4]: Print Individual Error Message
def print_individual_error(
    basename: str,
    error: str,
    progress_console=None
) -> None:
    """
    Print error message for an individual item.
    
    Args:
        basename: Item basename/name
        error: Error message
        progress_console: Optional Progress console for proper rendering
    """
    output = f"[red]✗ {basename} failed: {error}[/red]"
    
    if progress_console:
        progress_console.print(output)
    else:
        console.print(output)

# In[5]: Print Batch Summary
def print_batch_summary(
    total_processed: int,
    failed_items: list,
    operation: str = "processed"
) -> None:
    """
    Print summary of a batch operation.
    
    Args:
        total_processed: Total number of items processed
        failed_items: List of failed item names/basenames
        operation: Description of operation (e.g., "indexed", "processed", "counted")
    """
    if failed_items:
        console.print(f"[red]✗ {len(failed_items)} files failed to be {operation}[/red]")
    else:
        console.print(f"[bold green]✓ Successfully {operation} {total_processed} files[/bold green]")