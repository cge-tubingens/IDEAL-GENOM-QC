"""Shared utility functions for genomic analysis pipelines."""

import os
import psutil
from pathlib import Path


def get_optimal_threads(reserve: int = 2, default: int = 10) -> int:
    """
    Calculate optimal thread count for genomic analysis operations.
    
    Determines the number of threads to use based on available CPU cores,
    reserving some cores for system operations.
    
    Parameters
    ----------
    reserve : int, default=2
        Number of cores to reserve for system operations
    default : int, default=10
        Default thread count if CPU detection fails
        
    Returns
    -------
    int
        Optimal number of threads to use
        
    Examples
    --------
    >>> threads = get_optimal_threads()  # On 16-core system, returns 14
    >>> threads = get_optimal_threads(reserve=4)  # Returns 12
    """
    cpu_count = os.cpu_count()
    if cpu_count is not None:
        return max(1, cpu_count - reserve)
    # Fallback: use half of logical cores or default
    return max(1, (psutil.cpu_count(logical=True) or default) // 2)


def get_available_memory(fraction: float = 2/3) -> int:
    """
    Calculate available memory for genomic analysis operations.
    
    Determines the amount of memory to allocate based on currently available
    system memory, using a configurable fraction to avoid system instability.
    
    Parameters
    ----------
    fraction : float, default=2/3
        Fraction of available memory to use (should be between 0 and 1)
        
    Returns
    -------
    int
        Memory in MB to allocate
        
    Raises
    ------
    ValueError
        If fraction is not between 0 and 1
        
    Examples
    --------
    >>> memory_mb = get_available_memory()  # Uses 2/3 of available memory
    >>> memory_mb = get_available_memory(fraction=0.5)  # Uses half
    """
    if not 0 < fraction <= 1:
        raise ValueError(f"fraction must be between 0 and 1, got {fraction}")
    
    memory_info = psutil.virtual_memory()
    available_memory_mb = memory_info.available / (1024 * 1024)
    return int(round(fraction * available_memory_mb, 0))


def count_file_lines(file_path: Path) -> int:
    """
    Count lines in a file efficiently.
    
    Uses a generator expression for memory-efficient line counting,
    suitable for large genomic data files.
    
    Parameters
    ----------
    file_path : Path
        Path to the file to count
        
    Returns
    -------
    int
        Number of lines in the file
        
    Raises
    ------
    FileNotFoundError
        If the file does not exist
    IOError
        If the file cannot be read
        
    Examples
    --------
    >>> from pathlib import Path
    >>> count = count_file_lines(Path('variants.bim'))
    """
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    with open(file_path, 'r') as f:
        return sum(1 for _ in f)


def validate_output_dir(output_dir: Path, create: bool = True) -> Path:
    """
    Validate and optionally create output directory.
    
    Parameters
    ----------
    output_dir : Path
        Path to the output directory
    create : bool, default=True
        If True, create the directory if it doesn't exist
        
    Returns
    -------
    Path
        Validated output directory path
        
    Raises
    ------
    FileNotFoundError
        If directory doesn't exist and create=False
    PermissionError
        If directory cannot be created due to permissions
        
    Examples
    --------
    >>> from pathlib import Path
    >>> output = validate_output_dir(Path('/data/results'))
    """
    if not output_dir.exists():
        if create:
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(f"Output directory does not exist: {output_dir}")
    
    if not output_dir.is_dir():
        raise NotADirectoryError(f"Path exists but is not a directory: {output_dir}")
    
    return output_dir


def format_memory_size(bytes_size: int) -> str:
    """
    Format byte size into human-readable string.
    
    Parameters
    ----------
    bytes_size : int
        Size in bytes
        
    Returns
    -------
    str
        Formatted size string (e.g., '1.5 GB', '256 MB')
        
    Examples
    --------
    >>> format_memory_size(1536 * 1024 * 1024)
    '1.50 GB'
    >>> format_memory_size(512 * 1024)
    '512.00 KB'
    """
    size = float(bytes_size)
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return f"{size:.2f} {unit}"
        size /= 1024.0
    return f"{size:.2f} PB"
