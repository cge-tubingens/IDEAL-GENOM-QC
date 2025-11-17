"""Shared utility functions for genomic analysis pipelines."""

import os
import psutil
import logging
import contextlib
import tempfile
import shutil
import weakref
from pathlib import Path
from typing import Optional, Union, List, Set, Generator, Iterator, Any, IO
from collections.abc import Generator as GeneratorABC

logger = logging.getLogger(__name__)


def get_optimal_threads(reserve: int = 2, default: int = 10, max_threads: Optional[int] = None) -> int:
    """
    Calculate optimal thread count for genomic analysis operations.
    
    Determines the number of threads to use based on available CPU cores,
    reserving some cores for system operations with robust fallback handling.
    
    Parameters
    ----------
    reserve : int, default=2
        Number of cores to reserve for system operations
    default : int, default=10
        Default thread count if CPU detection fails
    max_threads : int, optional
        Maximum number of threads to use (caps the result)
        
    Returns
    -------
    int
        Optimal number of threads to use (always >= 1)
        
    Examples
    --------
    >>> threads = get_optimal_threads()  # On 16-core system, returns 14
    >>> threads = get_optimal_threads(reserve=4)  # Returns 12
    >>> threads = get_optimal_threads(max_threads=8)  # Never exceeds 8
    """
    if reserve < 0:
        raise ValueError(f"reserve must be >= 0, got {reserve}")
    if default < 1:
        raise ValueError(f"default must be >= 1, got {default}")
    if max_threads is not None and max_threads < 1:
        raise ValueError(f"max_threads must be >= 1, got {max_threads}")
    
    # Try multiple methods to get CPU count
    cpu_count = None
    
    # Method 1: os.cpu_count() - preferred
    try:
        cpu_count = os.cpu_count()
        if cpu_count and cpu_count > 0:
            logger.debug(f"CPU count from os.cpu_count(): {cpu_count}")
        else:
            cpu_count = None
    except Exception as e:
        logger.warning(f"Failed to get CPU count from os.cpu_count(): {e}")
        cpu_count = None
    
    # Method 2: psutil fallback
    if cpu_count is None:
        try:
            cpu_count = psutil.cpu_count(logical=True)
            if cpu_count and cpu_count > 0:
                logger.debug(f"CPU count from psutil.cpu_count(): {cpu_count}")
            else:
                cpu_count = None
        except Exception as e:
            logger.warning(f"Failed to get CPU count from psutil: {e}")
            cpu_count = None
    
    # Calculate optimal threads
    if cpu_count is not None and cpu_count > 0:
        optimal = max(1, cpu_count - reserve)
    else:
        logger.warning(f"Could not detect CPU count, using default: {default}")
        optimal = default
    
    # Apply max_threads cap if specified
    if max_threads is not None:
        optimal = min(optimal, max_threads)
    
    logger.debug(f"Calculated optimal threads: {optimal} (CPU: {cpu_count}, reserve: {reserve})")
    return optimal


def get_available_memory(
    fraction: float = 2/3, 
    min_mb: int = 512, 
    max_mb: Optional[int] = None,
    safety_buffer_mb: int = 1024
) -> int:
    """
    Calculate available memory for genomic analysis operations with safety checks.
    
    Determines the amount of memory to allocate based on currently available
    system memory, using a configurable fraction with minimum and maximum limits
    to avoid system instability.
    
    Parameters
    ----------
    fraction : float, default=2/3
        Fraction of available memory to use (should be between 0 and 1)
    min_mb : int, default=512
        Minimum memory to allocate in MB
    max_mb : int, optional
        Maximum memory to allocate in MB (None for no limit)
    safety_buffer_mb : int, default=1024
        Safety buffer to always leave available for system (MB)
        
    Returns
    -------
    int
        Memory in MB to allocate
        
    Raises
    ------
    ValueError
        If fraction is not between 0 and 1, or other parameters are invalid
    RuntimeError
        If insufficient memory is available
        
    Examples
    --------
    >>> memory_mb = get_available_memory()  # Uses 2/3 of available memory
    >>> memory_mb = get_available_memory(fraction=0.5, max_mb=8192)  # Uses half, max 8GB
    """
    if not 0 < fraction <= 1:
        raise ValueError(f"fraction must be between 0 and 1, got {fraction}")
    if min_mb < 0:
        raise ValueError(f"min_mb must be >= 0, got {min_mb}")
    if max_mb is not None and max_mb < min_mb:
        raise ValueError(f"max_mb ({max_mb}) must be >= min_mb ({min_mb})")
    if safety_buffer_mb < 0:
        raise ValueError(f"safety_buffer_mb must be >= 0, got {safety_buffer_mb}")
    
    try:
        memory_info = psutil.virtual_memory()
        
        # Log memory information for debugging
        total_mb = memory_info.total / (1024 * 1024)
        available_mb = memory_info.available / (1024 * 1024)
        used_mb = memory_info.used / (1024 * 1024)
        
        logger.debug(f"System memory - Total: {total_mb:.0f}MB, Available: {available_mb:.0f}MB, Used: {used_mb:.0f}MB")
        
        # Calculate usable memory after safety buffer
        usable_memory_mb = max(0, available_mb - safety_buffer_mb)
        
        # Apply fraction
        target_memory_mb = usable_memory_mb * fraction
        
        # Apply minimum constraint
        if target_memory_mb < min_mb:
            if usable_memory_mb < min_mb:
                raise RuntimeError(
                    f"Insufficient memory available. Need at least {min_mb}MB, "
                    f"but only {usable_memory_mb:.0f}MB usable after safety buffer."
                )
            target_memory_mb = min_mb
            logger.warning(f"Using minimum memory {min_mb}MB (would be {target_memory_mb:.0f}MB with fraction)")
        
        # Apply maximum constraint if specified
        if max_mb is not None and target_memory_mb > max_mb:
            target_memory_mb = max_mb
            logger.debug(f"Capped memory allocation to maximum {max_mb}MB")
        
        result = int(round(target_memory_mb, 0))
        logger.debug(f"Allocated memory: {result}MB (fraction: {fraction}, available: {available_mb:.0f}MB)")
        
        return result
        
    except Exception as e:
        logger.error(f"Failed to get memory information: {e}")
        raise RuntimeError(f"Could not determine available memory: {e}")


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

def get_system_resource_info() -> dict:
    """
    Get comprehensive system resource information.
    
    Returns
    -------
    dict
        Dictionary containing CPU, memory, and disk information
        
    Examples
    --------
    >>> info = get_system_resource_info()
    >>> print(f"Available memory: {info['memory']['available_mb']:.0f} MB")
    """
    try:
        # CPU information
        cpu_count_physical = psutil.cpu_count(logical=False)
        cpu_count_logical = psutil.cpu_count(logical=True)
        cpu_freq = psutil.cpu_freq()
        
        # Memory information
        memory = psutil.virtual_memory()
        swap = psutil.swap_memory()
        
        # Disk information for current working directory
        disk = psutil.disk_usage('.')
        
        return {
            'cpu': {
                'physical_cores': cpu_count_physical,
                'logical_cores': cpu_count_logical,
                'frequency_mhz': cpu_freq.current if cpu_freq else None,
            },
            'memory': {
                'total_mb': memory.total / (1024 * 1024),
                'available_mb': memory.available / (1024 * 1024),
                'used_mb': memory.used / (1024 * 1024),
                'percent_used': memory.percent,
            },
            'swap': {
                'total_mb': swap.total / (1024 * 1024),
                'used_mb': swap.used / (1024 * 1024),
                'percent_used': swap.percent,
            },
            'disk': {
                'total_gb': disk.total / (1024 * 1024 * 1024),
                'free_gb': disk.free / (1024 * 1024 * 1024),
                'used_gb': disk.used / (1024 * 1024 * 1024),
                'percent_used': (disk.used / disk.total) * 100,
            }
        }
    except Exception as e:
        logger.error(f"Failed to get system resource information: {e}")
        return {}
