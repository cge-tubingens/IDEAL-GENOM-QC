"""Command execution utilities for external genomic tools."""

import subprocess
import logging
import shlex
from typing import Optional, List, Union
from pathlib import Path

logger = logging.getLogger(__name__)


class CommandExecutionError(Exception):
    """Raised when a shell command fails."""
    pass


def shell_do(
    command: Union[str, List[str]],
    cwd: Optional[str] = None,
    log_file: Optional[str] = None,
    capture_output: bool = False,
    check: bool = True
) -> subprocess.CompletedProcess:
    """
    Execute a shell command for genomic analysis tools.
    
    This is a wrapper around subprocess.run with logging and error handling
    tailored for genomic analysis pipelines (PLINK, GCTA, bcftools, etc.).
    
    Parameters
    ----------
    command : str or list of str
        Command to execute. Can be a string or list of arguments.
    cwd : str, optional
        Working directory for command execution
    log_file : str, optional
        Path to file where stdout/stderr should be logged
    capture_output : bool, default=False
        If True, capture stdout and stderr in returned object
    check : bool, default=True
        If True, raise CommandExecutionError on non-zero exit code
        
    Returns
    -------
    subprocess.CompletedProcess
        Completed process with returncode, stdout, stderr
        
    Raises
    ------
    CommandExecutionError
        If command fails and check=True
        
    Examples
    --------
    >>> # Execute PLINK command
    >>> shell_do("plink --bfile input --maf 0.01 --make-bed --out output")
    
    >>> # Execute with working directory
    >>> shell_do(
    ...     ["bcftools", "view", "-Oz", "input.vcf"],
    ...     cwd="/data/work",
    ...     log_file="/data/logs/bcftools.log"
    ... )
    """
    # Convert string command to list if needed
    if isinstance(command, str):
        cmd_list = shlex.split(command)
        cmd_str = command
    else:
        cmd_list = command
        cmd_str = ' '.join(command)
    
    logger.info(f"Executing: {cmd_str}")
    if cwd:
        logger.debug(f"Working directory: {cwd}")
    
    # Prepare subprocess arguments
    subprocess_kwargs = {
        'args': cmd_list,
        'cwd': cwd,
        'text': True,
    }
    
    # Handle output capture/logging
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(log_file, 'w') as f:
            subprocess_kwargs['stdout'] = f
            subprocess_kwargs['stderr'] = subprocess.STDOUT
            
            try:
                result = subprocess.run(**subprocess_kwargs, check=False)
            except Exception as e:
                logger.error(f"Command execution error: {e}")
                if check:
                    raise CommandExecutionError(f"Failed to execute: {cmd_str}\nError: {e}")
                # Create a mock CompletedProcess for error case
                result = subprocess.CompletedProcess(cmd_list, returncode=1)
    else:
        if capture_output:
            subprocess_kwargs['capture_output'] = True
        
        try:
            result = subprocess.run(**subprocess_kwargs, check=False)
        except Exception as e:
            logger.error(f"Command execution error: {e}")
            raise CommandExecutionError(f"Failed to execute: {cmd_str}\nError: {e}")
    
    # Check return code
    if result.returncode != 0:
        error_msg = f"Command failed with exit code {result.returncode}: {cmd_str}"
        logger.error(error_msg)
        
        if capture_output and result.stderr:
            logger.error(f"stderr: {result.stderr}")
        
        if check:
            raise CommandExecutionError(error_msg)
    else:
        logger.info(f"Command completed successfully")
    
    return result


def run_plink(
    args: List[str],
    log_file: Optional[str] = None,
    cwd: Optional[str] = None
) -> subprocess.CompletedProcess:
    """
    Execute PLINK command.
    
    Parameters
    ----------
    args : list of str
        PLINK arguments (without 'plink' command itself)
    log_file : str, optional
        Path to log file
    cwd : str, optional
        Working directory
        
    Returns
    -------
    subprocess.CompletedProcess
        Completed process
        
    Examples
    --------
    >>> run_plink([
    ...     '--bfile', 'input',
    ...     '--maf', '0.01',
    ...     '--make-bed',
    ...     '--out', 'output'
    ... ])
    """
    command = ['plink'] + args
    return shell_do(command, cwd=cwd, log_file=log_file)


def run_plink2(
    args: List[str],
    log_file: Optional[str] = None,
    cwd: Optional[str] = None
) -> subprocess.CompletedProcess:
    """
    Execute PLINK2 command.
    
    Parameters
    ----------
    args : list of str
        PLINK2 arguments (without 'plink2' command itself)
    log_file : str, optional
        Path to log file
    cwd : str, optional
        Working directory
        
    Returns
    -------
    subprocess.CompletedProcess
        Completed process
    """
    command = ['plink2'] + args
    return shell_do(command, cwd=cwd, log_file=log_file)


def run_gcta(
    args: List[str],
    log_file: Optional[str] = None,
    cwd: Optional[str] = None
) -> subprocess.CompletedProcess:
    """
    Execute GCTA command.
    
    Parameters
    ----------
    args : list of str
        GCTA arguments (without 'gcta64' command itself)
    log_file : str, optional
        Path to log file
    cwd : str, optional
        Working directory
        
    Returns
    -------
    subprocess.CompletedProcess
        Completed process
    """
    command = ['gcta64'] + args
    return shell_do(command, cwd=cwd, log_file=log_file)


def run_bcftools(
    args: List[str],
    log_file: Optional[str] = None,
    cwd: Optional[str] = None
) -> subprocess.CompletedProcess:
    """
    Execute bcftools command.
    
    Parameters
    ----------
    args : list of str
        bcftools arguments (without 'bcftools' command itself)
    log_file : str, optional
        Path to log file
    cwd : str, optional
        Working directory
        
    Returns
    -------
    subprocess.CompletedProcess
        Completed process
    """
    command = ['bcftools'] + args
    return shell_do(command, cwd=cwd, log_file=log_file)
