"""Environment verification utilities for genomic analysis software."""

import subprocess
import shutil
import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class SoftwareInfo:
    """Information about installed software."""
    name: str
    installed: bool
    path: Optional[str] = None
    version: Optional[str] = None
    error: Optional[str] = None


class EnvironmentError(Exception):
    """Raised when required software is not available."""
    pass

def check_software_installed(software: str) -> bool:
    """
    Check if a software tool is installed and available in PATH.
    
    Parameters
    ----------
    software : str
        Name of the software to check
        
    Returns
    -------
    bool
        True if software is installed, False otherwise
    """
    return shutil.which(software) is not None

def get_software_version(
    software: str,
    version_flag: str = '--version'
) -> Optional[str]:
    """
    Get version string of installed software.
    
    Parameters
    ----------
    software : str
        Name of the software
    version_flag : str, default='--version'
        Flag to get version information
        
    Returns
    -------
    str or None
        Version string if available, None otherwise
    """
    try:
        result = subprocess.run(
            [software, version_flag],
            capture_output=True,
            text=True,
            timeout=5
        )
        
        if result.returncode == 0:
            # Try to extract version from output
            output = result.stdout + result.stderr
            return output.strip().split('\n')[0]
        
    except (subprocess.SubprocessError, FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    return None

def verify_software(
    software: str,
    version_flag: str = '--version',
    required: bool = True
) -> SoftwareInfo:
    """
    Verify that a software tool is installed and get its info.
    
    Parameters
    ----------
    software : str
        Name of the software to verify
    version_flag : str, default='--version'
        Flag to get version information
    required : bool, default=True
        If True, log as error when not found; otherwise log as warning
        
    Returns
    -------
    SoftwareInfo
        Information about the software installation
    """
    path = shutil.which(software)
    
    if path is None:
        error_msg = f"{software} not found in PATH"
        
        if required:
            logger.error(error_msg)
        else:
            logger.warning(error_msg)
        
        return SoftwareInfo(
            name=software,
            installed=False,
            error=error_msg
        )
    
    version = get_software_version(software, version_flag)
    
    logger.info(f"✓ {software} found: {path}")
    if version:
        logger.info(f"  Version: {version}")
    
    return SoftwareInfo(
        name=software,
        installed=True,
        path=path,
        version=version
    )

def verify_plink_installation() -> Tuple[bool, bool]:
    """
    Verify PLINK installation and check for both PLINK 1.9 and PLINK2.
    
    PLINK 1.9 and PLINK2 are typically installed together, with commands
    'plink' (PLINK 1.9) and 'plink2' (PLINK2).
    
    Returns
    -------
    tuple of (bool, bool)
        (plink19_installed, plink2_installed)
    """
    plink_installed = check_software_installed('plink')
    plink2_installed = check_software_installed('plink2')
    
    if plink_installed:
        version = get_software_version('plink', '--version')
        logger.info(f"✓ PLINK 1.9 found")
        if version:
            logger.info(f"  Version: {version}")
    
    if plink2_installed:
        version = get_software_version('plink2', '--version')
        logger.info(f"✓ PLINK2 found")
        if version:
            logger.info(f"  Version: {version}")
    
    return plink_installed, plink2_installed

def verify_genomic_environment(
    required_tools: Optional[List[str]] = None,
    optional_tools: Optional[List[str]] = None
) -> Dict[str, SoftwareInfo]:
    """
    Verify that all required genomic analysis tools are installed.
    
    Parameters
    ----------
    required_tools : list of str, optional
        List of required software tools. 
        Default: ['plink', 'plink2', 'bcftools']
    optional_tools : list of str, optional
        List of optional software tools. 
        Default: ['gcta64', 'king']
        
    Returns
    -------
    dict
        Dictionary mapping software names to SoftwareInfo objects
        
    Raises
    ------
    EnvironmentError
        If any required tool is not installed
        
    Examples
    --------
    >>> verify_genomic_environment()
    >>> # Or check specific tools
    >>> verify_genomic_environment(
    ...     required_tools=['plink', 'bcftools'],
    ...     optional_tools=['gcta64', 'king']
    ... )
    """
    if required_tools is None:
        required_tools = ['plink', 'plink2', 'bcftools']
    
    if optional_tools is None:
        optional_tools = ['gcta64', 'king']
    
    logger.info("Verifying genomic analysis environment...")
    logger.info("=" * 60)
    
    results = {}
    missing_required = []
    
    # Check required tools
    logger.info("\nRequired tools:")
    for tool in required_tools:
        info = verify_software(tool, required=True)
        results[tool] = info
        
        if not info.installed:
            missing_required.append(tool)
    
    # Check optional tools
    if optional_tools:
        logger.info("\nOptional tools:")
        for tool in optional_tools:
            info = verify_software(tool, required=False)
            results[tool] = info
    
    logger.info("\n" + "=" * 60)
    
    # Raise error if required tools are missing
    if missing_required:
        error_msg = (
            f"Missing required software: {', '.join(missing_required)}\n"
            f"Please install the missing tools and ensure they are in your PATH."
        )
        logger.error(error_msg)
        raise EnvironmentError(error_msg)
    
    logger.info("✓ All required tools are installed")
    
    return results

def check_plink_versions() -> Tuple[Optional[str], Optional[str]]:
    """
    Check versions of both PLINK 1.9 and PLINK2.
    
    Returns
    -------
    tuple of (str or None, str or None)
        (plink19_version, plink2_version)
    """
    plink_version = get_software_version('plink', '--version')
    plink2_version = get_software_version('plink2', '--version')
    
    return plink_version, plink2_version

# Quick verification function for use in pipelines
def quick_verify(tools: Optional[List[str]] = None) -> None:
    """
    Quick verification of required tools with minimal output.
    
    Parameters
    ----------
    tools : list of str, optional
        Tools to verify. Default: ['plink', 'plink2', 'bcftools']
        
    Raises
    ------
    EnvironmentError
        If any tool is not installed
        
    Examples
    --------
    >>> # Verify default tools
    >>> quick_verify()
    
    >>> # Verify specific tools
    >>> quick_verify(['plink', 'gcta64'])
    """
    if tools is None:
        tools = ['plink', 'plink2', 'bcftools']
    
    missing = []
    for tool in tools:
        if not check_software_installed(tool):
            missing.append(tool)
    
    if missing:
        raise EnvironmentError(
            f"Missing required software: {', '.join(missing)}"
        )
