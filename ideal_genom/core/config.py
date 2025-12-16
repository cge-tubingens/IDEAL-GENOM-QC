"""Configuration loading and validation for ideal_genom pipelines."""

import yaml
from pathlib import Path
from typing import Dict, Any, List
import logging

logger = logging.getLogger(__name__)


class ConfigurationError(Exception):
    """Raised when configuration is invalid."""
    pass


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load pipeline configuration from YAML file.
    
    Parameters
    ----------
    config_path : str
        Path to YAML configuration file
        
    Returns
    -------
    dict
        Parsed configuration dictionary
        
    Raises
    ------
    ConfigurationError
        If configuration is invalid or file not found
    """
    config_file = Path(config_path)
    
    if not config_file.exists():
        raise ConfigurationError(f"Configuration file not found: {config_path}")
    
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ConfigurationError(f"Failed to parse YAML: {e}")
    
    validate_config(config)
    return config


def validate_config(config: Dict[str, Any]) -> None:
    """
    Validate pipeline configuration structure.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary to validate
        
    Raises
    ------
    ConfigurationError
        If configuration structure is invalid
    """
    if 'pipeline' not in config:
        raise ConfigurationError("Configuration must contain 'pipeline' key")
    
    pipeline = config['pipeline']
    
    # Check required pipeline fields
    required_fields = ['name', 'base_output_dir', 'steps']
    for field in required_fields:
        if field not in pipeline:
            raise ConfigurationError(f"Pipeline must contain '{field}' key")
    
    # Validate steps
    if not isinstance(pipeline['steps'], list):
        raise ConfigurationError("'steps' must be a list")
    
    if len(pipeline['steps']) == 0:
        raise ConfigurationError("Pipeline must contain at least one step")
    
    for i, step in enumerate(pipeline['steps']):
        validate_step(step, i)


def validate_step(step: Dict[str, Any], index: int) -> None:
    """
    Validate a single pipeline step configuration.
    
    Parameters
    ----------
    step : dict
        Step configuration dictionary
    index : int
        Step index in pipeline (for error messages)
        
    Raises
    ------
    ConfigurationError
        If step configuration is invalid
    """
    required_fields = ['name', 'module', 'class', 'init_params']
    
    for field in required_fields:
        if field not in step:
            raise ConfigurationError(
                f"Step {index} ('{step.get('name', 'unnamed')}') "
                f"missing required field: '{field}'"
            )
    
    # Validate init_params
    init_params = step['init_params']
    required_init = ['input_path', 'input_name', 'output_path']
    
    for param in required_init:
        if param not in init_params:
            raise ConfigurationError(
                f"Step {index} ('{step['name']}') init_params "
                f"missing required parameter: '{param}'"
            )