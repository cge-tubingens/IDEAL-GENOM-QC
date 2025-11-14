"""Pipeline orchestration engine for ideal_genom."""

import os
import re
import importlib
import logging

from typing import Dict, Any, Optional
from pathlib import Path


logger = logging.getLogger(__name__)


class PipelineExecutor:
    """
    Orchestrates execution of sub-pipeline classes based on configuration.
    
    Attributes
    ----------
    config : dict
        Pipeline configuration dictionary
    steps : dict
        Dictionary of instantiated sub-pipeline objects
    base_output_dir : str
        Base output directory for all pipeline steps
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize pipeline executor.
        
        Parameters
        ----------
        config : dict
            Pipeline configuration dictionary (from config.load_config)
        """
        self.config = config
        self.steps = {}  # Store instantiated sub-pipeline objects
        self.base_output_dir = config['pipeline']['base_output_dir']
        self.pipeline_name = config['pipeline']['name']
        
        # Create base output directory
        os.makedirs(self.base_output_dir, exist_ok=True)
        
        # Setup logging
        self.logger = self._setup_logging()
        
    def execute(self) -> None:
        """Execute all pipeline steps sequentially."""
        pipeline_steps = self.config['pipeline']['steps']
        
        self.logger.info(f"Starting pipeline: {self.pipeline_name}")
        self.logger.info(f"Output directory: {self.base_output_dir}")
        self.logger.info(f"Total steps: {len(pipeline_steps)}")
        
        for i, step_config in enumerate(pipeline_steps, 1):
            step_name = step_config['name']
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"Step {i}/{len(pipeline_steps)}: {step_name}")
            self.logger.info(f"{'='*60}")
            
            try:
                self._execute_step(step_config)
                self.logger.info(f"✓ Completed step: {step_name}")
            except Exception as e:
                self.logger.error(f"✗ Failed step: {step_name}")
                self.logger.error(f"Error: {str(e)}", exc_info=True)
                raise RuntimeError(f"Pipeline failed at step '{step_name}': {str(e)}")
        
        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"Pipeline completed successfully: {self.pipeline_name}")
        self.logger.info(f"{'='*60}")
    
    def _execute_step(self, step_config: Dict[str, Any]) -> None:
        """
        Execute a single pipeline step.
        
        Parameters
        ----------
        step_config : dict
            Configuration for the step to execute
        """
        step_name = step_config['name']
        
        # Import the class dynamically
        module_path = step_config['module']
        class_name = step_config['class']
        
        self.logger.info(f"Loading {class_name} from {module_path}")
        
        try:
            module = importlib.import_module(module_path)
            pipeline_class = getattr(module, class_name)
        except (ImportError, AttributeError) as e:
            raise ImportError(
                f"Failed to import {class_name} from {module_path}: {str(e)}"
            )
        
        # Resolve parameters (handle ${...} references)
        init_params = self._resolve_params(step_config['init_params'])
        execute_params = self._resolve_params(step_config.get('execute_params', {}))
        
        self.logger.info(f"Initializing {class_name}")
        self.logger.debug(f"Init params: {init_params}")
        self.logger.debug(f"Execute params: {execute_params}")
        
        # Instantiate the sub-pipeline class
        pipeline_instance = pipeline_class(**init_params)
        
        # Determine the execute method name
        # Convention: execute_<step_name>_pipeline
        execute_method_name = f'execute_{step_name}_pipeline'
        
        if not hasattr(pipeline_instance, execute_method_name):
            raise AttributeError(
                f"{class_name} does not have method '{execute_method_name}'"
            )
        
        execute_method = getattr(pipeline_instance, execute_method_name)
        
        self.logger.info(f"Executing {execute_method_name}")
        execute_method(execute_params)
        
        # Store the instance for reference by subsequent steps
        self.steps[step_name] = pipeline_instance
        
        self.logger.info(f"Step output stored as: steps.{step_name}")
    
    def _resolve_params(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Resolve parameter values, including references to previous steps.
        
        Parameters
        ----------
        params : dict
            Parameters dictionary potentially containing references
            
        Returns
        -------
        dict
            Parameters with all references resolved
        """
        resolved = {}
        
        for key, value in params.items():
            if isinstance(value, str):
                resolved[key] = self._resolve_string_value(value)
            elif isinstance(value, dict):
                resolved[key] = self._resolve_params(value)
            elif isinstance(value, list):
                resolved[key] = [
                    self._resolve_string_value(v) if isinstance(v, str) else v
                    for v in value
                ]
            else:
                resolved[key] = value
        
        return resolved
    
    def _resolve_string_value(self, value: str) -> Any:
        """
        Resolve a string value that may contain references.
        
        Supports:
        - ${base_output_dir} - pipeline base output directory
        - ${steps.step_name.attribute} - attribute from previous step
        
        Parameters
        ----------
        value : str
            String value potentially containing references
            
        Returns
        -------
        Any
            Resolved value
        """
        # Find all ${...} patterns
        pattern = r'\$\{([^}]+)\}'
        matches = re.findall(pattern, value)
        
        if not matches:
            return value
        
        resolved_value = value
        
        for match in matches:
            ref_value = self._resolve_reference(match)
            resolved_value = resolved_value.replace(f'${{{match}}}', str(ref_value))
        
        return resolved_value
    
    def _resolve_reference(self, reference: str) -> Any:
        """
        Resolve a reference string to its actual value.
        
        Parameters
        ----------
        reference : str
            Reference string (without ${ })
            
        Returns
        -------
        Any
            Resolved value
            
        Raises
        ------
        ValueError
            If reference is invalid or step not found
        """
        # Handle base_output_dir
        if reference == 'base_output_dir':
            return self.base_output_dir
        
        # Handle steps.step_name.attribute
        if reference.startswith('steps.'):
            parts = reference.split('.')
            
            if len(parts) < 3:
                raise ValueError(
                    f"Invalid step reference: {reference}. "
                    f"Expected format: steps.step_name.attribute"
                )
            
            step_name = parts[1]
            attr_path = parts[2:]
            
            if step_name not in self.steps:
                raise ValueError(
                    f"Step '{step_name}' not found. "
                    f"Available steps: {list(self.steps.keys())}"
                )
            
            # Navigate through nested attributes
            obj = self.steps[step_name]
            for attr in attr_path:
                if not hasattr(obj, attr):
                    raise ValueError(
                        f"Step '{step_name}' does not have attribute '{attr}'"
                    )
                obj = getattr(obj, attr)
            
            return obj
        
        raise ValueError(f"Unknown reference: {reference}")
    
    def _setup_logging(self) -> logging.Logger:
        """
        Setup logging for the pipeline.
        
        Returns
        -------
        logging.Logger
            Configured logger instance
        """
        log_dir = os.path.join(self.base_output_dir, 'pipeline_logs')
        os.makedirs(log_dir, exist_ok=True)
        
        # Create logger
        pipeline_logger = logging.getLogger(f'ideal_genom.pipeline.{self.pipeline_name}')
        pipeline_logger.setLevel(logging.DEBUG)
        
        # Avoid duplicate handlers
        if pipeline_logger.handlers:
            return pipeline_logger
        
        # File handler
        log_file = os.path.join(log_dir, f'{self.pipeline_name}.log')
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # Formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)
        
        pipeline_logger.addHandler(file_handler)
        pipeline_logger.addHandler(console_handler)
        
        return pipeline_logger
    
    def get_step_output(self, step_name: str, attribute: str = 'output_path') -> Any:
        """
        Get output from a completed step.
        
        Parameters
        ----------
        step_name : str
            Name of the step
        attribute : str, optional
            Attribute to retrieve (default: 'output_path')
            
        Returns
        -------
        Any
            Value of the requested attribute
            
        Raises
        ------
        ValueError
            If step not found or attribute doesn't exist
        """
        if step_name not in self.steps:
            raise ValueError(f"Step '{step_name}' not found")
        
        step_instance = self.steps[step_name]
        
        if not hasattr(step_instance, attribute):
            raise ValueError(
                f"Step '{step_name}' does not have attribute '{attribute}'"
            )
        
        return getattr(step_instance, attribute)
