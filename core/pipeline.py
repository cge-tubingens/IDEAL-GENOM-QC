"""Pipeline orchestration engine for ideal_genom."""

import os
import re
import importlib
import logging
from pathlib import Path

from typing import Dict, Any, List

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
    
    def __init__(self, config: Dict[str, Any], dry_run: bool = False):
        """
        Initialize pipeline executor.
        
        Parameters
        ----------
        config : dict
            Pipeline configuration dictionary (from config.load_config)
        dry_run : bool
            If True, skip directory creation and actual execution
        """
        self.config = config
        self.steps = {}  # Store instantiated sub-pipeline objects
        self.base_output_dir = config['pipeline']['base_output_dir']
        self.pipeline_name = config['pipeline']['name']
        self.dry_run = dry_run
        
        # Create base output directory (skip if dry run)
        if not dry_run:
            os.makedirs(self.base_output_dir, exist_ok=True)
        
        # Setup logging
        self.logger = self._setup_logging()
        
    def execute(self) -> None:
        """Execute all pipeline steps sequentially."""
        pipeline_steps = self.config['pipeline']['steps']
        
        # Filter enabled steps and validate dependencies
        enabled_steps = self._filter_enabled_steps(pipeline_steps)
        self._validate_step_dependencies(enabled_steps)
        
        self.logger.info(f"Starting pipeline: {self.pipeline_name}")
        self.logger.info(f"Output directory: {self.base_output_dir}")
        self.logger.info(f"Total steps: {len(pipeline_steps)}")
        self.logger.info(f"Enabled steps: {len(enabled_steps)}")
        
        # First pass: Instantiate ALL classes (enabled and disabled) for reference resolution
        self.logger.info(f"\nInstantiating all pipeline classes for reference resolution...")
        for step_config in pipeline_steps:
            step_name = step_config['name']
            is_enabled = step_config.get('enabled', True)
            
            if not is_enabled:
                self.logger.info(f"Instantiating disabled step for references: {step_name}")
                try:
                    self._instantiate_step(step_config)
                except Exception as e:
                    self.logger.warning(f"⚠️  Failed to instantiate disabled step '{step_name}': {e}")
                    # Continue anyway - this step won't be available for references
        
        # Second pass: Execute only enabled steps
        for i, step_config in enumerate(enabled_steps, 1):
            step_name = step_config['name']
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"Step {i}/{len(enabled_steps)}: {step_name}")
            self.logger.info(f"{'='*60}")
            
            try:
                # If already instantiated (from first pass), just execute
                if step_name in self.steps:
                    self._execute_existing_step(step_config)
                else:
                    # Instantiate and execute
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
        
        # Convert string paths to Path objects for parameters ending with '_path'
        init_params = self._convert_paths_to_path_objects(init_params)
        execute_params = self._convert_paths_to_path_objects(execute_params)
        
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
        
        # Perform cleanup if configured
        self._perform_cleanup(step_name, pipeline_instance)
        
        self.logger.info(f"Step output stored as: steps.{step_name}")
    
    def _instantiate_step(self, step_config: Dict[str, Any]) -> None:
        """
        Instantiate a pipeline step class without executing it.
        Used for disabled steps to enable reference resolution.
        
        Parameters
        ----------
        step_config : dict
            Configuration for the step to instantiate
        """
        step_name = step_config['name']
        
        # Import the class dynamically
        module_path = step_config['module']
        class_name = step_config['class']
        
        self.logger.debug(f"Loading {class_name} from {module_path} (instantiate only)")
        
        try:
            module = importlib.import_module(module_path)
            pipeline_class = getattr(module, class_name)
        except (ImportError, AttributeError) as e:
            raise ImportError(
                f"Failed to import {class_name} from {module_path}: {str(e)}"
            )
        
        # Resolve parameters (handle ${...} references)
        init_params = self._resolve_params(step_config['init_params'])
        
        # Convert string paths to Path objects for parameters ending with '_path'
        init_params = self._convert_paths_to_path_objects(init_params)
        
        self.logger.debug(f"Instantiating {class_name} (disabled step)")
        self.logger.debug(f"Init params: {init_params}")
        
        # Instantiate the sub-pipeline class
        pipeline_instance = pipeline_class(**init_params)
        
        # Store the instance for reference by subsequent steps
        self.steps[step_name] = pipeline_instance
        
        self.logger.debug(f"Disabled step instantiated and stored as: steps.{step_name}")
    
    def _execute_existing_step(self, step_config: Dict[str, Any]) -> None:
        """
        Execute a step that has already been instantiated.
        
        Parameters
        ----------
        step_config : dict
            Configuration for the step to execute
        """
        step_name = step_config['name']
        
        # Get the already instantiated pipeline object
        pipeline_instance = self.steps[step_name]
        
        # Resolve execute parameters
        execute_params = self._resolve_params(step_config.get('execute_params', {}))
        execute_params = self._convert_paths_to_path_objects(execute_params)
        
        # Determine the execute method name
        execute_method_name = f'execute_{step_name}_pipeline'
        
        if not hasattr(pipeline_instance, execute_method_name):
            raise AttributeError(
                f"{pipeline_instance.__class__.__name__} does not have method '{execute_method_name}'"
            )
        
        execute_method = getattr(pipeline_instance, execute_method_name)
        
        self.logger.info(f"Executing {execute_method_name}")
        self.logger.debug(f"Execute params: {execute_params}")
        execute_method(execute_params)
        
        # Perform cleanup if configured
        self._perform_cleanup(step_name, pipeline_instance)
    
    def _perform_cleanup(self, step_name: str, pipeline_instance: Any) -> None:
        """
        Perform cleanup of intermediate files for Sample QC and Variant QC steps.
        
        Parameters
        ----------
        step_name : str
            Name of the pipeline step
        pipeline_instance : Any
            Instance of the pipeline class that was executed
        """
        # Check if cleanup is disabled globally
        keep_intermediate = self.config.get('settings', {}).get('files', {}).get('keep_intermediate', True)
        
        self.logger.info(f"🔍 Cleanup check for {step_name}:")
        self.logger.info(f"   - keep_intermediate setting: {keep_intermediate}")
        self.logger.info(f"   - Config path: settings.files.keep_intermediate")
        
        if keep_intermediate:
            self.logger.info(f"⏭️  Skipping cleanup for {step_name} - keep_intermediate is True")
            return
        
        # Only handle sample_qc and variant_qc steps
        if step_name not in ['sample_qc', 'variant_qc']:
            self.logger.info(f"⏭️  No cleanup configured for step: {step_name}")
            return
        
        self.logger.info(f"🧹 Initiating cleanup for {step_name}...")
        
        try:
            if step_name == 'sample_qc':
                self._cleanup_sample_qc(pipeline_instance)
            elif step_name == 'variant_qc':
                self._cleanup_variant_qc(pipeline_instance)
                
            self.logger.info(f"✅ Cleanup completed for {step_name}")
            
        except Exception as e:
            self.logger.warning(f"⚠️  Cleanup failed for {step_name}: {e}")
            # Don't fail the pipeline for cleanup issues
    
    def _cleanup_sample_qc(self, pipeline_instance: Any) -> None:
        """Cleanup intermediate files for Sample QC step."""
        from ideal_genom.qc.sample_qc import SampleQCCleanUp
        
        # Get required paths from pipeline instance
        if not hasattr(pipeline_instance, 'output_path') or not hasattr(pipeline_instance, 'input_path'):
            self.logger.warning("❌ SampleQC instance missing required paths for cleanup")
            self.logger.warning(f"   Available attributes: {[attr for attr in dir(pipeline_instance) if not attr.startswith('_')]}")
            return
        
        output_path = pipeline_instance.output_path
        input_path = pipeline_instance.input_path
        
        self.logger.info(f"🧹 Running Sample QC cleanup...")
        self.logger.info(f"   - Output path: {output_path}")
        self.logger.info(f"   - Input path: {input_path}")
        
        cleanup = SampleQCCleanUp(
            output_path=output_path,
            input_path=input_path
        )
        cleanup.clean_all()
    
    def _cleanup_variant_qc(self, pipeline_instance: Any) -> None:
        """Cleanup intermediate files for Variant QC step.""" 
        from ideal_genom.qc.variant_qc import VariantQCCleanUp
        
        # Get required paths from pipeline instance
        if not hasattr(pipeline_instance, 'output_path'):
            self.logger.warning("❌ VariantQC instance missing output_path for cleanup")
            self.logger.warning(f"   Available attributes: {[attr for attr in dir(pipeline_instance) if not attr.startswith('_')]}")
            return
        
        output_path = pipeline_instance.output_path
        
        self.logger.info(f"🧹 Running Variant QC cleanup...")
        self.logger.info(f"   - Output path: {output_path}")
        
        cleanup = VariantQCCleanUp(output_path=output_path)
        cleanup.clean_results_files()
    
    
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
    
    def _convert_paths_to_path_objects(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Convert string parameters ending with '_path' or '_file' to Path objects.
        
        Parameters
        ----------
        params : dict
            Parameters dictionary potentially containing path strings
            
        Returns
        -------
        dict
            Parameters with path strings converted to Path objects
        """
        converted = {}
        
        for key, value in params.items():
            if (key.endswith('_path') or key.endswith('_file')) and isinstance(value, str):
                converted[key] = Path(value)
            else:
                converted[key] = value
        
        return converted
    
    def _filter_enabled_steps(self, pipeline_steps: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Filter pipeline steps to include only enabled ones.
        
        Parameters
        ----------
        pipeline_steps : list
            List of all pipeline step configurations
            
        Returns
        -------
        list
            List of enabled pipeline steps
        """
        enabled_steps = []
        
        for step in pipeline_steps:
            # Default to enabled if not specified
            if step.get('enabled', True):
                enabled_steps.append(step)
            else:
                self.logger.info(f"Skipping disabled step: {step['name']}")
        
        return enabled_steps
    
    def _validate_step_dependencies(self, enabled_steps: List[Dict[str, Any]]) -> None:
        """
        Validate step dependencies and issue warnings for potential issues.
        
        Parameters
        ----------
        enabled_steps : list
            List of enabled pipeline steps
        """
        enabled_step_names = {step['name'] for step in enabled_steps}
        all_step_names = {step['name'] for step in self.config['pipeline']['steps']}
        
        # Define step dependencies for genomic QC workflow
        dependencies = {
            'variant_qc': ['sample_qc'],
            'ancestry_qc': ['sample_qc'],  
            'population_analysis': []  # Can run standalone or after any QC step
        }
        
        # Check for missing dependencies and issue warnings
        for step in enabled_steps:
            step_name = step['name']
            
            if step_name in dependencies:
                missing_deps = []
                unavailable_deps = []
                
                for dep in dependencies[step_name]:
                    if dep not in enabled_step_names:
                        missing_deps.append(dep)
                        if dep not in all_step_names:
                            unavailable_deps.append(dep)
                
                if missing_deps:
                    if unavailable_deps:
                        self.logger.warning(
                            f"⚠️  Step '{step_name}' is enabled but dependency "
                            f"step(s) {unavailable_deps} are not defined in configuration."
                        )
                    else:
                        self.logger.info(
                            f"ℹ️  Step '{step_name}' is enabled but dependency "
                            f"step(s) {missing_deps} are disabled. "
                            f"Disabled steps will be instantiated for reference resolution."
                        )
        
        # Check for proper step ordering
        self._validate_step_ordering(enabled_steps)
    
    def _validate_step_ordering(self, enabled_steps: List[Dict[str, Any]]) -> None:
        """
        Validate that steps are in proper execution order.
        
        Parameters
        ----------
        enabled_steps : list
            List of enabled pipeline steps
        """
        # Define preferred step order
        preferred_order = ['sample_qc', 'ancestry_qc', 'variant_qc', 'population_analysis']
        
        # Get positions of enabled steps in preferred order
        step_positions = {}
        for step in enabled_steps:
            step_name = step['name']
            if step_name in preferred_order:
                step_positions[step_name] = preferred_order.index(step_name)
        
        # Check if steps are in correct order
        if len(step_positions) > 1:
            current_order = [step['name'] for step in enabled_steps if step['name'] in preferred_order]
            sorted_order = sorted(current_order, key=lambda x: preferred_order.index(x))
            
            if current_order != sorted_order:
                self.logger.warning(
                    f"⚠️  Steps may be out of optimal order. "
                    f"Current: {current_order}, Recommended: {sorted_order}. "
                    f"This may cause issues with step dependencies."
                )

    def _setup_logging(self) -> logging.Logger:
        """
        Setup logging for the pipeline.
        
        Returns
        -------
        logging.Logger
            Configured logger instance
        """
        # Create logger
        pipeline_logger = logging.getLogger(f'ideal_genom.pipeline.{self.pipeline_name}')
        pipeline_logger.setLevel(logging.DEBUG)
        
        # Avoid duplicate handlers
        if pipeline_logger.handlers:
            return pipeline_logger
        
        # Console handler (always available)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # File handler (only if not dry run)
        if not self.dry_run:
            log_dir = os.path.join(self.base_output_dir, 'pipeline_logs')
            os.makedirs(log_dir, exist_ok=True)
            log_file = os.path.join(log_dir, f'{self.pipeline_name}.log')
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(logging.DEBUG)
        # Formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # Add file handler only if not dry run
        if not self.dry_run:
            file_handler.setFormatter(formatter)
            pipeline_logger.addHandler(file_handler)
        
        console_handler.setFormatter(formatter)
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
    
    def get_pipeline_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the pipeline configuration and status.
        
        Returns
        -------
        dict
            Pipeline summary including enabled steps, dependencies, and configuration
        """
        pipeline_steps = self.config['pipeline']['steps']
        enabled_steps = self._filter_enabled_steps(pipeline_steps)
        
        summary = {
            'pipeline_name': self.pipeline_name,
            'base_output_dir': self.base_output_dir,
            'total_steps': len(pipeline_steps),
            'enabled_steps': len(enabled_steps),
            'steps': []
        }
        
        for step in enabled_steps:
            step_info = {
                'name': step['name'],
                'module': step['module'],
                'class': step['class'],
                'enabled': step.get('enabled', True)
            }
            summary['steps'].append(step_info)
        
        return summary
