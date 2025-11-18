import subprocess
import sys
import argparse
import os
    
def arg_parser() -> dict:

    # define parser
    parser = argparse.ArgumentParser(description='Addresses to configuration files')

    ## Required config paths
    parser.add_argument('--path-params', type=str, required=True, help='Path to JSON file with genotype QC parameters.')
    parser.add_argument('--file-folders', type=str, required=True, help='Path to JSON file with folder names and locations for QC data.')
    parser.add_argument('--steps', type=str, required=True, help='Path to JSON file specifying pipeline steps.')

    parser.add_argument('--recompute-merge', type=str, nargs='?', default=True, const=None, help='boolean that determines if the merge of the reference data and study data must be recomputed.')

    parser.add_argument('--build', type=str, default='38', choices=['37', '38'], help='Genome build to use (37 or 38).')

    # parse args and turn into dict
    args = parser.parse_args()

    return vars(args)

def validate_config(data_path, params_path, steps_path, built):
    """Validate config files and genome build."""
    
    if not os.path.exists(data_path):
        raise FileNotFoundError(
            "Configuration file with path to data and analysis results cannot be found."
        )

    if not os.path.exists(params_path):
        raise FileNotFoundError(
            "Configuration file with pipeline parameters cannot be found."
        )

    if not os.path.exists(steps_path):
        raise FileNotFoundError(
            "Configuration file with pipeline steps cannot be found."
        )

    if built not in ['37', '38']:
        raise ValueError("Built of the human genome must be 37 or 38.")

    print("✅ Configuration files and genome build validated successfully.")
