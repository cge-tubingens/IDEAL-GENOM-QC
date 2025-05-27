import subprocess
import sys
import argparse
import os
import shutil

def shell_do(command, print_cmd=False, log=False, return_log=False, err=False):
    """
    From GenoTools
    """
    
    if print_cmd:
        print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res = subprocess.run(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Check if the command failed
    if res.returncode != 0:
        error_message = f"Command failed with return code {res.returncode}: {res.stderr.decode('utf-8')}"
        print(error_message, file=sys.stderr)
        raise RuntimeError(error_message)

    output = res.stdout.decode('utf-8') + res.stderr.decode('utf-8')

    if log:
        print(output)
    if return_log:
        return output
    if err:
        return res.stderr.decode('utf-8')
    
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

def delete_temp_files(files_to_keep: list, path_to_folder: str) -> None:

    """
    Function to delete temporary files that were created during the pipeline execution. Moreover, it creates a directory called 'log_files' to save al `.log` files originated from the pipeline execution.

    Parameters
    ----------
    files_to_keep: list
        list of strings where its elements are the names of files and folders that should be kept.
    path_to_folder: str
        full path to the folder where the temporary files are located.
    """

    for file in os.listdir(path_to_folder):
        file_split = file.split('.')
        if file_split[-1]!='log' and file not in files_to_keep and file_split[-1]!='hh':
            if os.path.isfile(os.path.join(path_to_folder, file)):
                os.remove(
                    os.path.join(path_to_folder, file)
                )
        
    # create log folder for dependables
    logs_dir = os.path.join(path_to_folder, 'log_files')
    if not os.path.exists(logs_dir):
        os.mkdir(logs_dir)

    for file in os.listdir(path_to_folder):
        if file.split('.')[-1]=='log' or file.split('.')[-1]=='hh':
            shutil.move(
                os.path.join(path_to_folder, file),
                os.path.join(logs_dir, file)
            )

    pass
