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

    output = res.stdout.decode('utf-8') + res.stderr.decode('utf-8')

    if log:
        print(output)
    if return_log:
        return output
    if err:
        return res.stderr.decode('utf-8')
    
def arg_parser()->dict:

    # define parser
    parser = argparse.ArgumentParser(description='Adresses to configuration files')

    # parameters of quality control
    parser.add_argument('--path_params', type=str, nargs='?', default=None, const=None, help='Full path to the JSON file containing genotype quality control parameters.')

    # path to data and names of files
    parser.add_argument('--file_folders', type=str, nargs='?', default=None, const=None, help='Full path to the JSON file containing folder names and locations for genotype quality control data.')

    # parse args and turn into dict
    args = parser.parse_args()

    return args

def delete_temp_files(files_to_keep:list, path_to_folder:str)->None:

    for file in os.listdir(path_to_folder):
            file_split = file.split('.')
            if file_split[-1]!='log' and file not in files_to_keep:
                os.remove(
                    os.path.join(path_to_folder, file)
                )
        
    # create log folder for dependables
    logs_dir = os.path.join(path_to_folder, 'log_files')
    os.mkdir(logs_dir)

    for file in os.listdir(path_to_folder):
        if file.split('.')[-1]=='log':
            shutil.move(
                os.path.join(path_to_folder, file),
                os.path.join(logs_dir, file)
            )

        pass
