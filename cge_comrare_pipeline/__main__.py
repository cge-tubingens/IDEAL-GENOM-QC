import os
import json
import pandas as pd

from cge_comrare_pipeline.Helpers import arg_parser

def execute_main():

    args = arg_parser()
    args_dict = vars(args)

    from cge_comrare_pipeline.SampleQC import SampleQC
    from cge_comrare_pipeline.VariantQC import VariantQC

    params_path = args_dict['path_params']
    data_path = args_dict['file_folders']

    # check path to config files
    if not os.path.exists(data_path):
        raise FileNotFoundError("Configuration file with path to data and analysis results cannot be found.")
    
    if not os.path.exists(params_path):
        raise FileNotFoundError("Configuration file with pipeline parameters cannot be found.")

    # open config file
    with open(data_path, 'r') as file:
        data_dict = json.load(file)

    if params_path is None:
        params_dict = {
                "maf" : 0.05,
                "geno": 0.1,
                "mind": 0.1,
                "hwe" : 0.00000005,
                "sex_check": [0.2, 0.8],
                "indep-pairwise": [50, 5, 0.2],
                "chr": 24,
                "outlier_threshold": 6,
                "pca": 10
                }
    else:
        with open(params_path, 'r') as file:
            params_dict = json.loads(file)

    # class instances
    sample_qc = SampleQC(
        input_path      =data_dict['input_directory'],
        input_name      =data_dict['input_prefix'],
        output_path     =data_dict['output_directory'],
        output_name     =data_dict['output_prefix'],
        config_dict     =params_dict,
        dependables_path=data_dict['dependables_directory']
    )

    variant_qc = VariantQC(
        input_path      =data_dict['input_directory'],
        input_name      =data_dict['input_prefix'],
        output_path     =data_dict['output_directory'],
        output_name     =data_dict['output_prefix'],
        config_path     =params_path,
        dependables_path=data_dict['dependables_directory']
    )

    # pipeline steps
    steps = {
        'ld_prune'      : sample_qc.run_ld_prune,
        'heterozygosity': sample_qc.run_heterozygosity_rate,
        'sex_check'     : sample_qc.run_sex_check,
        'relatedness'   : sample_qc.run_relatedness_prune,
        'delete_samples': sample_qc.delete_failing_QC,
        'ancestry_one'  : sample_qc.divergent_ancestry_step_one,
        'run_pca'       : sample_qc.run_pca_analysis,
        'miss_data'     : variant_qc.missing_data_rate,
        'call_rate'     : variant_qc.different_genotype_call_rate,
        'delete_markers': variant_qc.remove_markers
    }

    # execute pipeline
    for step in steps.keys():
        steps[step]()

    return None

if __name__ == "__main__":
    execute_main()
