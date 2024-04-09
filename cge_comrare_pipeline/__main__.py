import os
import json
import pandas as pd

from cge_comrare_pipeline.Helpers import arg_parser

def execute_main()->None:

    args = arg_parser()
    args_dict = vars(args)

    from cge_comrare_pipeline.SampleQC import SampleQC
    from cge_comrare_pipeline.VariantQC import VariantQC
    from cge_comrare_pipeline.PCA import PCA

    params_path = args_dict['path_params']
    data_path = args_dict['file_folders']
    steps_path = args_dict['steps']

    # check path to config files
    if not os.path.exists(data_path):
        raise FileNotFoundError("Configuration file with path to data and analysis results cannot be found.")
    
    if not os.path.exists(params_path):
        raise FileNotFoundError("Configuration file with pipeline parameters cannot be found.")
    
    if not os.path.exists(steps_path):
        raise FileNotFoundError("Configuration file with pipeline steps cannot be found.")

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
            params_dict = json.load(file)

    with open(steps_path, 'r') as file:
        steps_dict = json.load(file)

    # execute step by step
    if steps_dict['pca']:

        # instantiate PCA class 
        pca_qc = PCA(
            input_path      =data_dict['input_directory'],
            input_name      =data_dict['input_prefix'],
            output_path     =data_dict['output_directory'],
            output_name     =data_dict['output_prefix'],
            config_dict     =params_dict,
            dependables_path=data_dict['dependables_directory']
        )

        # pipeline steps
        pca_steps = {
            'filter_snps'              : pca_qc.filter_problematic_snps,
            'LD_pruning'               : pca_qc.ld_pruning,
            'reference_pruning'        : pca_qc.prune_reference_panel,
            'chr_missmatch'            : pca_qc.chromosome_missmatch,
            'pos_missmatch_allele_flip': pca_qc.position_missmatch_allele_flip,
            'remove_missmatch'         : pca_qc.remove_missmatch,
            'merging'                  : pca_qc.merge_with_reference,
            'pca_analysis'             : pca_qc.run_pca_analysis,
            'pca_plot'                 : pca_qc.pca_plot
        }

        for step in pca_steps.keys():
            pca_steps[step]()

    if steps_dict['sample']:
        # instantiate SampleQC class
        sample_qc = SampleQC(
            input_path      =os.path.join(data_dict['output_directory'], 'pca_results'),
            input_name      =data_dict['output_prefix']+'.clean',
            output_path     =data_dict['output_directory'],
            output_name     =data_dict['output_prefix'],
            config_dict     =params_dict,
            dependables_path=data_dict['dependables_directory']
        )

        # pipeline steps
        smpl_steps = {
        'sex_check'     : sample_qc.run_sex_check,
        'heterozygosity': sample_qc.run_heterozygosity_rate,
        'relatedness'   : sample_qc.run_relatedness_prune,
        'delete_samples': sample_qc.delete_failing_QC,
        }

        for step in smpl_steps.keys():
            smpl_steps[step]()

    if steps_dict['variant']:
        variant_qc = VariantQC(
            input_path      =os.path.join(data_dict['output_directory'], 'sample_qc_results'),
            input_name      =data_dict['output_prefix']+'.clean',
            output_path     =data_dict['output_directory'],
            output_name     =data_dict['output_prefix'],
            config_dict     =params_dict,
            dependables_path=data_dict['dependables_directory']
        )

        vrnt_steps = {
            'miss_data'     : variant_qc.missing_data_rate,
            'call_rate'     : variant_qc.different_genotype_call_rate,
            'delete_markers': variant_qc.remove_markers
        }

        for step in vrnt_steps.keys():
            vrnt_steps[step]()

    return None

if __name__ == "__main__":
    execute_main()
