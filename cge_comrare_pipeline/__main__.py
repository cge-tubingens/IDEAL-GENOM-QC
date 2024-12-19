import os
import json
import pandas as pd

from cge_comrare_pipeline.Helpers import arg_parser

from cge_comrare_pipeline.SampleQC import SampleQC
from cge_comrare_pipeline.VariantQC import VariantQC
from cge_comrare_pipeline.AncestryQC import AncestryQC
from cge_comrare_pipeline.UMAPplot import UMAPplot

def pipe_SamOutVar(params_dict:dict, data_dict:dict, steps_dict:dict, use_kingship:str)->None:

    if steps_dict['sample']:
        # instantiate SampleQC class
        sample_qc = SampleQC(
            input_path      =data_dict['input_directory'],
            input_name      =data_dict['input_prefix'],
            output_path     =data_dict['output_directory'],
            output_name     =data_dict['output_prefix'],
            dependables_path=data_dict['dependables_directory'],
        )

        # pipeline steps
        sample_qc_steps = {
            'rename SNPs'           : (sample_qc.execute_rename_snps, (True,)),
            'hh_to_missing'         : (sample_qc.execute_haploid_to_missing, ()),
            'ld_pruning'            : (sample_qc.execute_ld_pruning, (sample_params['ind_par'],)),
            'miss_genotype'         : (sample_qc.execute_miss_genotype, (sample_params['mind'],)),
            'sex_check'             : (sample_qc.execute_sex_check, (sample_params['sex_check'],)),
            'heterozygosity'        : (sample_qc.execute_heterozygosity_rate, (sample_params['maf'],)),
            'duplicates_relatedness': (sample_qc.execute_duplicate_relatedness, (sample_params['kingship'], use_kingship,)),
            'get_fail_samples'      : (sample_qc.get_fail_samples, (sample_params['mind'], sample_params['het_deviation'], sample_params['maf'], sample_params['ibd_threshold'],)),
            'drop_fail_samples'     : (sample_qc.execute_drop_samples, ())
        }

        step_description = {
            'rename SNPs'           : 'Rename SNPs to chr:pos:ref:alt',
            'hh_to_missing'         : 'Solve hh warnings by setting to missing',
            'ld_pruning'            : 'Perform LD pruning',
            'miss_genotype'         : 'Get samples with high missing rate',
            'sex_check'             : 'Get samples with discordant sex information',
            'heterozygosity'        : 'Get samples with high heterozygosity rate',
            'duplicates_relatedness': 'Get samples with high relatedness rate or duplicates',
            'get_fail_samples'      : 'Get samples that failed quality control',
            'drop_fail_samples'     : 'Drop samples that failed quality control'
        }

        for name, (func, params) in sample_qc_steps.items():
            print(f"\033[1m{step_description[name]}.\033[0m")
            func(*params)

        print("\033[92mSample quality control done.\033[0m")

    # execute step by step
    if steps_dict['pca']:

        # instantiate PCA class 
        pca_qc = PCA(
            input_path      =os.path.join(data_dict['output_directory'], 'sample_qc_results'),
            input_name      =data_dict['output_prefix']+'.clean',
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
            'umap_plot'                : pca_qc.reference_umap_plot,
            'pca_plot'                 : pca_qc.pca_plot
        }

        for step in pca_steps.keys():
            print(step)
            pca_steps[step]()

        print("\033[92mEthnicity outliers analysis done.\033[0m")

    if steps_dict['variant']:
        variant_qc = VariantQC(
            input_path      =os.path.join(data_dict['output_directory'], 'pca_results'),
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

        print("\033[92mVariant quality control done.\033[0m")

    if steps_dict['umap_plots']:

        # instantiate umap class
        umap_plots = UMAPplot(
            input_path      =os.path.join(data_dict['output_directory'], 'variant_qc_results'), 
            input_name      =data_dict['output_prefix']+'.vrnt_clean', 
            dependables_path=data_dict['dependables_directory'],
            config_dict     =params_dict,
            output_path     =data_dict['output_directory']
        )
        umap_steps = {
            'ld_pruning': umap_plots.ld_pruning,
            'comp_pca'  : umap_plots.compute_pcas,
            'draw_plots': umap_plots.generate_plots
        }

        for step in umap_steps.keys():
            print(step)
            umap_steps[step]()

        print("\033[92mUMAP plots done.\033[0m")       

    pass

def execute_main()->str:

    args = arg_parser()
    args_dict = vars(args)

    params_path = args_dict['path_params']
    data_path   = args_dict['file_folders']
    steps_path  = args_dict['steps']
    use_kingship= args_dict['use_kingship'].lower()

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

    pipe_SamOutVar(
        params_dict =params_dict,
        data_dict   =data_dict,
        steps_dict  =steps_dict,
        use_kingship=use_kingship
    )

    return "Pipeline is completed"

if __name__ == "__main__":
    execute_main()
