import os
import json
import pandas as pd

from ideal_genom_qc.Helpers import arg_parser

from ideal_genom_qc.SampleQC import SampleQC
from ideal_genom_qc.VariantQC import VariantQC
from ideal_genom_qc.AncestryQC import AncestryQC
from ideal_genom_qc.UMAPplot import UMAPplot

def qc_pipeline(params_dict:dict, data_dict:dict, steps_dict:dict, use_kingship:str)->None:

    sample_params = params_dict['sample_qc']
    ancestry_params = params_dict['ancestry_qc']
    variant_qc_params = params_dict['variant_qc']
    umap_params = params_dict['umap_plot']

    use_kingship = use_kingship.lower()
    if use_kingship == 'true':
        use_kingship = True
    else:
        use_kingship = False

    print('use kingship', type(use_kingship))

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
            'ld_pruning'            : (sample_qc.execute_ld_pruning, (sample_params['ind_pair'],)),
            'miss_genotype'         : (sample_qc.execute_miss_genotype, (sample_params['mind'],)),
            'sex_check'             : (sample_qc.execute_sex_check, (sample_params['sex_check'],)),
            'heterozygosity'        : (sample_qc.execute_heterozygosity_rate, (sample_params['maf'],)),
            'duplicates_relatedness': (sample_qc.execute_duplicate_relatedness, (sample_params['kingship'], use_kingship,)),
            'get_fail_samples'      : (sample_qc.get_fail_samples, (sample_params['mind'], sample_params['het_deviation'], sample_params['maf'], sample_params['ibd_threshold'],)),
            'drop_fail_samples'     : (sample_qc.execute_drop_samples, ()),
            'recover_SNPs_names'    : (sample_qc.execute_recover_snp_names, (True,),)
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
            'drop_fail_samples'     : 'Drop samples that failed quality control',
            'recover_SNPs_names'    : 'Recover SNPs names'
        }

        for name, (func, params) in sample_qc_steps.items():
            print(f"\033[34m{step_description[name]}.\033[0m")
            func(*params)

        print("\033[92mSample quality control done.\033[0m")

    # execute step by step
    if steps_dict['ancestry']:

        # instantiate PCA class 
        ancestry_qc = AncestryQC(
            input_path      =os.path.join(data_dict['output_directory'], 'sample_qc_results', 'clean_files'),
            input_name      =data_dict['output_prefix']+'-clean-samples',
            output_path     =data_dict['output_directory'],
            output_name     =data_dict['output_prefix'],
            dependables_path=data_dict['dependables_directory']
        )

        ancestry_qc_steps = {
            'filter_snps'              : (ancestry_qc.execute_filter_prob_snps, ()),
            'LD_pruning'               : (ancestry_qc.execute_ld_pruning, (ancestry_params['ind_pair'],)),
            'reference_pruning'        : (ancestry_qc.execute_ld_prune_ref_panel, ()),
            'chr_missmatch'            : (ancestry_qc.execute_fix_chromosome_missmatch,()),
            'pos_missmatch_allele_flip': (ancestry_qc.execute_fix_position_missmatch_allele_flip, ()),
            'remove_missmatch'         : (ancestry_qc.execute_remove_missmatch, ()),
            'merging'                  : (ancestry_qc.execute_merge_ref_study, ()),
            'pca_analysis'             : (ancestry_qc.execute_pc_decomposition, (ancestry_params['pca'], ancestry_params['maf'],)),
            'get_outliers'             : (ancestry_qc.get_ancestry_outliers, (ancestry_params['ref_threshold'], ancestry_params['stu_threshold'], 'SAS', ancestry_params['num_pca'],)),
            'pca_plot'                 : (ancestry_qc.pca_plot, ()),
            'drop_fail_samples'        : (ancestry_qc.execute_drop_ancestry_outliers, ())
        }

        step_description = {
            'filter_snps'              : 'filter problematic snps',
            'LD_pruning'               : 'LD prune study_population',
            'reference_pruning'        : 'LD prune reference panel',
            'chr_missmatch'            : 'fix any chromosome missmatch',
            'pos_missmatch_allele_flip': 'fix any issue with position missmatch and allele flip',
            'remove_missmatch'         : 'remove missmatchesd SNPs',
            'merging'                  : 'merge reference panel and study population',
            'pca_analysis'             : 'execute PC decomposition',
            'get_outliers'             : 'find samples with discordant ancestry',
            'pca_plot'                 : 'generate PCA plot',
            'drop_fail_samples'        : 'drop samples with discordant ancestry'
        }

        for name, (func, params) in ancestry_qc_steps.items():
            print(f"\033[34m{step_description[name]}.\033[0m")
            func(*params)

        print("\033[92mAncestry outliers analysis done.\033[0m")

    if steps_dict['variant']:
        variant_qc = VariantQC(
            input_path      =os.path.join(data_dict['output_directory'], 'ancestry_results', 'clean_files'),
            input_name      =data_dict['output_prefix']+'-ancestry-clean',
            output_path     =data_dict['output_directory'],
            output_name     =data_dict['output_prefix']
        )

        variant_qc_steps = {
            'Missing data rate'         : (variant_qc.execute_missing_data_rate, (variant_qc_params['chr-y'],)),
            'Different genotype'        : (variant_qc.execute_different_genotype_call_rate, ()),
            'Get fail variants'         : (variant_qc.get_fail_variants, (variant_qc_params['miss_data_rate'], variant_qc_params['diff_genotype_rate'],)),
            'Drop fail variants'        : (variant_qc.execute_drop_variants, (variant_qc_params['geno'], variant_qc_params['hwe'], variant_qc_params['maf'],)),
        }

        step_description = {
            'Missing data rate'         : 'Compute missing data rate for males and females',
            'Different genotype'        : 'Case/control nonrandom missingness test',
            'Get fail variants'         : 'Get variants that failed quality control',
            'Drop fail variants'        : 'Drop variants that failed quality control'
        }

        for name, (func, params) in variant_qc_steps.items():
            print(f"\033[34m{step_description[name]}.\033[0m")
            func(*params)

        print("\033[92mVariant quality control done.\033[0m")

    if steps_dict['umap']:

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

    with open(params_path, 'r') as file:
        params_dict = json.load(file)

    with open(steps_path, 'r') as file:
        steps_dict = json.load(file)

    qc_pipeline(
        params_dict =params_dict,
        data_dict   =data_dict,
        steps_dict  =steps_dict,
        use_kingship=use_kingship
    )

    return "Pipeline is completed"

if __name__ == "__main__":
    execute_main()
