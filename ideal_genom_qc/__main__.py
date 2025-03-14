import os
import json

from pathlib import Path

from ideal_genom_qc.Helpers import arg_parser

from ideal_genom_qc.SampleQC import SampleQC
from ideal_genom_qc.VariantQC import VariantQC
from ideal_genom_qc.AncestryQC import AncestryQC
from ideal_genom_qc.UMAPplot import UMAPplot

def qc_pipeline(params_dict: dict, data_dict: dict, steps_dict: dict, use_kingship: str, recompute_merge: str)->None:

    sample_params     = params_dict['sample_qc']
    ancestry_params   = params_dict['ancestry_qc']
    variant_qc_params = params_dict['variant_qc']
    umap_params       = params_dict['umap_plot']

    use_kingship = use_kingship.lower()
    if use_kingship == 'true':
        use_kingship = True
    else:
        use_kingship = False

    recompute_merge = recompute_merge.lower()
    if recompute_merge == 'true':
        recompute_merge = True
    else:
        recompute_merge = False

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

        sample_step_description = {
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
            print(f"\033[34m{sample_step_description[name]}.\033[0m")
            func(*params)

        print("\033[92mSample quality control done.\033[0m")

    # execute step by step
    if steps_dict['ancestry']:

        # instantiate AncestryQC class
        ancestry_qc = AncestryQC(
            input_path = Path(data_dict['output_directory']) / 'sample_qc_results' / 'clean_files', 
            input_name = data_dict['output_prefix']+'-clean-samples', 
            output_path= Path(data_dict['output_directory']), 
            output_name= data_dict['output_prefix'], 
            high_ld_regions= Path(data_dict['dependables_directory']) / 'high-LD-regions.txt',
            recompute_merge=recompute_merge
        )

        ancestry_qc_steps = {
            'merge_study_reference'    : (ancestry_qc.merge_reference_study, {"ind_pair":ancestry_params['indep']}),
            'delete_intermediate_files': (ancestry_qc._clean_merging_dir, {}),
            'pca_analysis'             : (ancestry_qc.run_pca, 
                {
                    "ref_population": ancestry_params['reference_pop'],
                    "pca":ancestry_params['pca'],
                    "maf":ancestry_params['maf'],
                    "num_pca":ancestry_params['num_pcs'],
                    "ref_threshold":ancestry_params['ref_threshold'],
                    "stu_threshold":ancestry_params['stu_threshold'],
                }
            ),
        }

        step_description = {
            'merge_study_reference'    : "Merge reference genome with study genome",
            'delete_intermediate_files': "Delete intermediate files generated during merging",
            'pca_analysis'             : "Run a PCA analysis to perfom ancestry QC"
        }

        for name, (func, params) in ancestry_qc_steps.items():
            print(f"\033[1m{step_description[name]}.\033[0m")
            func(**params)

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
            'Drop fail variants'        : (variant_qc.execute_drop_variants, (variant_qc_params['maf'], variant_qc_params['geno'], variant_qc_params['hwe'],)),
        }

        variant_step_description = {
            'Missing data rate'         : 'Compute missing data rate for males and females',
            'Different genotype'        : 'Case/control nonrandom missingness test',
            'Get fail variants'         : 'Get variants that failed quality control',
            'Drop fail variants'        : 'Drop variants that failed quality control'
        }

        for name, (func, params) in variant_qc_steps.items():
            print(f"\033[34m{variant_step_description[name]}.\033[0m")
            func(*params)

        print("\033[92mVariant quality control done.\033[0m")

    if steps_dict['umap']:

        # instantiate umap class
        umap_plots = UMAPplot(
            input_path      =os.path.join(data_dict['output_directory'], 'variant_qc_results', 'clean_files'), 
            input_name      =data_dict['output_prefix']+'-variantQCed', 
            dependables_path=data_dict['dependables_directory'],
            output_path     =data_dict['output_directory']
        )

        umap_steps = {
            'ld_pruning': (umap_plots.ld_pruning, (umap_params['maf'], umap_params['geno'], umap_params['mind'], umap_params['hwe'], umap_params['ind_pair'],)),
            'comp_pca'  : (umap_plots.compute_pcas, (umap_params['pca'],)),
            'draw_plots': (umap_plots.generate_plots, (umap_params['n_neighbors'], umap_params['min_dist'], umap_params['metric'], ))
        }

        umap_step_description = {
            'ld_pruning': 'LD pruning',
            'comp_pca'  : 'Compute PCAs',
            'draw_plots': 'Generate UMAP plots'
        }

        for name, (func, params) in umap_steps.items():
            print(f"\033[34m{umap_step_description[name]}.\033[0m")
            func(*params)

        print("\033[92mUMAP plots done.\033[0m")       

    pass

def execute_main()->str:

    args = arg_parser()
    args_dict = vars(args)

    params_path = args_dict['path_params']
    data_path   = args_dict['file_folders']
    steps_path  = args_dict['steps']
    use_kingship= args_dict['use_kingship'].lower()
    recompute_merge = args_dict['recompute_merge'].lower()
    built      = args_dict['built']

    # check path to config files
    if not os.path.exists(data_path):
        raise FileNotFoundError("Configuration file with path to data and analysis results cannot be found.")
    
    if not os.path.exists(params_path):
        raise FileNotFoundError("Configuration file with pipeline parameters cannot be found.")
    
    if not os.path.exists(steps_path):
        raise FileNotFoundError("Configuration file with pipeline steps cannot be found.")
    
    if built not in ['37', '38']:
        raise ValueError("Built of the human genome must be 37 or 38.")

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
        use_kingship=use_kingship,
        recompute_merge=recompute_merge
    )

    return "Pipeline is completed"

if __name__ == "__main__":
    execute_main()
