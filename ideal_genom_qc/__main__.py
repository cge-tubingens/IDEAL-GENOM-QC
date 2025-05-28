import os
import json
import logging

from pathlib import Path

from ideal_genom_qc.Helpers import arg_parser

from ideal_genom_qc.SampleQC import SampleQC
from ideal_genom_qc.VariantQC import VariantQC
from ideal_genom_qc.AncestryQC import AncestryQC
from ideal_genom_qc.PopStructure import UMAPplot, FstSummary

from ideal_genom_qc.get_references import FetcherLDRegions
from ideal_genom_qc.check_tools import check_required_tools, get_tool_version, ToolNotFoundError

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def qc_pipeline(params_dict: dict, data_dict: dict, steps_dict: dict, recompute_merge: str, built: str) -> None:

    sample_params     = params_dict['sample_qc']
    ancestry_params   = params_dict['ancestry_qc']
    variant_qc_params = params_dict['variant_qc']
    umap_params       = params_dict['umap_plot']

    input_path = Path(data_dict['input_directory'])
    output_path = Path(data_dict['output_directory'])

    recompute_merge_str = recompute_merge.lower()
    recompute_merge_bool = (recompute_merge_str == 'true')

    high_ld_file = Path(data_dict['high_ld_file'])

    if not high_ld_file.exists() or not high_ld_file.is_file():
        logger.info("LD regions file not found.")
        logger.info("Downloading LD regions file.")

        fetcher = FetcherLDRegions(built=built)
        high_ld_file = fetcher.get_ld_regions()

        logger.info(f"LD regions file downloaded to {high_ld_file}.")


    if steps_dict['sample']:
        # instantiate SampleQC class
        sample_qc = SampleQC(
            input_path      =input_path,
            input_name      =data_dict['input_prefix'],
            output_path     =output_path,
            output_name     =data_dict['output_prefix'],
            high_ld_file    =high_ld_file
        )

        sample_qc_steps = {
            'rename SNPs'           : (sample_qc.execute_rename_snpid, {"rename": sample_params['rename_snp']}),
            'hh_to_missing'         : (sample_qc.execute_haploid_to_missing, {"hh_to_missing": sample_params['hh_to_missing']}),
            'ld_pruning'            : (sample_qc.execute_ld_pruning, {"ind_pair": sample_params['ind_pair']}),
            'miss_genotype'         : (sample_qc.execute_miss_genotype, { "mind": sample_params['mind']}),
            'sex_check'             : (sample_qc.execute_sex_check, {"sex_check": sample_params['sex_check']}),
            'heterozygosity'        : (sample_qc.execute_heterozygosity_rate, {"maf": sample_params['maf']}),
            'duplicates_relatedness': (sample_qc.execute_duplicate_relatedness, {"kingship": sample_params['kingship'], "use_king": sample_params['use_kingship']}),
            'get_fail_samples'      : (sample_qc.get_fail_samples, {"call_rate_thres": sample_params['mind'], "std_deviation_het": sample_params['het_deviation'], "maf_het": sample_params['maf'], "ibd_threshold": sample_params['ibd_threshold']}),
            'drop_fail_samples'     : (sample_qc.execute_drop_samples, {}),
            'clean_input_files'     : (sample_qc.clean_input_folder, {}),
            'clean_results_files'   : (sample_qc.clean_result_folder, {}),
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
            'clean_input_files'     : 'Clean input folder',
            'clean_results_files'   : 'Clean results folder',
        }

        for name, (func, params) in sample_qc_steps.items():
            print(f"\033[1m{step_description[name]}.\033[0m")
            func(**params)

        print("\033[92mSample quality control done.\033[0m")

    # execute step by step
    if steps_dict['ancestry']:

        # instantiate AncestryQC class
        ancestry_qc = AncestryQC(
            input_path = output_path / 'sample_qc_results' / 'clean_files', 
            input_name = data_dict['output_prefix']+'-clean-samples', 
            output_path= output_path, 
            output_name= data_dict['output_prefix'], 
            high_ld_file= high_ld_file,
            recompute_merge=recompute_merge_bool,
            built=built
        )

        ancestry_qc_steps = {
            'merge_study_reference'    : (ancestry_qc.merge_reference_study, {"ind_pair":ancestry_params['ind_pair']}),
            'delete_intermediate_files': (ancestry_qc._clean_merging_dir, {}),
            'pca_analysis'             : (ancestry_qc.run_pca, 
                {
                    "ref_population": ancestry_params['reference_pop'],
                    "pca":ancestry_params['pca'],
                    "maf":ancestry_params['maf'],
                    "num_pca":ancestry_params['num_pcs'],
                    "ref_threshold":ancestry_params['ref_threshold'],
                    "stu_threshold":ancestry_params['stu_threshold'],
                    "aspect_ratio":ancestry_params['aspect_ratio'],
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
            input_path      =output_path / 'ancestry_qc_results' / 'clean_files',
            input_name      =data_dict['output_prefix']+'-ancestry-cleaned',
            output_path     =output_path,
            output_name     =data_dict['output_prefix']
        )

        variant_qc_steps = {
            'Missing data rate'         : (variant_qc.execute_missing_data_rate, {'chr_y': variant_qc_params['chr-y']}),
            'Different genotype'        : (variant_qc.execute_different_genotype_call_rate, {}),
            'Hardy-Weinberg equilibrium': (variant_qc.execute_hwe_test, {}),
            'Get fail variants'         : (variant_qc.get_fail_variants, {'marker_call_rate_thres': variant_qc_params['miss_data_rate'], 'case_controls_thres': variant_qc_params['diff_genotype_rate'], 'hwe_threshold':variant_qc_params['hwe']}),
            'Drop fail variants'        : (variant_qc.execute_drop_variants, {'maf': variant_qc_params['maf'], 'geno': variant_qc_params['geno'], 'hwe': variant_qc_params['hwe']}),
        }

        variant_step_description = {
            'Missing data rate'         : 'Compute missing data rate for males and females',
            'Different genotype'        : 'Case/control nonrandom missingness test',
            'Hardy-Weinberg equilibrium': 'Hardy-Weinberg equilibrium test',
            'Get fail variants'         : 'Get variants that failed quality control',
            'Drop fail variants'        : 'Drop variants that failed quality control'
        }

        for name, (func, params) in variant_qc_steps.items():
            print(f"\033[34m{variant_step_description[name]}.\033[0m")
            func(**params)

        print("\033[92mVariant quality control done.\033[0m")

    if steps_dict['umap']:

        # instantiate umap class
        umap_plots = UMAPplot(
            input_path      =output_path / 'variant_qc_results' / 'clean_files', 
            input_name      =data_dict['output_prefix']+'-variantQCed', 
            output_path     =output_path
        )

        umap_steps = {
            'ld_pruning': (umap_plots.ld_pruning, {
                'maf': umap_params['umap_maf'], 
                'mind': umap_params['umap_mind'], 
                'geno': umap_params['umap_geno'], 
                'hwe': umap_params['umap_hwe'], 
                'ind_pair': umap_params['umap_ind_pair']
            }),
            'comp_pca'  : (umap_plots.compute_pcas, {'pca': umap_params['umap_pca']}),
            'draw_plots': (umap_plots.generate_plots, {
                'color_hue_file': Path(umap_params['color_hue_file']),
                'case_control_markers': umap_params['case_control_marker'],
                'n_neighbors': umap_params['n_neighbors'],
                'metric': umap_params['metric'],
                'min_dist': umap_params['min_dist'],
                'random_state': umap_params['random_state'],
                'umap_kwargs': umap_params['umap_kwargs'],
            })
        }


        umap_step_description = {
            'ld_pruning': 'LD pruning',
            'comp_pca'  : 'Compute PCAs',
            'draw_plots': 'Generate UMAP plots'
        }

        for name, (func, params) in umap_steps.items():
            print(f"\033[34m{umap_step_description[name]}.\033[0m")
            func(**params)

        print("\033[92mUMAP plots done.\033[0m")

    if steps_dict['fst']:

        # instantiate umap class
        fst = FstSummary(
            input_path      =output_path / 'variant_qc_results' / 'clean_files', 
            input_name      =data_dict['output_prefix']+'-variantQCed', 
            output_path     =output_path,
            built           =built
        )

        fst_steps = {
            'merge_study_reference': (fst.merge_reference_study, {"ind_pair":ancestry_params['ind_pair']}),
            'add_population_tags'  : (fst.add_population_tags, {}),
            'compute_fst'          : (fst.compute_fst, {}),
            'report_fst'           : (fst.report_fst, {}),
        }

        fst_step_description = {
            'merge_study_reference': 'Merge cleaned data with reference panel',
            'add_population_tags'  : 'Add population tags to the data',
            'compute_fst'          : 'Compute Fst statistics',
            'report_fst'           : 'Generate Fst report'
        }

        for name, (func, params) in fst_steps.items():
            print(f"\033[34m{fst_step_description[name]}.\033[0m")
            func(**params)

        print("\033[92mFst computation is done.\033[0m") 

    pass

def main()->str:

    required = ['plink', 'plink2']

    try:
        check_required_tools(required)
        for tool in required:
            version = get_tool_version(tool)
            logger.info(f"{tool} version: {version}")
    except ToolNotFoundError as e:
        logger.error(e)
        
    args_dict = arg_parser()
    #args_dict = vars(args)

    params_path = args_dict['path_params']
    data_path   = args_dict['file_folders']
    steps_path  = args_dict['steps']
    recompute_merge = args_dict['recompute_merge'].lower()
    built      = args_dict['build']

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
        recompute_merge=recompute_merge,
        built    =built
    )

    return "Pipeline is completed"

if __name__ == "__main__":
    main()
