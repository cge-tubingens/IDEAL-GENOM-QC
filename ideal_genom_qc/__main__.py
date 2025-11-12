import json
import logging

from pathlib import Path

from ideal_genom_qc.Helpers import arg_parser, validate_config

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

        # Execute the complete sample QC pipeline
        sample_qc.execute_sample_qc_pipeline(sample_params)
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

        ancestry_qc.execute_ancestry_pipeline(ancestry_params)

        print("\033[92mAncestry outliers analysis done.\033[0m")

    if steps_dict['variant']:
        variant_qc = VariantQC(
            input_path      =output_path / 'ancestry_qc_results' / 'clean_files',
            input_name      =data_dict['output_prefix']+'-ancestry-cleaned',
            output_path     =output_path,
            output_name     =data_dict['output_prefix']
        )

        variant_qc.execute_variant_qc_pipeline(variant_qc_params)

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

    # check CLI config arguments
    validate_config(args_dict['file_folders'], args_dict['path_params'], args_dict['steps'], args_dict['build'])

    params_path = args_dict['path_params']
    data_path   = args_dict['file_folders']
    steps_path  = args_dict['steps']
    recompute_merge = args_dict['recompute_merge'].lower()
    built      = args_dict['build']

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

    return "✅ Pipeline is completed"

if __name__ == "__main__":
    main()
