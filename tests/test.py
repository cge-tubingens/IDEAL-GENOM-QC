from cge_comrare_pipeline.SampleQC import SampleQC
from cge_comrare_pipeline.VariantQC import VariantQC
from cge_comrare_pipeline.PCA import PCA

import pandas as pd

import json
import os


INPUT_PATH = '/mnt/0A2AAC152AABFBB7/PipeLine/data/inputData'
INPUT_NAME = 'test_1'
OUTPUT_PATH= '/mnt/0A2AAC152AABFBB7/PipeLine/data/outputData'
OUTPUT_NAME= 'results'
CONFIG_PATH= '/mnt/0A2AAC152AABFBB7/comrare-pipeline/config.JSON'
DEPEND_PATH= '/mnt/0A2AAC152AABFBB7/PipeLine/data/auxiliarData'

data_path = '/mnt/0A2AAC152AABFBB7/PipeLine/data/config_files/parameters.JSON'
with open(data_path, 'r') as file:
        parameters = json.load(file)

pca = PCA(
    input_path=INPUT_PATH,
    input_name=INPUT_NAME,
    output_path=OUTPUT_PATH,
    output_name=OUTPUT_NAME,
    config_dict=parameters,
    dependables_path=DEPEND_PATH
)

pca.ld_pruning()

pca.run_pca_analysis()

pca.plot_pca(label='ourG')

sample_QC = SampleQC(
    input_path=INPUT_PATH,
    input_name=INPUT_NAME,
    output_path=OUTPUT_PATH,
    output_name=OUTPUT_NAME,
    config_dict=parameters,
    dependables_path=DEPEND_PATH
)

sample_QC.run_heterozygosity_rate()

sample_QC.run_sex_check()

sample_QC.run_relatedness_prune()

sample_QC.delete_failing_QC()

variant_QC = VariantQC(
    input_path=os.path.join(OUTPUT_PATH, 'clean_samples'),
    input_name=OUTPUT_NAME+'.smpl_clean',
    output_path=OUTPUT_PATH,
    output_name=OUTPUT_NAME,
    config_dict=parameters,
    dependables_path=DEPEND_PATH
)

variant_QC.missing_data_rate()

variant_QC.different_genotype_call_rate()

variant_QC.remove_markers()
