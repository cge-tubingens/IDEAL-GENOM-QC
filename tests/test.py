from classes.SampleQC import SampleQC
from classes.VariantQC import VariantQC

INPUT_PATH = '/mnt/0A2AAC152AABFBB7/PipeLine/data/inputData'
INPUT_NAME = 'test_1'
OUTPUT_PATH= '/mnt/0A2AAC152AABFBB7/PipeLine/data/outputData'
OUTPUT_NAME= 'results'
CONFIG_PATH= '/mnt/0A2AAC152AABFBB7/comrare-pipeline/config.JSON'
DEPEND_PATH= '/mnt/0A2AAC152AABFBB7/PipeLine/data/auxiliarData'

sample_QC = SampleQC(
    input_path=INPUT_PATH,
    input_name=INPUT_NAME,
    output_path=OUTPUT_PATH,
    output_name=OUTPUT_NAME,
    config_path=CONFIG_PATH,
    dependables_path=DEPEND_PATH
)

sample_QC.run_ld_prune(ld_region_file='high-LD-regions.txt')

sample_QC.run_heterozygosity_rate()

sample_QC.run_sex_check()

sample_QC.run_relatedness_prune()

sample_QC.delete_failing_QC()

sample_QC.divergent_ancestry_step_one(ld_region_file='high-LD-regions.txt')

sample_QC.run_pca_analysis()

variant_QC = VariantQC(
    input_path=INPUT_PATH,
    input_name=INPUT_NAME,
    output_path=OUTPUT_PATH,
    output_name=OUTPUT_NAME,
    config_path=CONFIG_PATH,
    dependables_path=DEPEND_PATH
)

variant_QC.missing_data_rate()

variant_QC.different_genotype_call_rate()

variant_QC.remove_markers()
