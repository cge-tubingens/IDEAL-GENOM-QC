"""
Module to perform principal component analysis to identify ethnicity outliers.
"""

import os
import subprocess
import shutil
import umap
import psutil
import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from ideal_genom_qc.Helpers import shell_do, delete_temp_files
from ideal_genom_qc.get_references import Fetcher1000Genome

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class AncestryQC:
    
    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, dependables_path:str) -> None:

        """
        Initialize the PCA (Principal Component Analysis) object.

        This method initializes the PCA object by checking the paths and the existence of necessary files. It ensures that the input and output paths, as well as the path to dependables, are set upon initialization. Additionally, it verifies the existence of required input files (.bed, .fam, .bim) and reference data from the dependables folder. If any of the necessary files or paths are missing, it raises appropriate errors.

        Parameters:
        -----------
        - input_path (str): Path to the folder containing the input data files.
        - input_name (str): Name of the input data files (without extension).
        - output_path (str): Path to the folder where the output will be saved.
        - output_name (str): Name for the output files.
        - config_dict (str): Dictionary containing configuration settings.
        - dependables_path (str): Path to the folder containing reference data.

        Raises:
        -------
        - ValueError: If input_path, output_path, or dependables_path are not set upon initialization.
        - FileNotFoundError: If any of the required files or paths are missing.

        Attributes:
        -----------
        - input_path (str): Path to the folder containing the input data files.
        - output_path (str): Path to the folder where the output will be saved.
        - input_name (str): Name of the input data files (without extension).
        - output_name (str): Name for the output files.
        - dependables (str): Path to the folder containing reference data.
        - config_dict (str): Dictionary containing configuration settings.
        - dependables_to_keep (list): List of reference data files to keep.
        - results_to_keep (list): List of result files to keep.
        - results_dir (str): Path to the folder where PCA results will be saved.
        - fails_dir (str): Path to the folder where failed samples will be saved.
        - plots_dir (str): Path to the folder where plots will be saved.
        """

        # check if paths are set
        if input_path is None or output_path is None or dependables_path is None:
            raise ValueError("values for input_path, output_path and dependables_path must be set upon initialization.")

        # Check path validity of input data
        bed_path = os.path.join(input_path, input_name + '.bed')
        fam_path = os.path.join(input_path, input_name + '.fam')
        bim_path = os.path.join(input_path, input_name + '.bim')

        bed_check = os.path.exists(bed_path)
        fam_check = os.path.exists(fam_path)
        bim_check = os.path.exists(bim_path)

        if not os.path.exists(input_path) or not os.path.exists(output_path):
            raise FileNotFoundError("input_path or output_path is not a valid path")
        if not bed_check:
            raise FileNotFoundError(".bed file not found")
        if not fam_check:
            raise FileNotFoundError(".fam file not found")
        if not bim_check:
            raise FileNotFoundError(".bim file not found")
        
        # check path validity of reference data
        self.reference_files = {
            'bed': os.path.join(dependables_path, 'all_phase3.bed'),
            'fam': os.path.join(dependables_path, 'all_phase3.fam'),
            'bim': os.path.join(dependables_path, 'all_phase3.bim'),
            'psam': os.path.join(dependables_path, 'all_phase3.psam'),
        }

        ld_region = os.path.join(dependables_path, 'high-LD-regions.txt')

        bed_1000g_check = os.path.exists(self.reference_files['bed'])
        fam_1000g_check = os.path.exists(self.reference_files['fam'])
        bim_1000g_check = os.path.exists(self.reference_files['bim'])
        psam_1000g_check= os.path.exists(self.reference_files['psam'])
        ld_region_check = os.path.exists(ld_region)

        if not os.path.exists(dependables_path):
            raise FileNotFoundError("dependables_path is not a valid path")
        
        logger.info("Checking reference data files")
        
        if not bed_1000g_check or not fam_1000g_check or not bim_1000g_check or not psam_1000g_check:

            logger.info("Reference data files not found")
            logger.info("Downloading reference data")

            reference_fetcher = Fetcher1000Genome()

            reference_fetcher.get_1000genomes()
            logger.info(f"Reference data downloaded successfully to {reference_fetcher.destination}")
            
            reference_fetcher.get_1000genomes_binaries()
            logger.info(f"Reference data binaries generated successfully")

            self.reference_files = {
                'bed': str(reference_fetcher.destination / "all_phase3.bed"),
                'fam': str(reference_fetcher.destination / "all_phase3.fam"),
                'bim': str(reference_fetcher.destination / "all_phase3.bim"),
                'psam': str(reference_fetcher.destination / "all_phase3.psam")
            }

        if not ld_region_check:
            raise FileNotFoundError("high LD regions file not found")

        self.input_path = input_path
        self.output_path= output_path
        self.input_name = input_name
        self.output_name= output_name
        self.dependables= dependables_path

        self.dependables_to_keep = ['all_phase3.bed', 'all_phase3.fam','all_phase3.bim', 'all_phase3.psam', 'high-LD-regions.txt', 'geographic_info.txt']

        self.results_to_keep = ['fail_samples']

        # create results folder
        self.results_dir = os.path.join(output_path, 'ancestry_results')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # create fails folder
        self.fails_dir = os.path.join(self.results_dir, 'fail_samples')
        if not os.path.exists(self.fails_dir):
            os.mkdir(self.fails_dir)

        # create clean files folder
        self.clean_dir = os.path.join(self.results_dir, 'clean_files')
        if not os.path.exists(self.clean_dir):
            os.mkdir(self.clean_dir)
        
        # create figures folder
        self.plots_dir = os.path.join(self.results_dir, 'ancestryQC_plots')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)

    def execute_filter_prob_snps(self)->dict:

        """
        Filter out problematic SNPs (Single Nucleotide Polymorphisms) from the study and reference datasets.

        This method filters SNPs that are non-A-T or non-G-C from the study and reference datasets. It first filters these SNPs from the input dataset and the reference panel separately using the 'filter_non_AT_or_GC_snps' method. Then, it excludes these filtered SNPs from both datasets using PLINK commands, creating new datasets without the problematic SNPs.

        Returns
        -------
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir
        dependables= self.dependables

        reference_panel = 'all_phase3'

        step = "fiter non A-T or G-C snps"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # find A->T and C->G SNPs in study data
        self.filter_non_AT_or_GC_snps(input_dir=input_path, input_name=input_name, results_dir=results_dir)

        # find A->T and C->G SNPs in reference data
        self.filter_non_AT_or_GC_snps(input_dir=dependables, input_name=reference_panel, results_dir=dependables)

        # generate cleaned study data files
        plink_cmd1 = f"plink --bfile  {os.path.join(input_path, input_name)} --chr 1-22 --exclude {os.path.join(results_dir, input_name+'.ac_get_snps')} --threads {max_threads} --make-bed --out {os.path.join(results_dir, input_name+'.no_ac_gt_snps')}"

        # generate cleaned reference data files
        plink_cmd2 = f"plink --bfile  {os.path.join(dependables, reference_panel)} --chr 1-22 --exclude {os.path.join(dependables, reference_panel+'.ac_get_snps')} --allow-extra-chr --memory {memory} --threads {max_threads} --make-bed --out {os.path.join(dependables, reference_panel+'.no_ac_gt_snps')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'reference_out': dependables,
            'study_data': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def execute_ld_pruning(self, ind_pair:list)->dict:

        """
        Prune samples based on Linkage Disequilibrium (LD).

        This method performs LD-based sample pruning using PLINK commands. It filters samples based on Minor Allele Frequency (maf), genotype missingness (geno), individual missingness (mind), and Hardy-Weinberg Equilibrium (hwe). Additionally, it excludes SNPs located in high LD regions specified in the dependables path. The resulting pruned dataset is saved as a new binary file.

        Raises:
        -------
        - TypeError: If maf, geno, mind, or hwe is not of type float.
        - ValueError: If maf, geno, mind, or hwe is not within the specified range.
        - FileNotFoundError: If the file with high LD regions is not found.

        Returns:
        --------
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_name      = self.input_name
        dependables_path= self.dependables
        results_dir     = self.results_dir

        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        if not isinstance(ind_pair[0], int) or not isinstance(ind_pair[1], int):
            raise TypeError("The first two elements in ind_pair values should be integers (windows size and step size)")
        
        if not isinstance(ind_pair[2], float):
            raise TypeError("The third element in ind_pair should be a float (r^2 threshold)")
        
        # check existence of high LD regions file
        high_ld_regions_file = os.path.join(dependables_path, 'high-LD-regions.txt')
        if not os.path.exists(high_ld_regions_file):
            raise FileNotFoundError("File with high LD region was not found")

        step = "ld_prune"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # generates prune.in and prune.out files
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.no_ac_gt_snps')} --exclude {high_ld_regions_file} --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --threads {max_threads} --out {os.path.join(results_dir, input_name)}"

        # prune and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, input_name+'.no_ac_gt_snps')} --extract {os.path.join(results_dir, input_name+'.prune.in')} --threads {max_threads} --make-bed --out {os.path.join(results_dir, input_name+'.pruned')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def execute_ld_prune_ref_panel(self)->dict:

        """
        Prune the reference panel based on the pruned SNPs of the study data.

        This method prunes the reference panel dataset based on the SNPs pruned from the study data during LD-based pruning. It generates a new binary file for the pruned reference panel.

        Returns
        -------
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_name = self.input_name
        dependables= self.dependables
        results_dir= self.results_dir

        step = "prune reference panel"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # generates a pruned reference data files
        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.no_ac_gt_snps')} --extract {os.path.join(results_dir, input_name+'.prune.in')} --make-bed --threads {max_threads} --out {os.path.join(dependables, 'all_phase3.pruned')} --memory {memory}"

        # executes PLINK command
        shell_do(plink_cmd, log=True)
        
        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': dependables
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def execute_fix_chromosome_missmatch(self)->dict:

        """
        Correct chromosome mismatch between study data and reference panel.

        This method corrects any chromosome mismatch between the pruned study data and the pruned reference panel by updating the chromosome information in the reference panel to match that of the study data. It generates a new binary file for the corrected reference panel.

        Returns
        -------
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_name = self.input_name
        dependables= self.dependables
        results_dir= self.results_dir

        step = "chromosome missmatch"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # check that the variant IDs of the reference data have the same chr. ID as the study data
        awk_cmd = f"awk 'BEGIN {{OFS=\"\t\"}} FNR==NR {{a[$2]=$1; next}} ($2 in a && a[$2] != $1) {{print a[$2],$2}}' {os.path.join(results_dir, input_name+'.pruned.bim')} {os.path.join(dependables, 'all_phase3.pruned.bim')} | sed -n '/^[XY]/!p' > {os.path.join(dependables, 'all_phase3.toUpdateChr')}"

        result = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)

        logs = [result.stderr, result.stdout]

        # generates cleaned reference data files
        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.pruned')} --allow-extra-chr --update-chr {os.path.join(dependables, 'all_phase3.toUpdateChr')} 1 2 --threads {max_threads} --make-bed --out {os.path.join(dependables, 'all_phase3.updateChr')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': dependables
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict,
            'awk_logs': logs
        }

        return out_dict

    def execute_fix_position_missmatch_allele_flip(self)->dict:

        """
        Function to handle position mismatch and allele flips.

        Returns
        -------
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_name  = self.input_name
        dependables = self.dependables
        results_dir = self.results_dir

        step = "possition missmatch and allele flips"

        # position missmatch
        awk_cmd1 = f"awk 'BEGIN {{OFS=\"\t\"}} FNR==NR {{a[$2]=$4; next}} ($2 in a && a[$2] != $4)  {{print a[$2],$2}}' {os.path.join(results_dir, input_name+'.pruned.bim')} {os.path.join(dependables, 'all_phase3.pruned.bim')} > {os.path.join(dependables, 'all_phase3.toUpdatePos')}"

        # possible allele flips
        awk_cmd2 = f"awk 'BEGIN {{OFS=\"\t\"}} FNR==NR {{a[$1$2$4]=$5$6; next}} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {{print $2}}' {os.path.join(results_dir, input_name+'.pruned.bim')} {os.path.join(dependables, 'all_phase3.pruned.bim')} > {os.path.join(dependables, 'all_phase3.toFlip')}"

        # executes awk commands
        awks = [awk_cmd1, awk_cmd2]
        logs = []
        for awk in awks:
            result = subprocess.run(awk, shell=True, capture_output=True, text=True)
            logs.append([result.stderr, result.stdout])

        # update positions and flip alleles
        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.updateChr')} --update-map {os.path.join(dependables, 'all_phase3.toUpdatePos')} 1 2 --flip {os.path.join(dependables, 'all_phase3.toFlip')} --make-bed --out {os.path.join(dependables, 'all_phase3.flipped')}"

        # executes PLINK command
        shell_do(plink_cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': dependables,
            'other_files': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict,
            'awk_logs': logs
        }

        return out_dict

    def execute_remove_missmatch(self)->dict:

        """
        Function to remove mismatched alleles after allele flipping.

        Returns:
        --------
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_name  = self.input_name
        dependables = self.dependables
        results_dir = self.results_dir

        step = "remove missmatch"

        # identify alleles that do not match after allele flipping
        awk_cmd = f"awk 'BEGIN {{OFS=\"\t\"}} FNR==NR {{a[$1$2$4]=$5$6; next}} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {{print $2}}' {os.path.join(results_dir, input_name+'.pruned.bim')} {os.path.join(dependables, 'all_phase3.flipped.bim')} > {os.path.join(dependables, 'all_phase3.missmatch')}"

        # executes awk command
        result = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)
        log = [result.stderr, result.stdout]

        # generates cleaned binary files
        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.flipped')} --exclude {os.path.join(dependables, 'all_phase3.missmatch')} --make-bed --out {os.path.join(dependables, 'all_phase3.clean')}"

        # executes PLINK command
        shell_do(plink_cmd, log=True)

        self.dependables_to_keep.append('all_phase3.clean.bed')
        self.dependables_to_keep.append('all_phase3.clean.bim')
        self.dependables_to_keep.append('all_phase3.clean.fam')

        delete_temp_files(self.dependables_to_keep, dependables)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': dependables,
            'other_files': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict,
            'awk_logs': log
        }

        return out_dict
    
    def execute_merge_ref_study(self)->dict:

        """
        Function to merge reference panel with study data.

        Returns:
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_name       = self.input_name
        dependables = self.dependables
        results_dir = self.results_dir

        step = "merge reference panel with study data"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # merge refenrence and study data
        plink_cmd = f"plink --bfile {os.path.join(results_dir, input_name+'.pruned')} --bmerge {os.path.join(dependables, 'all_phase3.clean.bed')} {os.path.join(dependables, 'all_phase3.clean.bim')} {os.path.join(dependables, 'all_phase3.clean.fam')} --make-bed --out {os.path.join(results_dir, input_name+'.merged')} --memory {memory} --threads {max_threads}"

        # executes PLINK command
        shell_do(plink_cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': dependables,
            'other_files': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def execute_pc_decomposition(self, pca:int, maf:float)->dict:

        """
        Run Principal Component Analysis (PCA) on the study data and perform ancestry inference to filter samples based on population outliers.

        This function executes PCA analysis on the study data merged with the reference panel. It then performs ancestry inference to identify population outliers using predefined population tags and a specified threshold. Samples identified as population outliers are removed from the study data.

        Returns:
        -------
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.

        Raises:
        - TypeError: If the PCA parameter is not of type int.
        """

        input_name = self.input_name
        output_name= self.output_name
        results_dir= self.results_dir

        step = "pca_decomposition"

        # check `pca` type
        if not isinstance(pca, int):
            raise TypeError("pca should be an integer value")
        if pca < 1:
            raise ValueError("pca should be a positive integer")
        
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # runs pca analysis
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.merged')} --keep-allele-order --maf {maf} --out {os.path.join(results_dir, output_name+'.pca')} --pca {pca} --memory {memory} --threads {max_threads}"

        # executes PLINK command
        shell_do(plink_cmd1, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def get_ancestry_outliers(self, ref_threshold:float, stu_threshold:float, reference_pop:str, num_pcs:int=2)->dict:

        dependables= self.dependables
        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir
        fails_dir  = self.fails_dir
        output_name= self.output_name

        step = "get ancestry outliers"

        # add population tags to pca output
        df = self.population_tags(
            psam_path     =os.path.join(dependables, 'all_phase3.psam'),
            study_fam_path=os.path.join(input_path, input_name+'.fam')
        )
        df['ID1'] = df['ID1'].astype(str)

        # filter samples who are ethnicity outliers
        ancestry_fails = self.pca_fail(
            df_tags      =df, 
            results_dir  =results_dir,
            output_folder=fails_dir,
            output_name  =output_name, 
            ref_threshold=ref_threshold,
            stu_threshold=stu_threshold,
            reference_pop=reference_pop,
            num_pcs      =num_pcs
        )

        # report
        process_complete = True

        outfiles_dict = {
            'fail_samples': fails_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def execute_drop_ancestry_outliers(self)->dict:

        input_path = self.input_path
        input_name = self.input_name
        output_name= self.output_name
        fails_dir  = self.fails_dir
        clean_dir  = self.clean_dir

        step = "drop ancestry outliers"

        # create cleaned binary files
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --allow-no-sex --remove {os.path.join(fails_dir, output_name+'.fail-ancestry-qc.txt')} --make-bed --out {os.path.join(clean_dir, output_name+'-ancestry-clean')}"

        # execute PLINK command
        shell_do(plink_cmd2, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': clean_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def pca_plot(self)->dict:

        """
        Generate Principal Component Analysis (PCA) plots for the study data.

        This function generates PCA plots based on the results of PCA analysis previously     performed on the study data merged with the reference panel. Two types of plots are     generated: a 2D scatter plot showing the first two principal components and a 3D scatter plot showing the first three principal components. The plots are colored according to the population supergroup (SuperPop) obtained from the population tags.

        Returns:
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_path = self.input_path
        input_name = self.input_name
        output_name= self.output_name
        dependables= self.dependables
        results_dir= self.results_dir

        step = "generate pca plots"

        # add population tags to pca output
        df = self.population_tags(
            psam_path= os.path.join(dependables, 'all_phase3.psam'),
            study_fam_path=os.path.join(input_path, input_name+'.fam')
        )
        df['ID1'] = df['ID1'].astype(str)

        # load .eigenvec file and keep the first three principal components
        df_eigenvec = pd.read_csv(
            os.path.join(results_dir, output_name+'.pca.eigenvec'),
            header=None,
            sep   =r"\s+",
            engine='python'
        )
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:5]].copy()
        df_eigenvec.columns = ['ID1', 'ID2', 'pc_1', 'pc_2', 'pc_3']
        df_eigenvec['ID1'] = df_eigenvec['ID1'].astype(str)

        # merge to get data with tagged populations
        df = pd.merge(df_eigenvec, df, on=['ID1', 'ID2'])

        # generates a 2D scatter plot
        fig, ax = plt.subplots(figsize=(10,10))
        scatter_plot= sns.scatterplot(data=df, x='pc_1', y='pc_2', hue='SuperPop', ax=ax, marker='.', s=70)
        scatter_fig = scatter_plot.get_figure()
        scatter_fig.savefig(os.path.join(self.plots_dir, 'pca.jpeg'), dpi=400)

        # generates a 3D scatter plot
        fig2= plt.figure()
        ax  = fig2.add_subplot(111, projection='3d')

        for s in df['SuperPop'].unique():
            ax.scatter(
                xs=df.pc_1[df.SuperPop==s],
                ys=df.pc_2[df.SuperPop==s],
                zs=df.pc_3[df.SuperPop==s], 
                label=s
            )
        ax.legend()
        plt.savefig(os.path.join(self.plots_dir, 'pca_3d.jpeg'), dpi=400)
        plt.close()

        # report
        process_complete = True

        outfiles_dict = {
            'plots_out': self.plots_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    @staticmethod
    def filter_non_AT_or_GC_snps(input_dir:str, input_name:str, results_dir:str)->list:

        """
        Filters out single nucleotide polymorphisms (SNPs) that do not belong to the categories AT, TA, GC, or CG.

        This method filters SNPs from the input file based on their nucleotide sequences. It selects SNPs with sequences GC, CG, AT, or TA and writes their IDs to a file for further processing.

        Parameters:
        -----------
        - input_dir (str): The directory containing the input files.
        - input_name (str): The name of the input file.
        - results_dir (str): The directory where the results will be stored.

        Returns:
        - list: A list containing stdout and stderr outputs from the execution of the AWK command.
        """

        bim_target = os.path.join(input_dir, input_name+'.bim')
        ac_gt_snps = os.path.join(results_dir, input_name+'.ac_get_snps')

        awk_cmd = f"awk 'BEGIN {{OFS=\"\t\"}} ($5$6 == \"GC\" || $5$6 == \"CG\" || $5$6 == \"AT\" || $5$6 == \"TA\") {{print $2}}' {bim_target} > {ac_gt_snps}"

        result = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)

        logs = [result.stderr, result.stdout]

        return logs
    
    @staticmethod
    def population_tags(psam_path:str, study_fam_path:str)->pd.DataFrame:

        """
        Creates a DataFrame containing population tags by merging information from two input files.

        This method reads population information from a .psam file and individual IDs from a .fam file. It creates a DataFrame containing individual IDs and their corresponding population tags.

        Parameters:
        -----------
        - psam_path (str): The path to the .psam file containing population information.
        - study_fam_path (str): The path to the study .fam file containing individual IDs.

        Returns:
        --------
        - pd.DataFrame: A DataFrame containing population tags with individual IDs.
        """

        # Read population information from the .psam file
        df_psam = pd.read_csv(
            psam_path,
            sep='\t',
            usecols=['#IID', 'SuperPop']
        )

        # Set an ID column and rename columns for consistency
        df_psam['ID'] = 0
        df_psam = df_psam[['ID', '#IID', 'SuperPop']].copy()
        df_psam.columns = ['ID1', 'ID2', 'SuperPop']

        # read individual IDs from the study .fam file
        df_fam = pd.read_csv(
            study_fam_path,
            sep=' ',
            header=None,
            index_col=False
        )

        # select relevant columns, assign a placeholder population tag, and rename columns
        df_fam = df_fam[df_fam.columns[:2]].copy()
        df_fam['SuperPop'] = 'StPop'
        df_fam.columns = ['ID1', 'ID2', 'SuperPop']

        # concatenate the two DataFrames to merge the information
        return pd.concat([df_fam, df_psam], axis=0)

    @staticmethod
    def pca_fail(df_tags:pd.DataFrame, results_dir:str, output_folder:str, output_name:str, ref_threshold:int, stu_threshold:int, reference_pop:str, num_pcs:str=2)->str:

        """
        Identifies samples failing ancestry quality control based on principal component analysis (PCA).

        This method identifies samples failing ancestry quality control based on PCA results.
        It compares the PCA coordinates of study samples with reference population samples.
        Samples with PCA coordinates outside a certain threshold from the reference mean are considered outliers.

        Parameters:
        -----------
        - df_tags (pd.DataFrame): DataFrame containing population tags and PCA information.
        - results_dir (str): Directory path containing PCA results.
        - output_folder (str): Directory path to save output files.
        - output_name (str): Prefix for the output file names.
        - threshold (int): Threshold value for identifying outliers.

        Returns:
        --------
        - str: Path to the file containing the list of samples failing ancestry quality control.
        """

        if not isinstance(ref_threshold, (float, int)):
            raise TypeError("ref_threshold should be an integer or float value")
        if not isinstance(stu_threshold, (float, int)):
            raise TypeError("stu_threshold should be an integer or float value")
        if stu_threshold<=0:
            raise ValueError("stu_threshold should be a positive value")
        if ref_threshold<=0:
            raise ValueError("ref_threshold should be a positive value")
        if not isinstance(reference_pop, str):
            raise TypeError("reference_pop should be a string")
        if not isinstance(num_pcs, int):
            raise TypeError("num_pcs should be an integer value")
        if num_pcs<1:
            raise ValueError("num_pcs should be a positive integer")
        

        # filters South Asian subjects
        mask1 = (df_tags['SuperPop']==reference_pop)
        # filters subjects from study data
        mask2 = (df_tags['SuperPop']=='StPop')

        # generates two data frames with filtered subjects
        df_ref = df_tags[mask1].reset_index(drop=True)
        df_stu = df_tags[mask2].reset_index(drop=True)

        # read .eigenvec file
        df_eigenvec = pd.read_csv(
            os.path.join(results_dir, output_name+'.pca.eigenvec'),
            header=None,
            sep   =r"\s+",
            engine='python'
        )

        if num_pcs>df_eigenvec.shape[1]-2:
            raise ValueError("num_pcs should be less than or equal to the number of principal components in the .eigenvec file")
        
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:2+num_pcs]].copy()

        # renames columns for consistency
        new_col_names = []
        for k in range(2+num_pcs):
            if k<2:
                new_col_names.append(f"ID{k+1}")
            else:
                new_col_names.append(f"pc_{k-1}")
        df_eigenvec.columns = new_col_names

        # merge filtered subjects with its principal components
        df_ref = df_ref.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)
        df_stu = df_stu.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)

        # computes mean and standard deviation by columns in reference data
        mean_ref= df_ref[df_ref.columns[2:]].mean()
        std_ref = df_ref[df_ref.columns[2:]].std()

        # creates empty data frame
        outliers_1 = pd.DataFrame(columns=df_ref.columns)
        outliers_1[df_stu.columns[:2]] = df_stu[df_stu.columns[:2]]

        # identifies subjects with more than `ref_threshold` std deviations from the reference mean
        for col in outliers_1.columns[2:]:
            outliers_1[col] = (np.abs(df_stu[col] - mean_ref[col]) > ref_threshold*std_ref[col])

        outliers_1['is_out'] = (np.sum(outliers_1.iloc[:,2:], axis=1) >0)

        df_1 = outliers_1[outliers_1['is_out']].reset_index(drop=True)[['ID1', 'ID2']].copy()

        # computes mean and standard deviation by columns in study data
        mean_stu= df_stu[df_stu.columns[2:]].mean()
        std_stu = df_stu[df_stu.columns[2:]].std()

        # creates empty data frame
        outliers_2 = pd.DataFrame(columns=df_ref.columns)
        outliers_2[df_stu.columns[:2]] = df_stu[df_stu.columns[:2]]

        # identifies subjects with more than `stu_threshold` std deviation from the study mean
        for col in outliers_2.columns[2:]:
            outliers_2[col] = (np.abs(df_stu[col] - mean_stu[col]) > stu_threshold*std_stu[col])

        outliers_2['is_out'] = (np.sum(outliers_2.iloc[:,2:], axis=1) >0)

        df_2 = outliers_2[outliers_2['is_out']].reset_index(drop=True)[['ID1', 'ID2']].copy()

        df = pd.merge(df_1, df_2, on=['ID1', 'ID2'])

        # save samples considered as ethnicity outliers
        df.to_csv(
            os.path.join(output_folder, output_name+'.fail-ancestry-qc.txt'),
            header=None,
            index =False,
            sep   ='\t'
        )

        return os.path.join(output_folder, output_name+'.fail-ancestry-qc.txt')
