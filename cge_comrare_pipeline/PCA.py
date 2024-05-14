"""
Module to perform principal component analysis to identify ethnicity outliers.
"""

import os
import subprocess
import shutil
import umap

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from cge_comrare_pipeline.Helpers import shell_do, delete_temp_files

class PCA:
    
    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_dict:str, dependables_path:str) -> None:

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
        bed_1000g = os.path.join(dependables_path, 'all_phase3.bed')
        fam_1000g = os.path.join(dependables_path, 'all_phase3.fam')
        bim_1000g = os.path.join(dependables_path, 'all_phase3.bim')
        psam_1000g= os.path.join(dependables_path, 'all_phase3.psam')
        ld_region = os.path.join(dependables_path, 'high-LD-regions.txt')

        bed_1000g_check = os.path.exists(bed_1000g)
        fam_1000g_check = os.path.exists(fam_1000g)
        bim_1000g_check = os.path.exists(bim_1000g)
        psam_1000g_check= os.path.exists(psam_1000g)
        ld_region_check = os.path.exists(ld_region)

        if not os.path.exists(dependables_path):
            raise FileNotFoundError("dependables_path is not a valid path")
        if not bed_1000g_check:
            raise FileNotFoundError("all_phase3.bed file not found")
        if not fam_1000g_check:
            raise FileNotFoundError("all_phase3.fam file not found")
        if not bim_1000g_check:
            raise FileNotFoundError("all_phase3.bim file not found")
        if not psam_1000g_check:
            raise FileNotFoundError("all_phase3.psam file not found")
        if not ld_region_check:
            raise FileNotFoundError("high LD regions file not found")

        self.input_path = input_path
        self.output_path= output_path
        self.input_name = input_name
        self.output_name= output_name
        self.dependables= dependables_path
        self.config_dict= config_dict

        self.dependables_to_keep = ['all_phase3.bed', 'all_phase3.fam','all_phase3.bim', 'all_phase3.psam', 'high-LD-regions.txt']

        self.results_to_keep = ['fail_samples']

        # create results folder
        self.results_dir = os.path.join(output_path, 'pca_results')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # create fails folder
        self.fails_dir = os.path.join(self.results_dir, 'fail_samples')
        if not os.path.exists(self.fails_dir):
            os.mkdir(self.fails_dir)
        
        # create figures folder
        self.plots_dir = os.path.join(output_path, 'plots')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)

    def shorten_variant_id(self)->dict:

        """
        Function to deal with long variant IDs. It will be done at a later stage.
        """

        input_path      = self.input_path
        input_name      = self.input_name
        dependables_path= self.dependables

        reference_panel = 'all_phase3'

        step = "shorten length of variant IDs"

        awk_cmd1 = f"awk < {os.path.join(input_path, input_name+'.bim')} '{{print $1\":\"$4, $2}}' > {os.path.join(input_path, input_name+'.names')}"

        awk_cmd2 = f"awk < {os.path.join(input_path, input_name+'.bim')} '{{$2=$1\":\"$4;print $0}}' > {os.path.join(input_path, input_name+'_0.bim')}"

        awk_cmd3 = f"awk < {os.path.join(dependables_path, reference_panel+'.bim')} '{{print $1\":\"$4, $2}}' > {os.path.join(dependables_path, reference_panel+'.names')}"

        awk_cmd4 = f"awk < {os.path.join(dependables_path, reference_panel+'.bim')} '{{$2=$1\":\"$4;print $0}}' > {os.path.join(dependables_path, reference_panel+'_0.bim')}"

        shutil.copy(
            os.path.join(input_path, input_name+'.bed'), 
            os.path.join(input_path, input_name+'_0.bed')
        )
        shutil.copy(
            os.path.join(input_path, input_name+'.fam'), 
            os.path.join(input_path, input_name+'_0.fam')
        )
        shutil.copy(
            os.path.join(dependables_path, reference_panel+'.bed'), 
            os.path.join(dependables_path, reference_panel+'_0.bed')
        )
        shutil.copy(
            os.path.join(dependables_path, reference_panel+'.fam'), 
            os.path.join(dependables_path, reference_panel+'_0.fam')
        )

        logs = []
        cmds = [awk_cmd1, awk_cmd2, awk_cmd3, awk_cmd4]
        for cmd in cmds:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            logs.append([result.stderr, result.stdout])
        
        # report
        process_complete = True

        outfiles_dict = {
            'output': input_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def filter_problematic_snps(self)->dict:

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

        # find A->T and C->G SNPs in study data
        self.filter_non_AT_or_GC_snps(input_dir=input_path, input_name=input_name, results_dir=results_dir)

        # find A->T and C->G SNPs in reference data
        self.filter_non_AT_or_GC_snps(input_dir=dependables, input_name=reference_panel, results_dir=dependables)

        # generate cleaned study data files
        plink_cmd1 = f"plink --bfile  {os.path.join(input_path, input_name)} --chr 1-22 --exclude {os.path.join(results_dir, input_name+'.ac_get_snps')} --make-bed --out {os.path.join(results_dir, input_name+'.no_ac_gt_snps')}"

        # generate cleaned reference data files
        plink_cmd2 = f"plink --bfile  {os.path.join(dependables, reference_panel)} --chr 1-22 --exclude {os.path.join(dependables, reference_panel+'.ac_get_snps')} --allow-extra-chr --memory 10240 --make-bed --out {os.path.join(dependables, reference_panel+'.no_ac_gt_snps')}"

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

    def ld_pruning(self)->dict:

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

        maf      = self.config_dict['maf']
        geno     = self.config_dict['geno']
        mind     = self.config_dict['mind']
        hwe      = self.config_dict['hwe']
        ind_pair = self.config_dict['indep-pairwise']

        # Check type of maf
        if not isinstance(maf, float):
             raise TypeError("maf should be of type float.")

        # Check type of geno
        if not isinstance(geno, float):
            raise TypeError("geno should be of type float.")

        # Check type of mind
        if not isinstance(mind, float):
            raise TypeError("mind should be of type float.")
        
        # Check type of hwe
        if not isinstance(hwe, float):
            raise TypeError("hwe should be of type float.")
        
        # Check if maf is in range
        if maf < 0.05 or maf > 0.1:
            raise ValueError("maf should be between 0.05 and 0.1")
        
        # Check if geno is in range
        if geno < 0.05 or geno > 0.1:
            raise ValueError("geno should be between 0.05 and 0.1")
        
        # Check if mind is in range
        if mind < 0.1 or mind > 0.15:
            raise ValueError("mind should be between 0.1 and 0.15")
        
        # Check if hwe is in range
        if hwe < 0.00000001 or hwe > 0.001:
            raise ValueError("hwe should be between 0.00000001 and 0.001")
        
        # check existence of high LD regions file
        high_ld_regions_file = os.path.join(dependables_path, 'high-LD-regions.txt')
        if not os.path.exists(high_ld_regions_file):
            raise FileNotFoundError("File with high LD region was not found")

        step = "ld_prune"

        # generates prune.in and prune.out files
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.no_ac_gt_snps')} --maf {maf} --geno {geno} --mind {mind} --hwe {hwe} --exclude {high_ld_regions_file} --range --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --out {os.path.join(results_dir, input_name)}"

        # prune and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, input_name+'.no_ac_gt_snps')} --keep-allele-order --extract {os.path.join(results_dir, input_name+'.prune.in')} --make-bed --out {os.path.join(results_dir, input_name+'.pruned')}"

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
    
    def prune_reference_panel(self)->dict:

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

        # generates a pruned reference data files
        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.no_ac_gt_snps')} --keep-allele-order --allow-extra-chr --extract {os.path.join(results_dir, input_name+'.prune.in')} --make-bed --out {os.path.join(dependables, 'all_phase3.pruned')}"

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

    def chromosome_missmatch(self)->dict:

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

        # check that the variant IDs of the reference data have the same chr. ID as the study data
        awk_cmd = f"awk 'BEGIN {{OFS=\"\t\"}} FNR==NR {{a[$2]=$1; next}} ($2 in a && a[$2] != $1) {{print a[$2],$2}}' {os.path.join(results_dir, input_name+'.pruned.bim')} {os.path.join(dependables, 'all_phase3.pruned.bim')} | sed -n '/^[XY]/!p' > {os.path.join(dependables, 'all_phase3.toUpdateChr')}"

        result = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)

        logs = [result.stderr, result.stdout]

        # generates cleaned reference data files
        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.pruned')} --allow-extra-chr --update-chr {os.path.join(dependables, 'all_phase3.toUpdateChr')} 1 2 --make-bed --out {os.path.join(dependables, 'all_phase3.updateChr')}"

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

    def position_missmatch_allele_flip(self)->dict:

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

    def remove_missmatch(self)->dict:

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
    
    def merge_with_reference(self)->dict:

        """
        Function to merge reference panel with study data.

        Returns:
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_name       = self.input_name
        dependables = self.dependables
        results_dir = self.results_dir

        step = "merge reference panel with study data"

        # merge refenrence and study data
        plink_cmd = f"plink --bfile {os.path.join(results_dir, input_name+'.pruned')} --bmerge {os.path.join(dependables, 'all_phase3.clean.bed')} {os.path.join(dependables, 'all_phase3.clean.bim')} {os.path.join(dependables, 'all_phase3.clean.fam')} --make-bed --out {os.path.join(results_dir, input_name+'.merged')}"

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
    
    def run_pca_analysis(self)->dict:

        """
        Run Principal Component Analysis (PCA) on the study data and perform ancestry inference to filter samples based on population outliers.

        This function executes PCA analysis on the study data merged with the reference panel. It then performs ancestry inference to identify population outliers using predefined population tags and a specified threshold. Samples identified as population outliers are removed from the study data.

        Returns:
        -------
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.

        Raises:
        - TypeError: If the PCA parameter is not of type int.
        """

        input_path = self.input_path
        input_name = self.input_name
        output_name= self.output_name
        results_dir= self.results_dir
        fails_dir  = self.fails_dir
        dependables= self.dependables
        ref_threshold= self.config_dict['ref_threshold']
        stu_threshold= self.config_dict['stu_threshold']
        reference_pop= self.config_dict['reference_pop']
        pca        = self.config_dict['pca']

        step = "pca_analysis"

        # check `pca` type
        if not isinstance(pca, int):
            raise TypeError("pca should be an integer value")

        # runs pca analysis
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.merged')} --keep-allele-order --maf 0.01 --out {os.path.join(results_dir, output_name+'.pca')} --pca {pca}"

        # executes PLINK command
        shell_do(plink_cmd1, log=True)

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
            reference_pop=reference_pop
        )

        # create cleaned binary files
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --allow-no-sex --remove {ancestry_fails} --make-bed --out {os.path.join(self.results_dir, output_name+'.clean')}"

        self.results_to_keep.append(output_name+'.clean.bed')
        self.results_to_keep.append(output_name+'.clean.bim')
        self.results_to_keep.append(output_name+'.clean.fam')
        self.results_to_keep.append(output_name+'.pca.eigenvec')

        shell_do(plink_cmd2, log=True)

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
    
    def reference_umap_plot(self)->dict:

        input_path = self.input_path
        input_name = self.input_name
        output_name= self.output_name
        dependables= self.dependables
        results_dir= self.results_dir

        n_neighbors = self.config_dict["umap_n_neighbors"]
        min_dist    = self.config_dict["umap_min_dist"]
        metric      = self.config_dict["umap_metric"]
        pca         = self.config_dict['pca']

        step = "generate umap plots"

        # add population tags to pca output
        df = self.population_tags(
            psam_path= os.path.join(dependables, 'all_phase3.psam'),
            study_fam_path=os.path.join(input_path, input_name+'.fam')
        )
        df['ID1'] = df['ID1'].astype(str)

        # load .eigenvec file
        df_eigenvec = pd.read_csv(
            os.path.join(results_dir, output_name+'.pca.eigenvec'),
            header=None,
            sep=' '
        )
        new_cols = [f"pca_{k}" for k in range(1,pca+1)]
        df_eigenvec.columns = ['ID1', 'ID2'] + new_cols

        df_ids = df_eigenvec[['ID1', 'ID2']].copy()
        df_vals= df_eigenvec[new_cols].to_numpy()

        del df_eigenvec

        D2_redux = umap.UMAP(
            n_components=2,
            n_neighbors =n_neighbors,
            min_dist    =min_dist,
            metric      =metric
        )

        D3_redux = umap.UMAP(
            n_components=3,
            n_neighbors =n_neighbors,
            min_dist    =min_dist,
            metric      =metric
        )

        umap_2D_proj = D2_redux.fit_transform(df_vals)
        umap_3D_proj = D3_redux.fit_transform(df_vals)

        df_2D = pd.concat([df_ids, pd.DataFrame(data=umap_2D_proj, columns=['proj_1', 'proj_2'])], axis=1)
        df_3D = pd.concat([df_ids, pd.DataFrame(data=umap_3D_proj, columns=['proj_1', 'proj_2', 'proj_3'])], axis=1)
        
        df_2D = pd.merge(df_2D, df, on=['ID1', 'ID2'])
        df_3D = pd.merge(df_3D, df, on=['ID1', 'ID2'])

        # generates a 2D scatter plot
        fig, ax = plt.subplots(figsize=(10,10))
        scatter_plot= sns.scatterplot(
            data=df_2D, 
            x='proj_1', 
            y='proj_2', 
            hue='SuperPop',
            marker='.',
            alpha=0.6,
            ax=ax
        )
        scatter_fig = scatter_plot.get_figure()
        scatter_fig.savefig(os.path.join(self.plots_dir, 'umap_2d.pdf'))

        # generates a 3D scatter plot
        fig2= plt.figure(figsize=(10,10))
        ax  = fig2.add_subplot(111, projection='3d')

        for s in df['SuperPop'].unique():
            ax.scatter(
                xs=df_3D.proj_1[df_3D.SuperPop==s],
                ys=df_3D.proj_2[df_3D.SuperPop==s],
                zs=df_3D.proj_3[df_3D.SuperPop==s], 
                label=s,
                marker='.'
            )
        ax.legend()
        plt.savefig(os.path.join(self.plots_dir, 'umap_3d.pdf'), format='pdf')
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
            sep=' '
        )
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:5]].copy()
        df_eigenvec.columns = ['ID1', 'ID2', 'pc_1', 'pc_2', 'pc_3']
        df_eigenvec['ID1'] = df_eigenvec['ID1'].astype(str)

        # merge to get data with tagged populations
        df = pd.merge(df_eigenvec, df, on=['ID1', 'ID2'])

        # generates a 2D scatter plot
        fig, ax = plt.subplots(figsize=(10,10))
        scatter_plot= sns.scatterplot(data=df, x='pc_1', y='pc_2', hue='SuperPop', ax=ax, marker='.')
        scatter_fig = scatter_plot.get_figure()
        scatter_fig.savefig(os.path.join(self.plots_dir, 'pca.pdf'))

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
        plt.savefig(os.path.join(self.plots_dir, 'pca_3d.pdf'), format='pdf')
        plt.close()

        # delete temporary files
        # delete_temp_files(self.results_to_keep, results_dir)

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
    def pca_fail(df_tags:pd.DataFrame, results_dir:str, output_folder:str, output_name:str, ref_threshold:int, stu_threshold:int, reference_pop:str)->str:

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
            sep=' '
        )

        # renames columns for consistency
        new_col_names = []
        for k in range(df_eigenvec.shape[1]):
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
            index=False,
            sep=' '
        )

        return os.path.join(output_folder, output_name+'.fail-ancestry-qc.txt')
