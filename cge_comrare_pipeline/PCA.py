"""
Python module to perform PCA analysis and outliers removal
"""

import os
import json

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from cge_comrare_pipeline.Helpers import shell_do

class PCA:
    
    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_dict:str, dependables_path:str) -> None:

        # check if paths are set
        if input_path is None or output_path is None or dependables_path is None:
            raise ValueError("values for input_path, output_path and dependables_path must be set upon initialization.")

        # Check path validity
        bed_path = os.path.join(input_path, input_name + '.bed')
        fam_path = os.path.join(input_path, input_name + '.fam')
        bim_path = os.path.join(input_path, input_name + '.bim')

        bed_check = os.path.exists(bed_path)
        fam_check = os.path.exists(fam_path)
        bim_check = os.path.exists(bim_path)

        if not os.path.exists(input_path) or not os.path.exists(output_path):
            raise FileNotFoundError("input_path or output_path is not a valid path")
        if not os.path.exists(dependables_path):
            raise FileNotFoundError("dependables_path is not a valid path")
        if not bed_check:
            raise FileNotFoundError(".bed file not found")
        if not fam_check:
            raise FileNotFoundError(".fam file not found")
        if not bim_check:
            raise FileNotFoundError(".bim file not found")

        self.input_path     = input_path
        self.output_path    = output_path
        self.input_name     = input_name
        self.output_name    = output_name
        self.dependables    = dependables_path

        self.config_dict = config_dict

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

        pass

    def ld_pruning(self)->dict:

        """
        Funtion to prunes samples based on Linkage Disequilibrium

        Parameters:
        - ld_region_file: string
            file name with regions with high Linkage Distribution

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        input_path       = self.input_path
        input_name       = self.input_name
        output_path      = self.output_path
        output_name      = self.output_name
        dependables_path = self.dependables
        results_dir = self.results_dir

        maf = self.config_dict['maf']
        geno= self.config_dict['geno']
        mind= self.config_dict['mind']
        hwe = self.config_dict['hwe']
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

        # generates prune.in and prune.out
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --maf {maf} --geno {geno} --mind {mind} --hwe {hwe} --exclude {high_ld_regions_file} --range --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --out {os.path.join(results_dir, output_name)}"

        # prune and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --extract {os.path.join(results_dir, output_name+'.prune.in')} --make-bed --out {os.path.join(results_dir, output_name+'.pruned')}"

        # execute Plink commands
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
    
    def run_pca_analysis(self)->dict:

        """
        Function to prunes samples .....

        Parameters:
        - ld_region_file: string
            file name with regions with high Linkage Distribution

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        input_path = self.input_path
        input_name = self.input_name
        output_name= self.output_name
        results_dir= self.results_dir
        fails_dir  = self.fails_dir
        dependables= self.dependables
        threshold  = self.config_dict['outlier_threshold']
        pca        = self.config_dict['pca']
        maf        = self.config_dict['maf']
        mind       = self.config_dict['mind']

        step = "pca_analysis"

        # check `pca` type
        if not isinstance(pca, int):
            raise TypeError("pca should be an integer value")
        
        # merge target data with 1000genome data
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --autosome --allow-no-sex --maf {maf} --bmerge {os.path.join(dependables, 'clean.bed')} {os.path.join(dependables, 'clean.bim')} {os.path.join(dependables, 'clean.fam')} --extract {os.path.join(results_dir, output_name+'.prune.in')} --make-bed --out {os.path.join(results_dir, 'merged')}"

        # runs pca analysis
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, 'merged')} --keep-allele-order --maf 0.01 --out {os.path.join(results_dir, output_name+'.pca')} --pca {pca}"

        # executes Plink command
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        ancestry_fails = self.fail_pca(results_dir, output_name, fails_dir, threshold)

        # create cleaned binary files
        plink_cmd3 = f"plink --bfile {os.path.join(input_path, input_name)} --allow-no-sex --remove {ancestry_fails} --make-bed --out {os.path.join(self.results_dir, output_name+'.clean')}"

        shell_do(plink_cmd3, log=True)

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

    def plot_pca(self, label:str)->None:

        dependables = self.dependables
        results_dir = self.results_dir
        output_name = self.output_name

        step = "pca_analysis"

        # load .eigenvec file
        df_eigenvec = pd.read_csv(
            os.path.join(results_dir, output_name+'.pca.eigenval'),
            header=None,
            sep=' '
        )
        eigenvecs_mat = df_eigenvec[df_eigenvec.columns[2:4]].copy()
        eigenvecs_mat['group'] = label

        # runs pca analysis
        plink_cmd = f"plink --bfile {os.path.join(dependables, '1000g')} --keep-allele-order --out {os.path.join(dependables, '1000g.pca')} --pca 3"

        shell_do(plink_cmd, log=True)

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

    @staticmethod
    def fail_pca(folder_path:str, file_name:str, output_folder:str, threshold:int):

        # load .eigenvec file
        df_eigenvec = pd.read_csv(
            os.path.join(folder_path, file_name+'.pca.eigenvec'),
            header=None,
            sep='\s+'
        )
        eigenvecs_mat = df_eigenvec[df_eigenvec.columns[2:]].copy()

        means = eigenvecs_mat.mean()
        std   = eigenvecs_mat.std()

        for k in eigenvecs_mat.columns:
            eigenvecs_mat[k] = (np.abs(eigenvecs_mat[k] -means[k]) > threshold*std[k] )

        df_outs = df_eigenvec[df_eigenvec.columns[:2]].copy()
        df_outs['is_outlier'] = (np.sum(eigenvecs_mat, axis=1) >0)

        df_outs = df_outs[df_outs['is_outlier']].reset_index(drop=True).drop(columns='is_outlier')

        df_outs.to_csv(
            os.path.join(output_folder, file_name+'.fail-ancestry-qc.txt'),
            header=None,
            index=False,
            sep=' '
        )

        return os.path.join(output_folder, file_name+'.fail-ancestry-qc.txt')
    
    

#        pre_plot = 'ya vamos'
#        pre_plot
#
#        self.plot_pca(
#            results_dir, output_name, 'our G', self.dependables
#        )
