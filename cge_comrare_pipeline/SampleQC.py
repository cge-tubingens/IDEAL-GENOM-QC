"""
Python module to perform sample quality control
"""

import os
import json

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import colormaps

from cge_comrare_pipeline.Helpers import shell_do

class SampleQC:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_path:str, dependables_path:str) -> None:

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
        self.fst_pruned_dir = None

        # open config file
        with open(config_path, 'r') as file:
            self.config_dict = json.load(file)

        # create results folder
        self.results_dir = os.path.join(output_path, 'sample_qc_results')
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

    def run_ld_prune(self, ld_region_file:str)->dict:

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
        high_ld_regions_file = os.path.join(dependables_path, ld_region_file)
        if not os.path.exists(high_ld_regions_file):
            raise FileNotFoundError("File with high LD region was not found")

        step = "ld_prune"

        # create figures folder
        self.fst_pruned_dir = os.path.join(output_path, '1st_pruned')
        if not os.path.exists(self.fst_pruned_dir):
            os.mkdir(self.fst_pruned_dir)
        fst_pruned_dir= self.fst_pruned_dir

        # generates prune.in and prune.out
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --maf {maf} --geno {geno} --mind {mind} --hwe {hwe} --exclude {high_ld_regions_file} --range --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --out {os.path.join(fst_pruned_dir, output_name)}"

        # prune and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --extract {os.path.join(fst_pruned_dir, output_name+'.prune.in')} --make-bed --out {os.path.join(fst_pruned_dir, output_name)}"

        # execute Plink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': fst_pruned_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def run_heterozygosity_rate(self)->dict:

        """
        Function to identify individuals with elevated missing data rates or outlying heterozygosity rate.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        output_name= self.output_name
        result_path= self.results_dir
        plots_path = self.plots_dir
        fst_pruned_dir= self.fst_pruned_dir
        fails_dir= self.fails_dir

        step = "heterozygosity_rate"

        # create .imiss and .lmiss files
        plink_cmd1 = f"plink --bfile {os.path.join(fst_pruned_dir, output_name)} --keep-allele-order --missing --out {os.path.join(result_path, output_name+'_1')}"

        # create .het file
        plink_cmd2 = f"plink --bfile {os.path.join(fst_pruned_dir, output_name)} --keep-allele-order --het --autosome --extract {os.path.join(fst_pruned_dir, output_name+'.prune.in')} --out {os.path.join(result_path, output_name+'_1')}"

        # execute PLink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        fails_path = os.path.join(fails_dir, output_name+'.fail-imisshet-qc.txt')
        logFMISS, meanHET = self.fail_imiss_het(
            result_path, output_name+'_1',
            fails_path
        )

        # generate plot
        self.plot_imiss_het(
            logFMISS, meanHET, plots_path
        )

        # create cleaned binary file
        plink_cmd3 = f"plink --bfile {os.path.join(fst_pruned_dir, output_name)} --keep-allele-order --remove {fails_path} --make-bed --out {os.path.join(result_path, output_name+'_1')}"

        # execute PLink command
        shell_do(plink_cmd3, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': result_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def run_sex_check(self)->dict:

        """
        Function to identify individuals with discordant sex information.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        output_name= self.output_name
        result_path= self.results_dir
        fst_pruned_dir= self.fst_pruned_dir
        fails_dir = self.fails_dir

        sex_check = self.config_dict['sex_check']

        # check type sex_check
        if not isinstance(sex_check, list):
            raise TypeError("sex_check should be a list")
        if len(sex_check)!=2:
            raise ValueError("sex_check must be a list of length 2")
        
        for value in sex_check:
            if not isinstance(value, float):
                raise TypeError("sex_check values should be float")
            if 0 > value or value > 1:
                raise ValueError("sex_check values must be between 0 and 1")
        
        if sum(sex_check) != 1:
            raise ValueError("sex_check values should add to 1")
        
        step = "sex_check"

        # create .sexcheck file
        plink_cmd1 = f"plink --bfile {os.path.join(result_path, output_name+'_1')} --check-sex {sex_check[0]} {sex_check[1]} --keep-allele-order --extract {os.path.join(fst_pruned_dir, output_name+'.prune.in')} --out {os.path.join(result_path, output_name+'_2')}"

        # execute PLink command
        shell_do(plink_cmd1, log=True)

        # load file with sex analysis
        df = pd.read_csv(
            os.path.join(result_path, output_name+'_2.sexcheck'),
            sep='\s+'
        )

        # filter problematic samples and save file
        df_probs = df[df['STATUS']=='PROBLEM'].reset_index(drop=True)
        df_probs.to_csv(
            os.path.join(result_path, output_name+'_2.sexprobs'),
            index =False,
            sep   ='\t'
        )

        # save IDs of samples who failed sex check QC
        df_probs = df_probs.iloc[:,0:2].copy()
        df_probs.to_csv(
            os.path.join(fails_dir, output_name+'.fail-sexcheck-qc.txt'),
            index  =False,
            header =False,
            sep    =' '
        )

        # generate clean file
        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'_1')} --keep-allele-order --remove {os.path.join(fails_dir, output_name+'.fail-sexcheck-qc.txt')} --make-bed --out {os.path.join(result_path, output_name+'_2')}"

        # execute PLink command
        shell_do(plink_cmd2, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': result_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def run_relatedness_prune(self)->dict:

        """
        Function to identify duplicated or related individuals.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        result_path= self.results_dir
        output_name= self.output_name
        fails_dir  = self.fails_dir
        fst_pruned_dir=self.fst_pruned_dir

        step = "duplicate_relative_prune"

        to_remove = pd.DataFrame(columns=['FID', 'IID'])

        # run genome
        plink_cmd1 = f"plink --bfile {os.path.join(result_path, output_name+'_2')} --keep-allele-order --extract {os.path.join(fst_pruned_dir, output_name+'.prune.in')} --genome --out {os.path.join(result_path, output_name+'_3')}"

        # Generate new .imiss file
        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'_2')} --keep-allele-order --missing --out {os.path.join(result_path, output_name+'_3')}"

        # execute PLink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # load .imiss file
        df_imiss = pd.read_csv(
            os.path.join(result_path, output_name+'_3.imiss'),
            sep='\s+'
        )
        # load .genome file
        df_genome = pd.read_csv(
            os.path.join(result_path, output_name+'_3.genome'),
            sep='\s+'
        )

        # isolate duplicates or related samples
        df_dup = df_genome[df_genome['PI_HAT']>0.185].reset_index(drop=True)

        df_1 = pd.merge(
            df_dup[['FID1', 'IID1']], 
            df_imiss[['FID', 'IID', 'F_MISS']], 
            left_on =['FID1', 'IID1'],
            right_on=['FID', 'IID']
        ).drop(columns=['FID', 'IID'], inplace=False)

        df_2 = pd.merge(
            df_dup[['FID2', 'IID2']], 
            df_imiss[['FID', 'IID', 'F_MISS']], 
            left_on=['FID2', 'IID2'],
            right_on=['FID', 'IID']
        ).drop(columns=['FID', 'IID'], inplace=False)

        for k in range(len(df_dup)):

            if df_1.iloc[k,2]>df_2.iloc[k,2]:
                print('sample', df_1.iloc[k,0:2].to_list())
                to_remove.loc[k] = df_1.iloc[k,0:2].to_list()
            elif df_1.iloc[k,2]<df_2.iloc[k,2]:
                to_remove.loc[k] = df_2.iloc[k,0:2].to_list()
            else:
                to_remove.loc[k] = df_1.iloc[k,0:2].to_list()

        to_remove = to_remove.drop_duplicates(keep='last')
        to_remove.to_csv(
            os.path.join(fails_dir, output_name+'.fail-IBD1-qc.txt'),
            index=False,
            header=False,
            sep=" "
        )

        # create cleaned binary files
        plink_cmd3 = f"plink --bfile {os.path.join(result_path, output_name+'_2')} --keep-allele-order --remove {os.path.join(fails_dir, output_name+'.fail-IBD1-qc.txt')} --make-bed --out {os.path.join(result_path, output_name+'_3')}"

        # execute PLink command
        shell_do(plink_cmd3, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': result_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def delete_failing_QC(self)->None:

        """
        Function to remove samples that failed quality control.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path= self.results_dir
        output_name= self.output_name
        fails_dir  = self.fails_dir
        fst_pruned_dir=self.fst_pruned_dir

        step = "delete_sample_failed_QC"


        # load files with samples who failed one or several QC steps
        df_sex = pd.read_csv(
            os.path.join(fails_dir, output_name+'.fail-sexcheck-qc.txt'),
            sep      =' ',
            index_col=False,
            header   =None
        )

        df_imiss = pd.read_csv(
            os.path.join(fails_dir, output_name+'.fail-imisshet-qc.txt'),
            sep      =' ',
            index_col=False,
            header   =None
        )

        df_ibd1 = pd.read_csv(
            os.path.join(fails_dir, output_name+'.fail-IBD1-qc.txt'),
            sep      =' ',
            index_col=False,
            header   =None
        )

        # concatenate all failings samples
        df = pd.concat([df_sex, df_imiss, df_ibd1], axis=0)

        # sort values by the first column
        df.sort_values(by=df.columns[0], inplace=True)

        # drop duplicates
        df = df.drop_duplicates(keep='first')

        # save file with samples who failed QC
        df.to_csv(
            os.path.join(fails_dir, output_name+'.fail-qc_1-inds.txt'),
            sep   =' ',
            header=False,
            index =False
        )

        # generate cleaned binary files
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --remove {os.path.join(fails_dir, output_name+'.fail-qc_1-inds.txt')} --make-bed --out {os.path.join(result_path, output_name+'.pre_ind_clean')}"

        # execute PLink command
        shell_do(plink_cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': result_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def divergent_ancestry_step_one(self, ld_region_file)->dict:

        """
        Function to identify subject with divergent ancestry.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        result_path= self.results_dir
        output_name= self.output_name
        dependables_path = self.dependables

        maf = self.config_dict['maf']
        geno = self.config_dict['geno']
        hwe = self.config_dict['hwe']
        ind_pair = self.config_dict['indep-pairwise']
        mind = self.config_dict['mind']

        step = "preparation_for_pca"

        # generate prune.in and prune.out files
        plink_cmd1 = f"plink --bfile {os.path.join(result_path, output_name+'.pre_ind_clean')} --autosome   --exclude {os.path.join(dependables_path, ld_region_file)} --geno {geno} --hwe {hwe} --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --keep-allele-order --maf {maf} --mind {mind} --out {os.path.join(result_path, output_name+'.pre_ind_clean')} --range"

        # generate clean binary files
        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'.pre_ind_clean')} --extract {os.path.join(result_path, output_name+'.pre_ind_clean.prune.in')} --keep-allele-order --make-bed --out {os.path.join(result_path, output_name+'.pre_ind_clean.pca_ready')}"

        # execute PLink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': result_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def run_pca_analysis(self)->dict:

        """
        Funtion to prunes samples .....

        Parameters:
        - ld_region_file: string
            file name with regions with high Linkage Distribution

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        output_name = self.output_name
        result_path = self.results_dir
        fails_dir   = self.fails_dir
        threshold   = self.config_dict['outlier_threshold']
        pca         = self.config_dict['pca']

        step = "pca_analysis"

        # check `pca` type
        if not isinstance(pca, int):
            raise TypeError("pca should be an integer value")

        # runs pca analysis
        plink_cmd1 = f"plink --bfile {os.path.join(result_path, output_name+'.pre_ind_clean.pca_ready')}   --keep-allele-order --out {os.path.join(result_path, output_name+'.pca')} --pca {pca}"

        # executes Plink command
        shell_do(plink_cmd1, log=True)

        ancestry_fails = self.fail_pca(result_path, output_name, fails_dir, threshold)

        # create folder for cleaned files
        self.clean_samples_dir = os.path.join(self.output_path, 'clean_samples')
        if not os.path.exists(self.clean_samples_dir):
            os.mkdir(self.clean_samples_dir)

        # create cleaned binary files
        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'.pre_ind_clean')} --allow-no-sex --remove {ancestry_fails} --make-bed --out {os.path.join(self.clean_samples_dir, output_name+'.clean')}"

        shell_do(plink_cmd2, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': result_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    @staticmethod
    def plot_imiss_het(logFMISS, meanHET, figs_folder):

        # Calculate colors based on density
        norm = Normalize(vmin=min(logFMISS), vmax=max(logFMISS))
        colors = colormaps['viridis']

        fig_path = os.path.join(figs_folder, "imiss-vs-het.pdf")

        meanHet_low= np.mean(meanHET) - 2*np.std(meanHET)
        meanHet_up = np.mean(meanHET) + 2*np.std(meanHET)

        # Plotting
        plt.figure(figsize=(8, 6))
        plt.scatter(logFMISS, meanHET, cmap=colors, s=50, marker='o', norm=norm)
        plt.xlim(-3, 0)
        plt.ylim(0, 0.5)
        plt.xlabel("Proportion of missing genotypes")
        plt.ylabel("Heterozygosity rate")
        plt.yticks(np.arange(0, 0.51, 0.05))
        plt.xticks([-3, -2, -1, 0], [0.001, 0.01, 0.1, 1])
        plt.axhline(meanHet_low, color='red', linestyle='--')
        plt.axhline(meanHet_up, color='red', linestyle='--')
        plt.axvline(-1.522879, color='red', linestyle='--')
        plt.grid(True)
        plt.savefig(fig_path)

        return None

    @staticmethod
    def fail_imiss_het(folder_path:str, file_name:str, output_file:str):

        # load .het and .imiss files
        df_het = pd.read_csv(
            os.path.join(folder_path, file_name+'.het'),
            sep="\s+"
        )
        df_imiss = pd.read_csv(
            os.path.join(folder_path, file_name+'.imiss'),
            sep="\s+"
        )

        # compute Het mean
        df_het['meanHet'] = (df_het['N(NM)']-df_het['O(HOM)'])/df_het['N(NM)']

        df_imiss['logF_MISS'] = np.log10(df_imiss['F_MISS'])
    
        # compute the lower 2 standard deviation bound
        meanHet_lower = df_het['meanHet'].mean() - 2*df_het['meanHet'].std()

        # compute the upper 2 standard deviation bound
        meanHet_upper = df_het['meanHet'].mean() + 2*df_het['meanHet'].std()

        # filter samples
        mask = ((df_imiss['F_MISS']>=0.04) | (df_het['meanHet'] < meanHet_lower) | (df_het['meanHet'] > meanHet_upper))

        df = df_imiss[mask].reset_index(drop=True)
        df = df.iloc[:,0:2].copy()

        # save samples that failed imiss-het QC
        df.to_csv(
            path_or_buf =output_file, 
            sep         ='\t', 
            index       =False, 
            header      =False
        )

        return df_imiss['logF_MISS'], df_het['meanHet']

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
