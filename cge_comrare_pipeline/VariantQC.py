"""
Python module to perform variant quality control
"""

import os
import psutil

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from cge_comrare_pipeline.Helpers import shell_do, delete_temp_files

class VariantQC:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str) -> None:

        """
        Initialize the VariantQC object.

        Parameters:
        -----------
        - input_path (str): Path to input data.
        - input_name (str): Name of the input files.
        - output_path (str): Path to store output data.
        - output_name (str): Name of the output files.
        - config_dict (str): Configuration dictionary.
        - dependables_path (str): Path to dependent files.

        Raises:
        -------
        - ValueError: If values for input_path, output_path, and dependables_path are not provided upon initialization.

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
        if input_path is None or output_path is None:
            raise ValueError("values for input_path and output_path must be set upon initialization.")

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
        
        self.input_path = input_path
        self.output_path= output_path
        self.input_name = input_name
        self.output_name= output_name

        self.files_to_keep = ['fail_samples']

        # create results folder if not existent
        self.results_dir = os.path.join(output_path, 'variant_qc_results')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # create fails folder if not existent
        self.fails_dir = os.path.join(self.results_dir, 'fail_samples')
        if not os.path.exists(self.fails_dir):
            os.mkdir(self.fails_dir)

        # create clean files folder if not existent
        self.clean_dir = os.path.join(self.results_dir, 'clean_files')
        if not os.path.exists(self.clean_dir):
            os.mkdir(self.clean_dir)
        
        # create figures folder if not existent
        self.plots_dir = os.path.join(self.results_dir, 'variantQC_plots')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)

    def execute_missing_data_rate(self, chr_y:int=24)->dict:

        """
        Identify markers with an excessive missing rate.

        This function performs marker missing data analysis on input data using PLINK. It filters markers based on their missing rate.

        Returns:
        --------
        dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.

        Raises:
        -------
        TypeError: If 'chr_y' in config_dict is not an integer.
        ValueError: If 'chr_y' in config_dict is not between 0 and 26 (inclusive).
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path= self.results_dir
        output_name= self.output_name

        # check type for chr_y
        if not isinstance(chr_y, int):
            raise TypeError("chr_y should be of type integer.")
        
        if chr_y < 0 or chr_y > 26:
            raise ValueError("chr_y should be between 1 and 26")

        step = 'high_rate_missing_data'

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # generates  .lmiss and .imiss files for male subjects
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --missing --filter-males --chr {chr_y} --out {os.path.join(result_path, output_name+'-missing-males-only')} --memory {memory}"

        # generates .lmiss and. imiss files for female subjects
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --missing --not-chr {chr_y} --out {os.path.join(result_path, output_name+'-missing-not-y')} --memory {memory}"

        self.males_missing_data = os.path.join(result_path, output_name+'-missing-males-only.lmiss')
        self.females_missing_data = os.path.join(result_path, output_name+'-missing-not-y.lmiss')

        # execute PLINK commands
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

    def execute_different_genotype_call_rate(self)->dict:

        """
        Identify markers with different genotype call rates between cases and controls.

        This function performs a test for different genotype call rates between cases and controls using PLINK.
    
        Returns:
        --------
        dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path= self.results_dir
        output_name= self.output_name

        step = 'different_genotype_case_control'

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # generates .missing file
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name)} --test-missing --out {os.path.join(result_path, output_name+'-case-control-missing')} --memory {memory}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        self.case_control_missing = os.path.join(result_path, output_name+'-case-control-missing.missing')

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
    
    def get_fail_variants(self, marker_call_rate_thres:float=0.2, case_controls_thres:float=1e-5)->pd.DataFrame:
        
        """
        Identifies and reports variants that fail quality control checks based on missing data and genotype call rate.

        Parameters:
        ----------
        marker_call_rate_thres (float): Threshold for marker call rate to identify markers with missing data. Default is 0.2.
        case_controls_thres (float): Threshold for genotype call rate to identify markers with different genotype call rates between cases and controls. Default is 1e-5.
        
        Returns:
        --------
        pd.DataFrame: A DataFrame summarizing the counts of different failure types, including duplicated SNPs and total counts.
        """

        result_path= self.results_dir
        fails_dir  = self.fails_dir
        plots_dir  = self.plots_dir

        # ==========================================================================================================
        #                                             MARKERS WITH MISSING DATA 
        # ==========================================================================================================

        fail_missing_data = self.report_missing_data(
            directory      =result_path, 
            filename_male  =self.males_missing_data, 
            filename_female=self.females_missing_data, \
            threshold      =marker_call_rate_thres, 
            plots_dir      =plots_dir
        )

        # ==========================================================================================================
        #                                             MARKERS WITH DIFFERENT GENOTYPE CALL RATE
        # ==========================================================================================================

        fail_genotype = self.report_different_genotype_call_rate(
            directory=result_path, 
            filename =self.case_control_missing, 
            threshold=case_controls_thres, 
            plots_dir=plots_dir
        )

        fails = pd.concat([fail_missing_data, fail_genotype], axis=0, ignore_index=True)

        summary = fails['Failure'].value_counts().reset_index()
        num_dup = fails.duplicated(subset=['SNP']).sum()

        totals = summary.select_dtypes(include="number").sum() - num_dup
        dups_row = pd.DataFrame({'Failure':['Duplicated SNPs'], 'count':[-num_dup]})
        total_row = pd.DataFrame({col: [totals[col] if col in totals.index else "Total"] for col in summary.columns})

        fails = fails.drop_duplicates(subset='SNP', keep='first', inplace=False)

        fails = fails.drop(columns=['Failure'], inplace=False)

        fails.to_csv(os.path.join(fails_dir, 'fail_markers.txt'), sep='\t', header=False, index=False)

        return pd.concat([summary, dups_row, total_row], ignore_index=True)

    def execute_drop_variants(self, maf:float=5e-8, geno:float=0.1, hwe:float=5e-8)->dict:

        """
        Remove markers failing quality control.

        This function removes markers failing quality control based on specified thresholds for minor allele frequency (MAF), genotype call rate (geno), missing genotype rate (mind), and Hardy-Weinberg equilibrium (hwe).
    
        Returns:
        --------
        dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path= self.results_dir
        output_name= self.output_name
        fails_dir  = self.fails_dir
        clean_dir = self.clean_dir

        step = "remove_markers"

        # create cleaned binary files
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name)} --exclude {os.path.join(fails_dir, 'fail_markers.txt')} --autosome --maf {maf} --hwe {hwe} --geno {geno} --make-bed --out {os.path.join(clean_dir, output_name+'-variantQCed')}"

        # execute PLINK command
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

    def report_missing_data(self, directory:str, filename_male:str, filename_female:str, threshold:float, plots_dir:str, y_axis_cap:int=10):
   
        """
        Reports SNPs with missing data rates above a specified threshold for male and female subjects.
        This function reads .lmiss files for male and female subjects, filters SNPs with missing data rates
        above the given threshold, and generates histograms for the missing data rates. It then concatenates
        the filtered SNPs for both male and female subjects and returns them.

        Parameters:
        -----------
        directory (str): The directory where the .lmiss files are located.
        filename_male (str): The filename of the .lmiss file for male subjects.
        filename_female (str): The filename of the .lmiss file for female subjects.
        threshold (float): The threshold for the missing data rate. SNPs with missing data rates above this threshold will be reported.
        plots_dir (str): The directory where the histograms will be saved.
        y_axis_cap (int, optional): The maximum value for the y-axis in the histograms. Default is 10.

        Returns:
        --------
        pd.DataFrame: A DataFrame containing the SNPs that failed the missing data rate threshold for both male and female subjects.
        """

        # load .lmiss file for male subjects
        df_males = pd.read_csv(
            os.path.join(directory, filename_male),
            sep=r"\s+",
            engine='python'
        )
        
        ## filter male subjects
        fail_males = df_males[df_males['F_MISS']>=threshold].reset_index(drop=True)
        fail_males = fail_males[['SNP']].copy()
        fail_males['Failure'] = 'Missing data rate on males'

        # load .lmiss file for female subjects
        df_females = pd.read_csv(
            os.path.join(directory, filename_female),
            sep=r"\s+",
            engine='python'
        )
        
        ## filter female subjects
        fail_females = df_females[df_females['F_MISS']>=threshold].reset_index(drop=True)
        fail_females = fail_females[['SNP']].copy()
        fail_females['Failure'] = 'Missin data rate on females'

        #self.make_histogram(df_males['F_MISS'], fig_folder, 'missing_data_male.pdf')
        #self.make_histogram(df_females['F_MISS'], fig_folder, 'missing_data_female.pdf')

        # concatenate female and male subjects who failed QC
        fails = pd.concat([fail_females, fail_males], axis=0)

        return fails

    def report_different_genotype_call_rate(self, directory:str, filename:str, threshold:float, plots_dir:str):

        # load .missing file
        df_diffmiss = pd.read_csv(
            os.path.join(directory, filename),
            sep=r"\s+",
            engine='python'
        )

        # filter markers with different genotype call rate
        fail_diffmiss = df_diffmiss[df_diffmiss['P']<threshold].reset_index(drop=True)
        fail_diffmiss = fail_diffmiss[['SNP']].copy()
        fail_diffmiss['Failure'] = 'Different genotype call rate'

        return fail_diffmiss
    
    @staticmethod
    def make_histogram(F_MISS:pd.Series, figs_folder:str, output_name:str)->None:

        """
        Generate a histogram plot of missing data fraction.

        This static method generates a histogram plot of the missing data fraction (F_MISS) for Single Nucleotide Polymorphisms (SNPs).

        Parameters:
        -----------
        - F_MISS (array-like): Array-like object containing the fraction of missing data for each SNP.
        - figs_folder (str): Path to the folder where the histogram plot will be saved.
        - output_name (str): Name of the output histogram plot file.

        Returns:
        --------
        None
        """

        values = F_MISS.copy()

        # substitue 0 by machine epsilon
        for k in range(len(F_MISS)):
            if values[k] == 0:
                values[k] = np.finfo(np.float32).eps

        # log10 transform imput data
        Y = np.log10(values)

        fig_path = os.path.join(figs_folder, f"{output_name}.jpeg")

        plt.hist(Y, bins=50, color='red')
        plt.xlabel('Fraction of missing data')
        plt.ylabel('Number of SNPs')
        plt.title('All SNPs')
        plt.xlim(-4, 0)
        plt.ylim(0, 100000)

        # Label y-axis with the 'ylabels' values
        plt.yticks([])
        ylabels = ['0', '20000', '40000', '60000', '80000', '100000']
        plt.gca().set_yticks([int(label) for label in ylabels])
        plt.gca().set_yticklabels(ylabels)

        # Label x-axis with the 'xlabels' values
        plt.xticks([])
        xlabels = ['-4', '-3', '-2', '-1', '0']
        plt.gca().set_xticks([-4, -3, -2, -1, 0])
        plt.gca().set_xticklabels(xlabels)

        # Draw the vertical line indicating the cut off threshold
        plt.axvline(x=np.log10(0.2), linestyle='--', color='black')

        plt.savefig(fig_path)
        plt.close()

        return None
