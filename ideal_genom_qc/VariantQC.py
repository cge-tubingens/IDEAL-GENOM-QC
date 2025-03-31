"""
Python module to perform variant quality control
"""

import os
import psutil
import logging

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from pathlib import Path

from ideal_genom_qc.Helpers import shell_do, delete_temp_files

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class VariantQC:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name: str) -> None:


        if not isinstance(input_path, Path) or not isinstance(output_path, Path):
            raise TypeError("input_path and output_path should be of type Path")
        if not isinstance(input_name, str) or not isinstance(output_name, str):
            raise TypeError("input_name and output_name should be of type str")
        
        if not input_path.exists() or not output_path.exists():
            raise FileNotFoundError("input_path or output_path is not a valid path")
        if not (input_path / f"{input_name}.bed").exists():
            raise FileNotFoundError(".bed file not found")
        if not (input_path / f"{input_name}.fam").exists():
            raise FileNotFoundError(".fam file not found")
        if not (input_path / f"{input_name}.bim").exists():
            raise FileNotFoundError(".bim file not found")
        
        self.input_path = input_path
        self.output_path= output_path
        self.input_name = input_name
        self.output_name= output_name

        self.hwe_results = None

        # create results folder
        self.results_dir = self.output_path / 'variant_qc_results'
        self.results_dir.mkdir(parents=True, exist_ok=True)

        # create fails folder
        self.fails_dir = self.results_dir / 'fail_samples'
        self.fails_dir.mkdir(parents=True, exist_ok=True)

        # create clean files folder
        self.clean_dir = self.results_dir / 'clean_files'
        self.clean_dir.mkdir(parents=True, exist_ok=True)
        
        # create figures folder
        self.plots_dir = self.results_dir / 'plots'
        self.plots_dir.mkdir(parents=True, exist_ok=True)

    def execute_missing_data_rate(self, chr_y: int = 24) -> None:

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

        # check type for chr_y
        if not isinstance(chr_y, int):
            raise TypeError("chr_y should be of type integer.")
        
        if chr_y < 0 or chr_y > 26:
            raise ValueError("chr_y should be between 1 and 26")

        logger.info("Identifying markers with excessive missing rate...")

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # generates  .lmiss and .imiss files for male subjects
        plink_cmd1 = f"plink --bfile {self.input_path / self.input_name} --missing --filter-males --chr {chr_y} --out {self.results_dir / (self.output_name+'-missing-males-only')} --memory {memory}"

        # generates .lmiss and. imiss files for female subjects
        plink_cmd2 = f"plink --bfile {self.input_path / self.input_name} --missing --not-chr {chr_y} --out {self.results_dir / (self.output_name+'-missing-not-y')} --memory {memory}"

        self.males_missing_data = self.results_dir / (self.output_name+'-missing-males-only.lmiss')
        self.females_missing_data = self.results_dir / (self.output_name+'-missing-not-y.lmiss')

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        return

    def execute_different_genotype_call_rate(self) -> None:

        """
        Identify markers with different genotype call rates between cases and controls.

        This function performs a test for different genotype call rates between cases and controls using PLINK.
    
        Returns:
        --------
        dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        logger.info("Identifying markers with different genotype call rates between cases and controls...")

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # generates .missing file
        plink_cmd = f"plink --bfile {self.input_path / self.input_name} --test-missing --out {self.results_dir / (self.output_name+'-case-control-missing')} --memory {memory}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        self.case_control_missing = self.results_dir / (self.output_name+'-case-control-missing.missing')

        return
    
    def execute_hwe_test(self) -> None:

        logger.info('Computing Hardy-Weinberg Equilibrium test...')

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # PLINK command to compute HWE test
        plink_cmd = f"plink --bfile {self.input_path / self.input_name} --hardy --out {self.results_dir / (self.output_name+'-hwe')} --memory {memory}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)
        self.hwe_results = self.output_name+'-hwe.hwe'

        return
    
    def get_fail_variants(self, marker_call_rate_thres: float = 0.2, case_controls_thres: float = 1e-5, hwe_threshold: float = 5e-8) -> pd.DataFrame:
        
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

        # ==========================================================================================================
        #                                             MARKERS WITH MISSING DATA 
        # ==========================================================================================================

        fail_missing_data = self.report_missing_data(
            directory      =self.results_dir, 
            filename_male  =self.males_missing_data, 
            filename_female=self.females_missing_data,
            threshold      =marker_call_rate_thres, 
        )

        # ==========================================================================================================
        #                                             MARKERS WITH DIFFERENT GENOTYPE CALL RATE
        # ==========================================================================================================

        fail_genotype = self.report_different_genotype_call_rate(
            directory=self.results_dir, 
            filename =self.case_control_missing, 
            threshold=case_controls_thres, 
        )

        # ==========================================================================================================
        #                                             MARKERS FAILING HWE TEST
        # ==========================================================================================================

        fail_hwe = self.report_hwe(
            directory=self.results_dir,
            filename=self.hwe_results,
            hwe_threshold=hwe_threshold
        )

        fails = pd.concat([fail_missing_data, fail_genotype, fail_hwe], axis=0, ignore_index=True)

        summary = fails['Failure'].value_counts().reset_index()
        num_dup = fails.duplicated(subset=['SNP']).sum()

        totals = summary.select_dtypes(include="number").sum() - num_dup
        dups_row = pd.DataFrame({'Failure':['Duplicated SNPs'], 'count':[-num_dup]})
        total_row = pd.DataFrame({col: [totals[col] if col in totals.index else "Total"] for col in summary.columns})

        fails = fails.drop_duplicates(subset='SNP', keep='first', inplace=False)

        fails = fails.drop(columns=['Failure'], inplace=False)

        fails.to_csv(self.fails_dir / 'fail_markers.txt', sep='\t', header=False, index=False)

        return pd.concat([summary, dups_row, total_row], ignore_index=True)

    def execute_drop_variants(self, maf: float = 5e-8, geno: float = 0.1, hwe: float = 5e-8) -> None:

        logger.info("Removing markers failing quality control...")

        # create cleaned binary files
        plink_cmd = f"plink --bfile {self.input_path / self.input_name} --exclude {self.fails_dir / 'fail_markers.txt'} --autosome --maf {maf} --hwe {hwe} --geno {geno} --make-bed --out {self.clean_dir / (self.output_name+'-variantQCed')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return

    def report_missing_data(self, directory: str, filename_male: str, filename_female: str, threshold: float, y_axis_cap: int = 10) -> pd.DataFrame:
   
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
        fail_females['Failure'] = 'Missing data rate on females'

        self._make_histogram(df_males['F_MISS'], 'missing_data_male', threshold, 'Ratio of missing data', 'Missing data for males')
        self._make_histogram(df_females['F_MISS'], 'missing_data_female', threshold, 'Ratio of missing data', 'Missing data for females')

        # concatenate female and male subjects who failed QC
        fails = pd.concat([fail_females, fail_males], axis=0)

        return fails

    def report_different_genotype_call_rate(self, directory: str, filename: str, threshold: float) -> pd.DataFrame:
        """
        Reports markers with different genotype call rates based on a given threshold.
        This function reads a .missing file, filters markers with a different genotype call rate
        below the specified threshold, and returns a DataFrame containing these markers.

        Parameters:
        -----------
            directory (str): The directory where the .missing file is located.
            filename (str): The name of the .missing file.
            threshold (float): The threshold for filtering markers based on the P-value.
            plots_dir (str): The directory where plots will be saved (not used in this function).
        
        Returns:
        --------
            pd.DataFrame: A DataFrame containing markers with different genotype call rates
                          below the specified threshold. The DataFrame has two columns:
                          'SNP' and 'Failure', where 'Failure' is set to 'Different genotype call rate'.
        """

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
    
    def report_hwe(self, directory: Path, filename: str, hwe_threshold: float = 5e-8) -> pd.DataFrame:

        df_hwe = pd.read_csv(
            directory / filename,
            sep=r"\s+",
            engine='python'
        )

        fail_hwe = df_hwe[df_hwe['P']<hwe_threshold].reset_index(drop=True)
        fail_hwe = fail_hwe[['SNP']].copy()
        fail_hwe['Failure'] = 'HWE'

        df_all = df_hwe[df_hwe['TEST']=='ALL'].reset_index(drop=True)
        df_all['P'] = df_all['P'].replace(0, np.finfo(float).tiny)

        self._make_histogram(
            values=-np.log10(df_all['P']), 
            output_name='hwe-histogram', 
            threshold=-np.log10(hwe_threshold), 
            x_label='-log10(P) of HWE test', 
            title='HWE test',
            y_lim_cap=10000
        )

        return fail_hwe
    
    def _make_histogram(self, values:pd.Series, output_name:str, threshold: float, x_label: str, title: str, y_lim_cap: float=None)->None:

        """
        Generate a histogram plot of missing data fraction.

        This static method generates a histogram plot of the missing data fraction (F_MISS) for Single Nucleotide Polymorphisms (SNPs).

        Parameters:
        -----------
        - F_MISS (pandas.Series): Array-like object containing the fraction of missing data for each SNP.
        - figs_folder (str): Path to the folder where the histogram plot will be saved.
        - output_name (str): Name of the output histogram plot file.

        Returns:
        --------
        None
        """

        plt.clf()

        fig_path = self.plots_dir / f"{output_name}.pdf"

        plt.hist(values, bins=50, color='#1B9E77')
        plt.xlabel(x_label)
        plt.ylabel('Number of SNPs')
        plt.ylim(0, y_lim_cap if y_lim_cap else None)
        plt.title(title)

        # Draw the vertical line indicating the cut off threshold
        plt.axvline(x=threshold, linestyle='--', color='red')

        plt.savefig(fig_path, dpi=400)
        plt.show()
        plt.close()

        return None
