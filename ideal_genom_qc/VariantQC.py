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
from typing import Optional

from ideal_genom_qc.Helpers import shell_do, delete_temp_files

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class VariantQC:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name: str) -> None:
        """
        Initialize the VariantQC class.
        This class handles quality control for genetic variants data stored in PLINK binary format
        (.bed, .bim, .fam files).
        
        Parameters:
        -----------
        input_path: Path 
            Directory path containing input PLINK files
        input_name: str 
            Base name of input PLINK files (without extension)
        output_path: Path 
            Directory path where output files will be saved
        output_name: str 
            Base name for output files
        
        Raises:
        -------
        TypeError: 
            If input_path/output_path are not Path objects or if input_name/output_name are not strings
        FileNotFoundError: 
            If input_path/output_path don't exist or required PLINK files are not found
        
        Attributes:
        -----------
        input_path: Path 
            Path to input directory
        output_path: Path 
            Path to output directory  
        input_name: str 
            Base name of input files
        output_name: str 
            Base name for output files
        hwe_results: 
            Storage for Hardy-Weinberg equilibrium test results
        results_dir: Path 
            Directory for all QC results
        fails_dir: Path 
            Directory for failed samples
        clean_dir: Path 
            Directory for cleaned files
        plots_dir: Path 
            Directory for QC plots
        """

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
        Executes missing data rate analysis using PLINK for male and female subjects separately.
        This method performs two PLINK operations:
        1. Generates .lmiss and .imiss files for male subjects on chromosome Y
        2. Generates .lmiss and .imiss files for all subjects excluding chromosome Y
        
        Parameters
        ----------
        chr_y : int, default=24
            Chromosome Y number in the dataset. Must be between 0 and 26.
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If chr_y is not an integer
        ValueError
            If chr_y is not between 0 and 26
        
        Notes
        -----
        The method uses 2/3 of available system memory for PLINK operations.
        Output files are generated in the results directory with the following naming pattern:
        - {output_name}-missing-males-only.lmiss/.imiss : For male subjects
        - {output_name}-missing-not-y.lmiss/.imiss : For non-Y chromosome data
        The results are stored in self.males_missing_data and self.females_missing_data as Path objects.
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
        Execute test for different genotype call rates between cases and controls using PLINK.

        This method performs the following operations:
        1. Calculates available memory for PLINK execution
        2. Runs PLINK's --test-missing command to identify markers with significantly different
            missing rates between cases and controls
        3. Generates a .missing file with the results

        The method uses approximately 2/3 of available system memory for PLINK execution.

        Returns:
        --------
             None

        Side effects:
        -------------
            - Creates a .missing file in the results directory
            - Sets self.case_control_missing path attribute
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
        """
        Execute Hardy-Weinberg Equilibrium (HWE) test using PLINK.

        This method performs the following steps:
        1. Calculates available memory and allocates 2/3 for the test
        2. Runs PLINK command to compute HWE test on the input binary PLINK files
        3. Saves results to a .hwe output file

        The HWE test is used to assess whether genotype frequencies in a population remain constant 
        across generations under specific conditions.

        Returns:
        --------
            None

        Side effects:
        -------------
            - Creates a .hwe output file in the results directory
            - Sets self.hwe_results to the name of the output file
        """

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
    
    def get_fail_variants(self, marker_call_rate_thres: float = 0.2, case_controls_thres: float = 1e-5, hwe_threshold: float = 5e-8, male_female_y_cap: Optional[int] = None, hwe_y_cap: Optional[int] = None) -> pd.DataFrame:
        """
        Identify and consolidate failing variants based on multiple quality control criteria.
        This method combines the results of three QC checks:
        1. Variants with high missing data rates
        2. Variants with significantly different genotype call rates between cases and controls
        3. Variants failing Hardy-Weinberg equilibrium test
        
        Parameters
        ----------
        marker_call_rate_thres : float, optional
            Threshold for failing variants based on missing data rate (default: 0.2)
        case_controls_thres : float, optional
            P-value threshold for differential missingness between cases and controls (default: 1e-5)
        hwe_threshold : float, optional
            P-value threshold for Hardy-Weinberg equilibrium test (default: 5e-8)
        
        Returns
        -------
        pd.DataFrame
            A summary DataFrame containing:
            - Counts of variants failing each QC criterion
            - Number of variants failing multiple criteria (duplicates)
            - Total number of unique failing variants
        
        Notes
        -----
        - Results are also written to a tab-separated file 'fail_markers.txt'
        - Variants failing multiple criteria are only counted once in the final output file
        """
        

        # ==========================================================================================================
        #                                             MARKERS WITH MISSING DATA 
        # ==========================================================================================================

        fail_missing_data = self.report_missing_data(
            filename_male  =self.males_missing_data, 
            filename_female=self.females_missing_data,
            threshold      =marker_call_rate_thres,
            y_axis_cap     =male_female_y_cap
        )

        # ==========================================================================================================
        #                                             MARKERS WITH DIFFERENT GENOTYPE CALL RATE
        # ==========================================================================================================

        fail_genotype = self.report_different_genotype_call_rate(
            filename =self.case_control_missing, 
            threshold=case_controls_thres, 
        )

        # ==========================================================================================================
        #                                             MARKERS FAILING HWE TEST
        # ==========================================================================================================

        # Make sure self.hwe_results is not None before passing it to report_hwe
        if self.hwe_results is None:
            raise ValueError("HWE results not available. Run execute_hwe_test first.")
            
        fail_hwe = self.report_hwe(
            directory=self.results_dir,
            filename=self.hwe_results,
            hwe_threshold=hwe_threshold,
            y_lim_cap=hwe_y_cap
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
        """
        Execute variant filtering based on quality control parameters using PLINK.

        This method removes variants that fail quality control criteria including minor allele frequency (MAF),
        genotype missingness rate, and Hardy-Weinberg equilibrium (HWE) test.

        Parameters
        ----------
        maf : float, optional
            Minor allele frequency threshold. Variants with MAF below this value are removed.
            Default is 5e-8.
        geno : float, optional
            Maximum per-variant missing genotype rate. Variants with missing rate above this 
            value are removed. Default is 0.1 (10%).
        hwe : float, optional 
            Hardy-Weinberg equilibrium test p-value threshold. Variants with HWE p-value below
            this are removed. Default is 5e-8.

        Returns
        -------
        None
            Creates quality controlled PLINK binary files (.bed, .bim, .fam) in the clean directory
            with suffix '-variantQCed'.
        """

        logger.info("Removing markers failing quality control...")

        # create cleaned binary files
        plink_cmd = f"plink --bfile {self.input_path / self.input_name} --exclude {self.fails_dir / 'fail_markers.txt'} --autosome --maf {maf} --hwe {hwe} --geno {geno} --make-bed --out {self.clean_dir / (self.output_name+'-variantQCed')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return

    def report_missing_data(self, filename_male: Path, filename_female: Path, threshold: float, y_axis_cap: Optional[float] = None) -> pd.DataFrame:
        """
        Analyze and report missing data rates for male and female subjects.
        This method processes missing data information from separate files for male and female subjects,
        creates visualizations of missing data distributions, and identifies SNPs that fail the missing
        data threshold for each sex group.
        
        Parameters
        ----------
        directory : str
            Path to the directory containing the input files
        filename_male : str
            Name of the file containing missing data information for male subjects (.lmiss format)
        filename_female : str
            Name of the file containing missing data information for female subjects (.lmiss format)
        threshold : float
            Maximum allowed missing data rate (between 0 and 1)
        y_axis_cap : int, optional
            Upper limit for y-axis in histogram plots (default is 10)
        
        Returns
        -------
        pd.DataFrame
            A DataFrame containing SNPs that fail the missing data threshold for either sex,
            with columns ['SNP', 'Failure'] where 'Failure' indicates the failing category
        
        Notes
        -----
        The method generates two histogram plots saved as 'missing_data_male' and 
        'missing_data_female' showing the distribution of missing data rates for each sex.
        """
   

        # load .lmiss file for male subjects
        df_males = pd.read_csv(
            filename_male,
            sep=r"\s+",
            engine='python'
        )
        
        ## filter male subjects
        fail_males = df_males[df_males['F_MISS']>=threshold].reset_index(drop=True)
        fail_males = fail_males[['SNP']].copy()
        fail_males['Failure'] = 'Missing data rate on males'

        # load .lmiss file for female subjects
        df_females = pd.read_csv(
            filename_female,
            sep=r"\s+",
            engine='python'
        )
        
        ## filter female subjects
        fail_females = df_females[df_females['F_MISS']>=threshold].reset_index(drop=True)
        fail_females = fail_females[['SNP']].copy()
        fail_females['Failure'] = 'Missing data rate on females'

        self._make_histogram(df_males['F_MISS'], 'missing_data_male', threshold, 'Ratio of missing data', 'Missing data for males', y_lim_cap=y_axis_cap)
        self._make_histogram(df_females['F_MISS'], 'missing_data_female', threshold, 'Ratio of missing data', 'Missing data for females', y_lim_cap=y_axis_cap)

        # concatenate female and male subjects who failed QC
        fails = pd.concat([fail_females, fail_males], axis=0)

        return fails

    def report_different_genotype_call_rate(self, filename: Path, threshold: float) -> pd.DataFrame:
        """
        Reports markers with different genotype call rates based on a given threshold.
        This function reads a .missing file, filters markers with a different genotype call rate
        below the specified threshold, and returns a DataFrame containing these markers.

        Parameters:
        -----------
            directory (str): The directory where the .missing file is located.
            filename (str): The name of the .missing file.
            threshold (float): The threshold for filtering markers based on the P-value.
        
        Returns:
        --------
            pd.DataFrame: A DataFrame containing markers with different genotype call rates
                          below the specified threshold. The DataFrame has two columns:
                          'SNP' and 'Failure', where 'Failure' is set to 'Different genotype call rate'.
        """

        # load .missing file
        df_diffmiss = pd.read_csv(
            filename,
            sep=r"\s+",
            engine='python'
        )

        # filter markers with different genotype call rate
        fail_diffmiss = df_diffmiss[df_diffmiss['P']<threshold].reset_index(drop=True)
        fail_diffmiss = fail_diffmiss[['SNP']].copy()
        fail_diffmiss['Failure'] = 'Different genotype call rate'

        return fail_diffmiss
    
    def report_hwe(self, directory: Path, filename: str, hwe_threshold: float = 5e-8, y_lim_cap: Optional[float] = None) -> pd.DataFrame:
        """
        Generate Hardy-Weinberg Equilibrium (HWE) test report and visualization.

        This method reads HWE test results from a file, identifies variants that fail HWE,
        creates a histogram of -log10(P) values, and returns failed variants.

        Parameters
        ----------
        directory : Path
            Directory path where the HWE test results file is located
        filename : str
            Name of the file containing HWE test results
        hwe_threshold : float, optional
            P-value threshold for HWE test failure (default: 5e-8)

        Returns
        -------
        pd.DataFrame
            DataFrame containing variants that failed HWE test with columns:
            - SNP: variant identifier
            - Failure: reason for failure (always 'HWE')

        Notes
        -----
        The method creates a histogram plot saved as 'hwe-histogram' showing the
        distribution of -log10(P) values from HWE tests.
        """

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
            values=-np.log10(df_all['P']), # type: ignore
            output_name='hwe-histogram', 
            threshold=-np.log10(hwe_threshold), 
            x_label='-log10(P) of HWE test', 
            title='HWE test',
            y_lim_cap=y_lim_cap
        )

        return fail_hwe
    
    def _make_histogram(self, values: pd.Series, output_name: str, threshold: float, x_label: str, title: str, y_lim_cap: Optional[float] = None) -> None:
        """
        Creates a histogram plot with a vertical threshold line and saves it to a PDF file.
        
        Parameters
        ----------
        values : pd.Series
            The data values to plot in the histogram.
        values: pd.Series
            Series containing the numeric values to plot.
        output_name : str
            Name of the output file (without extension).
        threshold : float
            Value where to draw the vertical threshold line.
        x_label : str
            Label for the x-axis.
        title : str
            Title of the plot.
        y_lim_cap : float, optional
            Upper limit for y-axis. If None, automatically determined. Defaults to None.
        
        Returns
        -------
        None
            This function saves the plot to a file and displays it but does not return any value.
        
        Notes
        -----
        The plot is saved as a PDF file in the plots_dir directory defined in the class instance.
        The histogram uses 50 bins and a predefined color (#1B9E77).
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
        plt.show(block=False)
        plt.close()

        return None
