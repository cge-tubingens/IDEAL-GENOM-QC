"""
Python module to perform variant quality control
"""

import os

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from cge_comrare_pipeline.Helpers import shell_do

class VariantQC:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_dict:str, dependables_path:str) -> None:

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
        if input_path is None or output_path is None or dependables_path is None:
            raise ValueError("values for input_path, output_path and dependables_path must be set upon initialization.")

        self.input_path = input_path
        self.output_path= output_path
        self.input_name = input_name
        self.output_name= output_name
        self.dependables= dependables_path
        self.config_dict= config_dict

        # create results folder if not existent
        self.results_dir = os.path.join(output_path, 'variant_qc_results')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        # create fails folder if not existent
        self.fails_dir = os.path.join(self.results_dir, 'fail_samples')
        if not os.path.exists(self.fails_dir):
            os.mkdir(self.fails_dir)
        
        # create figures folder if not existent
        self.plots_dir = os.path.join(output_path, 'plots')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)

    def missing_data_rate(self)->dict:

        """
        Identify markers with an excessive missing rate.

        This function performs marker missing data analysis on input data using PLINK. It filters markers based on their missing rate.

        Returns:
        --------
        dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.

        Raises:
        -------
        TypeError: If 'chr' in config_dict is not an integer.
        ValueError: If 'chr' in config_dict is not between 0 and 26 (inclusive).
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path= self.results_dir
        output_name= self.output_name
        fails_dir  = self.fails_dir
        fig_folder = self.plots_dir

        chr = self.config_dict['chr']

        # check type for chr
        if not isinstance(chr, int):
            raise TypeError("chr should be of type integer.")
        
        if chr < 0 or chr > 26:
            raise ValueError("chr should be between 1 and 26")

        step = 'high_rate_missing_data'

        # generates  .lmiss and .imiss files for male subjects
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --missing --filter-males --chr {chr} --out {os.path.join(result_path, output_name+'.clean_m_only')}"

        # generates .lmiss and. imiss files for female subjects
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --missing --not-chr {chr} --out {os.path.join(result_path, output_name+'.clean_not_y')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # load .lmiss file for male subjects
        df_males = pd.read_csv(
            os.path.join(result_path, output_name+'.clean_m_only.lmiss'),
            sep="\s+"
        )
        self.make_histogram(df_males['F_MISS'], fig_folder, 'missing_data_male.pdf')

        # filter male subjects
        df_males = df_males[df_males['F_MISS']>0.2].reset_index(drop=True)
        df_males = df_males[['SNP']].copy()
        df_males.to_csv(
            os.path.join(fails_dir, output_name+'.clean_m_only-fail-lmiss-qc.txt'),
            sep=' ',
            header=False,
            index=False
        )

        # load .lmiss file for female subjects
        df_females = pd.read_csv(
            os.path.join(result_path, output_name+'.clean_not_y.lmiss'),
            sep="\s+"
        )
        self.make_histogram(df_females['F_MISS'], fig_folder, 'missing_data_female.pdf')

        # filter female subjects
        df_females = df_females[df_females['F_MISS']>0.2].reset_index(drop=True)
        df_females = df_females[['SNP']].copy()
        df_females.to_csv(
            os.path.join(fails_dir, output_name+'.clean_not_y-fail-lmiss-qc.txt'),
            sep=' ',
            header=False,
            index=False
        )

        # concatenate female and male subjects who failed QC
        df_fails = pd.concat([df_females, df_males], axis=0)
        df_fails.to_csv(
            os.path.join(fails_dir, output_name+'.clean-fail-lmiss-qc.txt'),
            sep=' ',
            header=False,
            index=False
        )

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

    def different_genotype_call_rate(self)->dict:

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
        fails_dir  = self.fails_dir

        step = 'different_genotype_case_control'

        # generates .missing file
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --test-missing --out {os.path.join(result_path, output_name+'.clean_1')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        # load .missing files
        df_missing = pd.read_csv(
            os.path.join(result_path, output_name+'.clean_1.missing'),
            sep='\s+'
        )

        # filter subjects
        df_missing = df_missing[df_missing['P']< 0.00001].reset_index(drop=True)
        df_missing = df_missing[['SNP']].copy()
        df_missing.to_csv(
            os.path.join(fails_dir, output_name+'.clean-fail-diffmiss-qc.txt'),
            header=False,
            index=False
        )

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

    def remove_markers(self)->dict:

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

        maf = self.config_dict['maf']
        geno= self.config_dict['geno']
        mind= self.config_dict['mind']
        hwe = self.config_dict['hwe']

        step = "remove_markers"

        # load markers with high missing rate
        lmiss_path = os.path.join(fails_dir, output_name+'.clean-fail-lmiss-qc.txt')
        if os.path.getsize(lmiss_path)==0:
            df_lmiss = pd.DataFrame()
        else:
            df_lmiss = pd.read_csv(
                lmiss_path,
                header=None,
                index_col=False,
                sep='\s+'
            )
        
        # load markers with different genotype call rate
        df_diffmiss = pd.read_csv(
            os.path.join(fails_dir, output_name+'.clean-fail-diffmiss-qc.txt'),
            header=None,
            index_col=False
        )

        # marge information from previous steps
        df_markers = pd.concat([df_lmiss, df_diffmiss], axis=0)
        df_markers = df_markers\
            .drop_duplicates(keep='first')\
            .sort_values(by=df_markers.columns[0], inplace=False)

        # save markers that failed variant quality control
        df_markers.to_csv(
            os.path.join(fails_dir, output_name+'.clean-fail-markers-qc.txt'),
            header=False,
            index=False
        )

        # create cleaned binary files
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --exclude {os.path.join(fails_dir, output_name+'.clean-fail-markers-qc.txt')} --maf {maf} --mind {mind} --hwe {hwe} --geno {geno} --make-bed --out {os.path.join(result_path, output_name+'.vrnt_clean')}"

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

        fig_path = os.path.join(figs_folder, f"{output_name}.pdf")

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

        return None
