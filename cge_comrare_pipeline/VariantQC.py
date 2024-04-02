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

        # check if paths are set
        if input_path is None or output_path is None or dependables_path is None:
            raise ValueError("values for input_path, output_path and dependables_path must be set upon initialization.")

        # Check path validity
#        bed_path = os.path.join(input_path, input_name + '.bed')
#        fam_path = os.path.join(input_path, input_name + '.fam')
#        bim_path = os.path.join(input_path, input_name + '.bim')
#
#        bed_check = os.path.exists(bed_path)
#        fam_check = os.path.exists(fam_path)
#        bim_check = os.path.exists(bim_path)
#
#        if not os.path.exists(input_path) or not os.path.exists(output_path):
#            raise FileNotFoundError("input_path or output_path is not a valid path")
#        if not os.path.exists(dependables_path):
#            raise FileNotFoundError("dependables_path is not a valid path")
#        if not bed_check:
#            raise FileNotFoundError(".bed file not found")
#        if not fam_check:
#            raise FileNotFoundError(".fam file not found")
#        if not bim_check:
#            raise FileNotFoundError(".bim file not found")

        self.input_path = input_path
        self.output_path= output_path
        self.input_name = input_name
        self.output_name= output_name
        self.dependables = dependables_path
        self.config_dict = config_dict

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
        Function to identify all markers with an excessive missing rate.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure.
            * 'output': Dictionary containing paths to the generated output files.
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path = self.results_dir
        output_name = self.output_name
        fails_dir   = self.fails_dir
        fig_folder = self.plots_dir

        chr = self.config_dict['chr']

        # check type for chr
        if not isinstance(chr, int):
            raise TypeError("chr should be of type integer.")
        
        if chr < 0 or chr > 26:
            raise ValueError("chr should be between 1 and 26")

        step = 'high_rate_missing_data'

        #
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --missing --filter-males --chr {chr} --out {os.path.join(result_path, output_name+'.clean_m_only')}"

        #
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --missing --not-chr {chr} --out {os.path.join(result_path, output_name+'.clean_not_y')}"

        # execute PLink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        df_males = pd.read_csv(
            os.path.join(result_path, output_name+'.clean_m_only.lmiss'),
            sep="\s+"
        )
        self.make_histogram(df_males['F_MISS'], fig_folder, 'missing_data_male.pdf')

        df_males = df_males[df_males['F_MISS']>0.2].reset_index(drop=True)
        df_males = df_males[['SNP']].copy()
        df_males.to_csv(
            os.path.join(fails_dir, output_name+'.clean_m_only-fail-lmiss-qc.txt'),
            sep=' ',
            header=False,
            index=False
        )

        df_females = pd.read_csv(
            os.path.join(result_path, output_name+'.clean_not_y.lmiss'),
            sep="\s+"
        )
        self.make_histogram(df_females['F_MISS'], fig_folder, 'missing_data_female.pdf')

        df_females = df_females[df_females['F_MISS']>0.2].reset_index(drop=True)
        df_females = df_females[['SNP']].copy()
        df_females.to_csv(
            os.path.join(fails_dir, output_name+'.clean_not_y-fail-lmiss-qc.txt'),
            sep=' ',
            header=False,
            index=False
        )

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
        Funtion to identify test markers for different genotype call rates between cases and controls.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure.
            * 'output': Dictionary containing paths to the generated output files.
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path = self.results_dir
        output_name = self.output_name
        fails_dir   = self.fails_dir

        step = 'different_genotype_case_control'

        # 
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --test-missing --out {os.path.join(result_path, output_name+'.clean_1')}"

        # execute PLink command
        shell_do(plink_cmd, log=True)

        df_missing = pd.read_csv(
            os.path.join(result_path, output_name+'.clean_1.missing'),
            sep='\s+'
        )

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
        Function to remove markers failing quality control.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure.
            * 'output': Dictionary containing paths to the generated output files.
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path= self.results_dir
        output_name= self.output_name
        output_path= self.output_path
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

    @staticmethod
    def make_histogram(F_MISS, figs_folder, output_name):

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
