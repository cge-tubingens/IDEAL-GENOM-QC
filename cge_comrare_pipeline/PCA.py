"""
Module to perform principal component analysis to identify ethnicity outliers.
"""

import os
import subprocess
import shutil

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from cge_comrare_pipeline.Helpers import shell_do, delete_temp_files

class PCA:
    
    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_dict:str, dependables_path:str) -> None:

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

        self.input_path     = input_path
        self.output_path    = output_path
        self.input_name     = input_name
        self.output_name    = output_name
        self.dependables    = dependables_path
        self.config_dict = config_dict

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

        pass

    def shorten_variant_id(self)->dict:

        """
        Function to deal with long variant IDs. It will be done at a later stage.
        """

        input_path       = self.input_path
        input_name       = self.input_name
        dependables_path = self.dependables

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

        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir
        dependables= self.dependables

        reference_panel = 'all_phase3'

        step = "fiter non A-T or G-C snps"

        self.filter_non_AT_or_GC_snps(input_dir=input_path, input_name=input_name, results_dir=results_dir)

        self.filter_non_AT_or_GC_snps(input_dir=dependables, input_name=reference_panel, results_dir=dependables)

        plink_cmd1 = f"plink --bfile  {os.path.join(input_path, input_name)} --chr 1-22 --exclude {os.path.join(results_dir, input_name+'.ac_get_snps')} --make-bed --out {os.path.join(results_dir, input_name+'.no_ac_gt_snps')}"

        plink_cmd2 = f"plink --bfile  {os.path.join(dependables, reference_panel)} --chr 1-22 --exclude {os.path.join(dependables, reference_panel+'.ac_get_snps')} --allow-extra-chr --memory 10240 --make-bed --out {os.path.join(dependables, reference_panel+'.no_ac_gt_snps')}"

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
        results_dir      = self.results_dir

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

        # generates prune.in and prune.out
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.no_ac_gt_snps')} --maf {maf} --geno {geno} --mind {mind} --hwe {hwe} --exclude {high_ld_regions_file} --range --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --out {os.path.join(results_dir, input_name)}"

        # prune and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, input_name+'.no_ac_gt_snps')} --keep-allele-order --extract {os.path.join(results_dir, input_name+'.prune.in')} --make-bed --out {os.path.join(results_dir, input_name+'.pruned')}"

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
    
    def prune_reference_panel(self)->dict:

        input_name = self.input_name
        dependables= self.dependables
        results_dir= self.results_dir

        step = "prune reference panel"

        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.no_ac_gt_snps')} --keep-allele-order --allow-extra-chr --extract {os.path.join(results_dir, input_name+'.prune.in')} --make-bed --out {os.path.join(dependables, 'all_phase3.pruned')}"

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

        input_name       = self.input_name
        dependables = self.dependables
        results_dir = self.results_dir

        step = "chromosome missmatch"

        awk_cmd = f"awk 'BEGIN {{OFS=\"\t\"}} FNR==NR {{a[$2]=$1; next}} ($2 in a && a[$2] != $1) {{print a[$2],$2}}' {os.path.join(results_dir, input_name+'.pruned.bim')} {os.path.join(dependables, 'all_phase3.pruned.bim')} | sed -n '/^[XY]/!p' > {os.path.join(dependables, 'all_phase3.toUpdateChr')}"

        result = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)

        logs = [result.stderr, result.stdout]

        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.pruned')} --allow-extra-chr --update-chr {os.path.join(dependables, 'all_phase3.toUpdateChr')} 1 2 --make-bed --out {os.path.join(dependables, 'all_phase3.updateChr')}"

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

        input_name  = self.input_name
        dependables = self.dependables
        results_dir = self.results_dir

        step = "possition missmatch and allele flips"

        # position missmatch
        awk_cmd1 = f"awk 'BEGIN {{OFS=\"\t\"}} FNR==NR {{a[$2]=$4; next}} ($2 in a && a[$2] != $4)  {{print a[$2],$2}}' {os.path.join(results_dir, input_name+'.pruned.bim')} {os.path.join(dependables, 'all_phase3.pruned.bim')} > {os.path.join(dependables, 'all_phase3.toUpdatePos')}"

        # possible allele flips
        awk_cmd2 = f"awk 'BEGIN {{OFS=\"\t\"}} FNR==NR {{a[$1$2$4]=$5$6; next}} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {{print $2}}' {os.path.join(results_dir, input_name+'.pruned.bim')} {os.path.join(dependables, 'all_phase3.pruned.bim')} > {os.path.join(dependables, 'all_phase3.toFlip')}"

        awks = [awk_cmd1, awk_cmd2]
        logs = []
        for awk in awks:
            result = subprocess.run(awk, shell=True, capture_output=True, text=True)
            logs.append([result.stderr, result.stdout])

        # update positions and flip alleles
        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.updateChr')} --update-map {os.path.join(dependables, 'all_phase3.toUpdatePos')} 1 2 --flip {os.path.join(dependables, 'all_phase3.toFlip')} --make-bed --out {os.path.join(dependables, 'all_phase3.flipped')}"

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

        input_name  = self.input_name
        dependables = self.dependables
        results_dir = self.results_dir

        step = "remove missmatch"

        awk_cmd = f"awk 'BEGIN {{OFS=\"\t\"}} FNR==NR {{a[$1$2$4]=$5$6; next}} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {{print $2}}' {os.path.join(results_dir, input_name+'.pruned.bim')} {os.path.join(dependables, 'all_phase3.flipped.bim')} > {os.path.join(dependables, 'all_phase3.missmatch')}"

        result = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)
        log = [result.stderr, result.stdout]

        plink_cmd = f"plink --bfile {os.path.join(dependables, 'all_phase3.flipped')} --exclude {os.path.join(dependables, 'all_phase3.missmatch')} --make-bed --out {os.path.join(dependables, 'all_phase3.clean')}"

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

        input_name       = self.input_name
        dependables = self.dependables
        results_dir = self.results_dir

        step = "merge reference panel with study data"

        plink_cmd = f"plink --bfile {os.path.join(results_dir, input_name+'.pruned')} --bmerge {os.path.join(dependables, 'all_phase3.clean.bed')} {os.path.join(dependables, 'all_phase3.clean.bim')} {os.path.join(dependables, 'all_phase3.clean.fam')} --make-bed --out {os.path.join(results_dir, input_name+'.merged')}"

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

        step = "pca_analysis"

        # check `pca` type
        if not isinstance(pca, int):
            raise TypeError("pca should be an integer value")

        # runs pca analysis
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.merged')} --keep-allele-order --maf 0.01 --out {os.path.join(results_dir, output_name+'.pca')} --pca {pca}"

        # executes Plink command
        shell_do(plink_cmd1, log=True)

        df = self.population_tags(
            psam_path= os.path.join(dependables, 'all_phase3.psam'),
            study_fam_path=os.path.join(input_path, input_name+'.fam')
        )
        df['ID1'] = df['ID1'].astype(str)

        ancestry_fails = self.pca_fail(
            df_tags      =df, 
            results_dir  =results_dir,
            output_folder=fails_dir,
            output_name  =output_name, 
            threshold    =threshold
        )

        # create cleaned binary files
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --allow-no-sex --remove {ancestry_fails} --make-bed --out {os.path.join(self.results_dir, output_name+'.clean')}"

        self.results_to_keep.append(output_name+'.clean.bed')
        self.results_to_keep.append(output_name+'.clean.bim')
        self.results_to_keep.append(output_name+'.clean.fam')

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

    def pca_plot(self)->dict:

        input_path  = self.input_path
        input_name  = self.input_name
        output_name = self.output_name
        dependables = self.dependables
        results_dir = self.results_dir

        step = "generate pca plots"

        df = self.population_tags(
            psam_path= os.path.join(dependables, 'all_phase3.psam'),
            study_fam_path=os.path.join(input_path, input_name+'.fam')
        )
        df['ID1'] = df['ID1'].astype(str)

        df_eigenvec = pd.read_csv(
            os.path.join(results_dir, output_name+'.pca.eigenvec'),
            header=None,
            sep=' '
        )
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:5]].copy()
        df_eigenvec.columns = ['ID1', 'ID2', 'pc_1', 'pc_2', 'pc_3']
        df_eigenvec['ID1'] = df_eigenvec['ID1'].astype(str)

        df = pd.merge(df_eigenvec, df, on=['ID1', 'ID2'])

        fig1 = sns.scatterplot(data=df, x='pc_1', y='pc_2', hue='SuperPop')
        plt.savefig(os.path.join(self.plots_dir, 'pca.pdf'), format='pdf')

        fig2 = plt.figure()
        ax = fig2.add_subplot(111, projection='3d')

        for s in df['SuperPop'].unique():
            ax.scatter(
                xs=df.pc_1[df.SuperPop==s],
                ys=df.pc_2[df.SuperPop==s],
                zs=df.pc_3[df.SuperPop==s], 
                label=s
            )
        ax.legend()
        plt.savefig(os.path.join(self.plots_dir, 'pca_3d.pdf'), format='pdf')

        delete_temp_files(self.results_to_keep, results_dir)

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

        bim_target = os.path.join(input_dir, input_name+'.bim')
        ac_gt_snps = os.path.join(results_dir, input_name+'.ac_get_snps')

        awk_cmd = f"awk 'BEGIN {{OFS=\"\t\"}} ($5$6 == \"GC\" || $5$6 == \"CG\" || $5$6 == \"AT\" || $5$6 == \"TA\") {{print $2}}' {bim_target} > {ac_gt_snps}"

        result = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)

        logs = [result.stderr, result.stdout]

        return logs
    
    @staticmethod
    def population_tags(psam_path:str, study_fam_path:str)->pd.DataFrame:

        df_psam = pd.read_csv(
            psam_path,
            sep='\t',
            usecols=['#IID', 'SuperPop']
        )

        df_psam['ID'] = 0
        df_psam = df_psam[['ID', '#IID', 'SuperPop']].copy()
        df_psam.columns = ['ID1', 'ID2', 'SuperPop']

        df_fam = pd.read_csv(
            study_fam_path,
            sep=' ',
            header=None,
            index_col=False
        )
        df_fam = df_fam[df_fam.columns[:2]].copy()
        df_fam['SuperPop'] = 'StPop'
        df_fam.columns = ['ID1', 'ID2', 'SuperPop']

        return pd.concat([df_fam, df_psam], axis=0)

    @staticmethod
    def pca_fail(df_tags:pd.DataFrame, results_dir:str, output_folder:str, output_name:str, threshold:int)->str:

        mask1 = (df_tags['SuperPop']=='SAS')
        mask2 = (df_tags['SuperPop']=='StPop')

        df_ref = df_tags[mask1].reset_index(drop=True)
        df_stu = df_tags[mask2].reset_index(drop=True)

        df_eigenvec = pd.read_csv(
            os.path.join(results_dir, output_name+'.pca.eigenvec'),
            header=None,
            sep=' '
        )

        new_col_names = []
        for k in range(df_eigenvec.shape[1]):
            if k<2:
                new_col_names.append(f"ID{k+1}")
            else:
                new_col_names.append(f"pc_{k-1}")
        df_eigenvec.columns = new_col_names

        df_ref = df_ref.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)
        df_stu = df_stu.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)

        mean_ref = df_ref[df_ref.columns[2:]].mean()
        std_ref = df_ref[df_ref.columns[2:]].std()

        outliers = pd.DataFrame(columns=df_ref.columns)
        outliers[df_stu.columns[:2]] = df_stu[df_stu.columns[:2]]

        for col in outliers.columns[2:]:
            outliers[col] = (np.abs(df_stu[col] - mean_ref[col]) > threshold*std_ref[col])

        outliers['is_out'] = (np.sum(outliers.iloc[:,2:], axis=1) >0)

        df = outliers[outliers['is_out']].reset_index(drop=True)[['ID1', 'ID2']].copy()

        df.to_csv(
            os.path.join(output_folder, output_name+'.fail-ancestry-qc.txt'),
            header=None,
            index=False,
            sep=' '
        )

        return os.path.join(output_folder, output_name+'.fail-ancestry-qc.txt')
