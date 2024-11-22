"""
Module to perform sample quality control
"""

import os
import psutil
import warnings

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import Normalize
from matplotlib import colormaps

from cge_comrare_pipeline.Helpers import shell_do, delete_temp_files

class SampleQC:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, dependables_path:str)->None:#, dependables_path:str, use_kingship:str) -> None:

        """
        Initialize the SampleQC object.

        Parameters:
        -----------
        - input_path (str): Path to the input data.
        - input_name (str): Name of the input data files (without extension).
        - output_path (str): Path to the folder where the output will be saved.
        - output_name (str): Name of the output data.
        - config_dict (str): Dictionary containing configuration settings.
        - dependables_path (str): Path to dependent files.

        Raises:
        ------
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
        

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.dependables = dependables_path
        #self.config_dict = config_dict
        #self.use_kingship= use_kingship

        self.files_to_keep = ['fail_samples']

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

    def execute_ld_pruning(self, ind_pair:list)->dict:
        
        """
        Executes LD pruning using PLINK.

        This method performs linkage disequilibrium (LD) pruning on genetic data using PLINK.
        It excludes high LD regions, performs an LD prune indep-pairwise test, and extracts
        pruned SNPs.

        This method runs three PLINK commands:
        1. Exclude complex regions according to LD.
        2. Run independent pairwise tests and generates prune.in and prune.out files.
        3. Extract pruned SNPs and generate a new dataset.

        Parameters:
        -----------
        ind_pair (list): A list containing three integers representing the window size,
                        step size, and r^2 threshold for the indep-pairwise test.
        
        Returns:
        --------
        dict: A dictionary containing the following keys:
                - 'pass' (bool): Indicates if the process completed successfully.
                - 'step' (str): The name of the step executed.
                - 'output' (dict): A dictionary with the key 'plink_out' pointing to the results directory.
        
        Raises:
        -------
            FileNotFoundError: If the file with high LD regions is not found.
            TypeError: If ind_pair is not a list.
            TypeError: If the first two elements in ind_pair values are not integers.
            TypeError: If the third element in ind_pair is not a float.
        """

        input_path      = self.input_path
        input_name      = self.input_name
        dependables_path= self.dependables
        results_dir     = self.results_dir
        
        # check existence of high LD regions file
        high_ld_regions_file = os.path.join(dependables_path, 'high-LD-regions.txt')
        if not os.path.exists(high_ld_regions_file):
            raise FileNotFoundError(f"File with high LD region was not found: {high_ld_regions_file}")
        
        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        if not isinstance(ind_pair[0], int) or not isinstance(ind_pair[1], int):
            raise TypeError("The first two elements in ind_pair values should be integers (windows size and step size)")
        
        if not isinstance(ind_pair[2], float):
            raise TypeError("The third element in ind_pair should be a float (r^2 threshold)")

        step = "ld_prune"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # exclude complex regions
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --exclude {high_ld_regions_file} --make-bed --out {os.path.join(results_dir, input_name+'-LDregionExcluded')}"

        # LD prune indep-pairwise test
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, input_name+'-LDregionExcluded')} --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --make-bed --out {os.path.join(results_dir, input_name+'-LDregionExcluded-prunning')}"

        plink_cmd3 = f"plink --bfile {os.path.join(results_dir, input_name+'-LDregionExcluded')} --extract {os.path.join(results_dir, input_name+'-LDregionExcluded-prunning.prune.in')} --make-bed --out {os.path.join(results_dir, input_name+'.LDpruned')} --memory {memory} --threads {max_threads}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2, plink_cmd3]
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

    def execute_miss_genotype(self, mind:float)->dict:
        
        """
        Executes PLINK commands to assess missing genotypes and filter samples based on a given missingness threshold.

        This method runs two PLINK commands:
        1. Computes missingness statistics genome-wide.
        2. Filters samples based on the specified missingness threshold (mind) and generates a new dataset.
        
        Parameters:
        -----------
        mind (float): The missingness threshold for filtering samples.
        
        Returns:
        --------
            dict: A dictionary containing the following keys:
                - 'pass' (bool): Indicates if the process completed successfully.
                - 'step' (str): The name of the step executed.
                - 'output' (dict): A dictionary with the key 'plink_out' pointing to the results directory.

        Raises:
        -------
            TypeError: If mind is not a float.
            ValueError: If mind is not between 0 and 1.
            UserWarning: If mind is outside the recommended range
        """


        input_name      = self.input_name
        output_name     = self.output_name
        results_dir     = self.results_dir

        if not isinstance(mind, float):
            raise TypeError("mind should be a float")
        
        # Check if mind is in range
        if mind < 0 or mind > 1:
            raise ValueError("mind should be between 0 and 1")
        
        # Check if mind is around typical values
        if mind <= 0.02 and mind >= 0.1:
            warnings.warn(f"The 'mind' value {mind} is outside the recommended range of 0.02 to 0.1.", UserWarning)

        step = "outlying_missing_genotype"

        # run mssingness across file genome-wide
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --missing --out {os.path.join(results_dir, output_name+'-missing')}"

        # produce a log file with samples excluded at CR 80% and generate plots
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --mind {mind} --make-bed --out {os.path.join(results_dir, output_name+'.-mind')}"

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
    
    def execute_sex_check(self, sex_check:list)->dict:

        """
        Executes a sex check using PLINK commands on genetic data.

        Parameters:
        -----------
        sex_check (list): A list containing two float elements that represent the sex check thresholds. 
                  The list should either be empty or contain exactly two float elements that sum to 1. 
                  If None or an empty list is provided, the samples will not be excluded and the estimates 
        
        Returns:
        --------
        dict: A dictionary containing the following keys:
              - 'pass': A boolean indicating if the process completed successfully.
              - 'step': A string describing the step performed.
              - 'output': A dictionary with the key 'plink_out' pointing to the results directory.
        
        Raises:
        -------
        TypeError: If sex_check is not a list or if its elements are not floats.
        ValueError: If sex_check does not have exactly two elements or if the elements do not sum to 1.
        """

        input_name = self.input_name
        output_name= self.output_name
        result_path= self.results_dir
        results_dir= self.results_dir

        # check type sex_check
        if sex_check is not None:
            if not isinstance(sex_check, list):
                TypeError('sex_check should be a list or None')
            if len(sex_check)>2:
                ValueError('sex_check should have a maximum of two elements')
            if len(sex_check)==1:
                ValueError('sex_check should have two elements or an empty list')
        
            if len(sex_check)>0:
                if not isinstance(sex_check[0], float) or not isinstance(sex_check[1], float):
                    TypeError('elements of sex_check should be floats')
                if sum(sex_check)!=1:
                    ValueError('elements of sex_check should sum to 1')
        
        step = "discordant sex information"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # run sex checking
        if sex_check is None or sex_check==[]:
            plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --check-sex --out {os.path.join(result_path, output_name+'-sexcheck')}"
        else:
            plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --check-sex {sex_check[0]} {sex_check[1]} --threads {max_threads} --out {os.path.join(result_path, output_name+'-sexcheck')}"

        # extract xchr SNPs
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --chr 23 --make-bed --out {os.path.join(result_path, output_name+'-xchr')}"

        # run missingness on xchr SNPs
        plink_cmd3 = f"plink --bfile {os.path.join(result_path, output_name+'-xchr')} --missing --out {os.path.join(result_path, output_name+'-xchr-missing')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2, plink_cmd3]
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

    def execute_heterozygosity_rate(self, maf:float)->dict:
        
        """
        Executes the heterozygosity rate calculation using PLINK commands.
        
        This method performs a series of PLINK commands to calculate the heterozygosity rate
        for a given minor allele frequency (MAF) threshold. It extracts autosomal SNPs, filters
        SNPs based on the MAF threshold, computes missingness, and converts files for heterozygosity computation.

        Finally, computes heterozygosity for each individual and writes the results to a summary file.
        
        Parameters:
        -----------
            maf (float): The minor allele frequency threshold. Must be a float between 0 and 0.5.
        
        Returns:
        --------
            dict: A dictionary containing the following keys:
                - 'pass' (bool): Indicates if the process completed successfully.
                - 'step' (str): The name of the step executed.
                - 'output' (dict): A dictionary with the key 'plink_out' pointing to the results directory.

        Raises:
        -------
            TypeError: If the provided MAF is not a float.
            ValueError: If the provided MAF is not between 0 and 0.5.
        """

        input_name = self.input_name
        output_name= self.output_name
        results_dir= self.results_dir

        if not isinstance(maf, float):
            raise TypeError("maf should be a float")
        if maf < 0 or maf >0.5:
            raise ValueError("maf should be between 0 and 0.5")

        step = "heterozygosity_rate"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # extract autosomal SNPS
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --autosome --make-bed --out {os.path.join(results_dir, output_name+'-chr1-22')}"

        # extract SNPs with minor allele frequency greater than threshold
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, output_name+'-chr1-22')} --maf {maf} --make-bed --out {os.path.join(results_dir, output_name+'-chr1-22-mafgreater')}"

        # extract SNPs with minor allele frequency less than threshold
        plink_cmd3 = f"plink --bfile {os.path.join(results_dir, output_name+'-chr1-22')} --exclude {os.path.join(results_dir, output_name+'-chr1-22-mafgreater.bim')} --make-bed --out {os.path.join(results_dir, output_name+'-chr1-22-mafless')}"

        # get missingness to plot against het
        plink_cmd4 = f"plink --bfile {os.path.join(results_dir, output_name+'-chr1-22-mafgreater')} --missing --out {os.path.join(results_dir, output_name+'-chr1-22-mafgreater-missing')}"
        plink_cmd5 = f"plink --bfile {os.path.join(results_dir, output_name+'-chr1-22-mafless')} --missing --out {os.path.join(results_dir, output_name+'-chr1-22-mafless-missing')}"

        # convert both to ped/map files for heterozigosity computation
        plink_cmd6 = f"plink --bfile {os.path.join(results_dir, output_name+'-chr1-22-mafgreater')} --recode --out {os.path.join(results_dir, output_name+'-chr1-22-mafgreater-recode')} --memory {memory} --threads {max_threads}"
        plink_cmd7 = f"plink --bfile {os.path.join(results_dir, output_name+'-chr1-22-mafless')} --recode --out {os.path.join(results_dir, output_name+'-chr1-22-mafless-recode')} --memory {memory} --threads {max_threads}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2, plink_cmd3, plink_cmd4, plink_cmd5, plink_cmd6, plink_cmd7]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self.compute_heterozigozity(
            ped_file=os.path.join(results_dir, output_name+'-chr1-22-mafgreater-recode.ped')
        )
        self.compute_heterozigozity(
            ped_file=os.path.join(results_dir, output_name+'-chr1-22-mafless-recode.ped')
        )

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
    
    def execute_dup_relatedness(self, kingship:float)->dict:
        
        """
        Executes the duplicates and relatedness analysis using PLINK2.

        This function computes the kinship-coefficient matrix for all samples and prunes for monozygotic twins or duplicates
        based on the provided kinship threshold. It uses PLINK2 commands to perform these operations and logs the execution.

        Parameters:
        -----------
            kingship (float): The kinship threshold value. Must be a float between 0 and 1.
        
        Returns:
        --------
            dict: A dictionary containing the following keys:
                - 'pass' (bool): Indicates if the process completed successfully.
                - 'step' (str): The name of the step executed.
                - 'output' (dict): A dictionary with the key 'plink_out' pointing to the results directory.

        Raises:
        -------
            TypeError: If the kingship parameter is not a float.
            ValueError: If the kingship parameter is not between 0 and 1.
        """

        input_name = self.input_name
        output_name= self.output_name
        results_dir= self.results_dir

        if not isinstance(kingship, float):
            raise TypeError("kingship should be a float")
        if kingship < 0 or kingship >1:
            raise ValueError("kingship should be between 0 and 1")
        
        step = "duplicates_and_relatedness"
        
        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)
        
        # Compute kinship-coefficient matrix for all samples
        plink2_cmd1 = f"plink2 --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --make-king triangle bin --out {os.path.join(results_dir, output_name+'-kinship-coefficient-matrix')} --memory {memory} --threads {max_threads}"

        # Prune for Monozygotic Twins OR Duplicates
        plink2_cmd2 = f"plink2 --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --king-cutoff {os.path.join(results_dir, output_name+'-kinship-coefficient-matrix')} {kingship} --out {os.path.join(results_dir, output_name+'-kinship-pruned-duplicates')} --memory {memory} --threads {max_threads}"

        # execute PLINK commands
        cmds = [plink2_cmd1, plink2_cmd2]
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

    
    @staticmethod
    def compute_heterozigozity(ped_file: str, map_file: str=None)->None:
        
        """
        Computes heterozygosity statistics for individuals in a PED file and writes the results to a summary file.

        The function reads genotype data from the specified PED file, calculates the total number of genotypes,
        the number of homozygous genotypes, the number of heterozygous genotypes, and their respective percentages
        for each individual. The results are written to a summary file in the same directory as the PED file.
        The summary file contains the following columns:
            - ID: Individual ID
            - total: Total number of genotypes
            - num_hom: Number of homozygous genotypes
            - num_het: Number of heterozygous genotypes
            - Percent_hom: Percentage of homozygous genotypes
            - Percent_het: Percentage of heterozygous genotypes

        Parameters:
        -----------
            ped_file (str): Path to the PED file containing genotype data.
            map_file (str, optional): Path to the MAP file. This parameter is currently not used. Defaults to None.

        Returns:
        --------
            None

        
        If the PED file is not found or an I/O error occurs, an error message is printed.
        """
        
        # Define output file name
        summary_file = f"Summary-{os.path.basename(ped_file)}"
        output_path = os.path.join(os.path.dirname(ped_file), summary_file)

        try:
            with open(ped_file, 'r') as ped, open(output_path, 'w') as output:
                # Write the header to the summary file
                output.write("ID\ttotal\tnum_hom\tnum_het\tPercent_hom\tPercent_het\n")

                for line in ped:
                    line = line.strip()
                    if not line:
                        continue
                    
                    # Split the line into columns
                    columns = line.split()
                    individual_id = columns[1]  # Individual ID (second column)
                    genotype_data = columns[6:]  # Genotype data starts at the 7th column

                    # Initialize counters
                    total = 0
                    hom = 0
                    het = 0

                    # Iterate through genotype pairs
                    for i in range(0, len(genotype_data), 2):
                        allele1 = genotype_data[i]
                        allele2 = genotype_data[i + 1]

                        if allele1 == allele2:
                            if allele1 not in ['0', 'N']:  # Exclude missing alleles
                                hom += 1
                                total += 1
                        elif allele1 != allele2:
                            het += 1
                            total += 1

                    # Calculate percentages
                    hom_percent = (hom / total) * 100 if total > 0 else 0.0
                    het_percent = (het / total) * 100 if total > 0 else 0.0

                    # Write the statistics to the output file
                    output.write(f"{individual_id}\t{total}\t{hom}\t{het}\t"
                                 f"{hom_percent:.2f}\t{het_percent:.2f}\n")

            print(f"Summary written to {summary_file}")
        except FileNotFoundError:
            print(f"Error: File {ped_file} not found.")
        except IOError as e:
            print(f"Error: {e}")

    def get_fail_samples(self)->None:

        input_path  = self.input_path
        input_name  = self.input_name
        result_path = self.results_dir
        output_name = self.output_name

        step = "get_fail samples"

        # load samples who failed call rate check
        imiss_file = os.path.join(result_path, output_name+'-missing.imiss')
        df_call_rate = pd.read_csv(
            imiss_file,
            sep='\s+',
            engine='python'
        )

        return df_call_rate

    @staticmethod
    def plot_imiss_het(logFMISS:pd.Series, meanHET:pd.Series, figs_folder:str)->None:

        """
        Plot missing genotypes proportion vs. heterozygosity rate.

        This static method plots the proportion of missing genotypes vs. heterozygosity rate for individuals. It visualizes the relationship between missing data and heterozygosity rate.

        Parameters:
        -----------
        - logFMISS (pd.Series): Pandas Series containing the logarithm of the proportion of missing genotypes for individuals.
        - meanHET (pd.Series): Pandas Series containing the mean heterozygosity rate for individuals.
        - figs_folder (str): Path to the folder where the plot will be saved.

        Returns:
        --------
        - None
        """

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
        plt.close()

        return None

    @staticmethod
    def fail_imiss_het(folder_path:str, file_name:str, output_file:str)->tuple:

        """
        Identify samples failing imiss-het quality control.

        This static method reads the .het and .imiss files generated by PLINK from the specified folder path and file name. It computes the mean heterozygosity rate and the logarithm of the proportion of missing genotypes for each sample. Then, it filters samples based on specific criteria (proportion of missing genotypes or mean heterozygosity rate). It saves the samples that failed imiss-het quality control to the specified output file.

        Parameters:
        -----------
        - folder_path (str): Path to the folder containing the .het and .imiss files.
        - file_name (str): Name of the PLINK file (without extension).
        - output_file (str): Path to save the samples that failed imiss-het quality control.

        Returns:
        --------
        - Tuple[pd.Series, pd.Series]: A tuple containing the logarithm of the proportion of missing genotypes (logF_MISS) and the mean heterozygosity rate (meanHet) for all samples.
        """

        # load .het and .imiss files
        df_het = pd.read_csv(
            os.path.join(folder_path, file_name+'.het'),
            sep="\s+",
            engine='python'
        )
        df_imiss = pd.read_csv(
            os.path.join(folder_path, file_name+'.imiss'),
            sep="\s+",
            engine='python'
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
            sep         =' ', 
            index       =False, 
            header      =False
        )

        return df_imiss['logF_MISS'], df_het['meanHet']

    @staticmethod
    def find_fail_ibd(output_prefix:str, results_folder:str, fails_folder:str, ibd_threshold:float)->None:

        # empty dataframe
        to_remove = pd.DataFrame(columns=['FID', 'IID'])

        # load .imiss file
        df_imiss = pd.read_csv(
            os.path.join(results_folder, output_prefix+'_3.imiss'),
            sep='\s+',
            engine='python'
        )
        # load .genome file
        df_genome = pd.read_csv(
            os.path.join(results_folder, output_prefix+'.genome'),
            sep='\s+',
            engine='python'
        )

        # isolate duplicates or related samples
        df_dup = df_genome[df_genome['PI_HAT']>ibd_threshold].reset_index(drop=True)

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
                to_remove.loc[k] = df_1.iloc[k,0:2].to_list()
            elif df_1.iloc[k,2]<df_2.iloc[k,2]:
                to_remove.loc[k] = df_2.iloc[k,0:2].to_list()
            else:
                to_remove.loc[k] = df_1.iloc[k,0:2].to_list()

        to_remove = to_remove.drop_duplicates(keep='last')
        to_remove.to_csv(
            os.path.join(fails_folder, output_prefix+'.fail-IBD1-qc.txt'),
            index=False,
            header=False,
            sep=" "
        )

        return None
    