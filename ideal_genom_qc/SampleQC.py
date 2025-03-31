"""
Module to perform sample quality control
"""

import os
import psutil
import warnings

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import colormaps
import seaborn as sns

from ideal_genom_qc.Helpers import shell_do
from ideal_genom_qc.get_references import FetcherLDRegions

from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class SampleQC:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name: str, high_ld_file: Path, built: str = '38') -> None:
        
        """
        Initialize SampleQC class for quality control of genetic data.
        This class handles quality control procedures for genetic data files in PLINK format
        (bed, bim, fam). It sets up the directory structure and validates input files.
        
        Parameters
        ----------
        input_path : Path
            Directory path containing the input PLINK files
        input_name : str
            Base name of the input PLINK files (without extension)
        output_path : Path
            Directory path where output files will be saved
        output_name : str
            Base name for output files (without extension)
        high_ld_file : Path
            Path to file containing high LD regions. If not found, will be fetched from package
        built : str, optional
            Genome build version, either '37' or '38' (default='38')
        
        Raises
        ------
        TypeError
            If input types are incorrect
        ValueError
            If genome build version is not '37' or '38'
        FileNotFoundError
            If input paths or required PLINK files are not found
        
        Attributes
        ----------
        renamed_snps : bool
            Flag indicating if SNPs should be renamed
        hh_to_missing : bool
            Flag indicating if heterozygous haploid genotypes should be set to missing
        pruned_file : None
            Placeholder for pruned file path
        results_dir : Path
            Directory for all QC results
        fails_dir : Path
            Directory for failed samples
        clean_dir : Path
            Directory for clean files
        plots_dir : Path
            Directory for QC plots
        """

        if not isinstance(input_path, Path) or not isinstance(output_path, Path):
            raise TypeError("input_path and output_path should be of type Path")
        if not isinstance(input_name, str) or not isinstance(output_name, str):
            raise TypeError("input_name and output_name should be of type str")
        if not isinstance(high_ld_file, Path):
            raise TypeError("high_ld_file should be of type Path")
        
        if not isinstance(built, str):
            raise TypeError("built should be of type str")
        if built not in ['37', '38']:
            raise ValueError("built should be either 37 or 38")
        
        if not input_path.exists() or not output_path.exists():
            raise FileNotFoundError("input_path or output_path is not a valid path")
        if not (input_path / f"{input_name}.bed").exists():
            raise FileNotFoundError(".bed file not found")
        if not (input_path / f"{input_name}.fam").exists():
            raise FileNotFoundError(".fam file not found")
        if not (input_path / f"{input_name}.bim").exists():
            raise FileNotFoundError(".bim file not found")
        if not high_ld_file.is_file():
            logger.info(f"High LD file not found at {high_ld_file}")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions()
            ld_fetcher.get_ld_regions()

            high_ld_file = ld_fetcher.ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_file}")
        
        self.input_path  = Path(input_path)
        self.output_path = Path(output_path)
        self.input_name  = input_name
        self.output_name = output_name
        self.high_ld_file = high_ld_file

        self.renamed_snps = False
        self.hh_to_missing= False
        self.pruned_file = None

        # create results folder
        self.results_dir = self.output_path / 'sample_qc_results'
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

    def execute_rename_snpid(self, rename: bool = True) -> None:

        if not isinstance(rename, bool):
            raise TypeError("rename should be a boolean")
        
        if not rename:
            logger.info(f"STEP: Rename SNPs. `rename` set to {rename}. Skipping renaming of SNPs in the study data")
            return
        else:
            logger.info(f"STEP: Rename SNPs. `rename` set to {rename}. Renaming SNPs in the study data to the format chr_pos_a1_a2")
            self.renamed_snps = True

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        plink2_cmd = f"plink2 --bfile {self.input_path / self.input_name} --set-all-var-ids @:#:$r:$a --threads {max_threads} --make-bed --out {self.input_path / (self.input_name+ '-renamed')}"

        # Execute PLINK2 command
        shell_do(plink2_cmd, log=True)

        return
    
    def execute_haploid_to_missing(self, hh_to_missing: bool = True) -> None:

        if not isinstance(hh_to_missing, bool):
            raise TypeError("renamed should be a boolean")
        
        if not hh_to_missing:
            logger.info(f"STEP: Convert haploid genotypes to missing values. `hh_to_missing` set to {hh_to_missing}. Skipping conversion of haploid genotypes to missing values")
            return
        else:
            logger.info(f"STEP: Convert haploid genotypes to missing values. `hh_to_missing` set to {hh_to_missing}. Converting haploid genotypes to missing values in the study data")
            self.hh_to_missing = True
        
        logger.info("STEP: Convert haploid genotypes to missing values")

        if self.renamed_snps:
            # PLINK command: convert haploid genotypes to missing
            plink_cmd = f"plink --bfile {self.input_path / (self.input_name+'-renamed')} --set-hh-missing --keep-allele-order --make-bed --out {self.input_path / (self.input_name+'-hh-missing')}"
        else:
            # PLINK command: convert haploid genotypes to missing
            plink_cmd = f"plink --bfile {self.input_path / self.input_name} --set-hh-missing --keep-allele-order --make-bed --out {self.input_path / (self.input_name+'-hh-missing')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return
    
    def execute_ld_pruning(self, ind_pair: list = [50, 5, 0.2]) -> None:
        
        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        if not isinstance(ind_pair[0], int) or not isinstance(ind_pair[1], int):
            raise TypeError("The first two elements in ind_pair values should be integers (windows size and step size)")
        
        if not isinstance(ind_pair[2], float):
            raise TypeError("The third element in ind_pair should be a float (r^2 threshold)")

        logger.info("STEP: LD pruning")

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        if self.hh_to_missing:
            ld_input = self.input_name+'-hh-missing'
        elif self.renamed_snps:
            ld_input = self.input_name+'-renamed'
        else:
            ld_input = self.input_name

        # exclude complex regions
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name+'-hh-missing')} --exclude {high_ld_regions_file} --make-bed --out {os.path.join(results_dir, input_name+'-LDregionExcluded')}"

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
            'pass'  : process_complete,
            'step'  : step,
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

        input_name  = self.input_name
        output_name = self.output_name
        results_dir = self.results_dir

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
        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --mind {mind} --make-bed --out {os.path.join(results_dir, output_name+'-mind')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self.call_rate_miss = os.path.join(results_dir, output_name+'-missing.imiss')

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass'  : process_complete,
            'step'  : step,
            'output': outfiles_dict
        }

        return out_dict
    
    def execute_sex_check(self, sex_check:list=[])->dict:

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

        self.sexcheck_miss = os.path.join(result_path, output_name+'-sexcheck.sexcheck')
        self.xchr_miss = os.path.join(result_path, output_name+'-xchr-missing.imiss')

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass'  : process_complete,
            'step'  : step,
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

        self.summary_greater = os.path.join(results_dir, 'Summary-'+output_name+'-chr1-22-mafgreater-recode.ped')
        self.summary_less    = os.path.join(results_dir, 'Summary-'+output_name+'-chr1-22-mafless-recode.ped')
        self.maf_greater_miss= os.path.join(results_dir, output_name+'-chr1-22-mafgreater-missing.imiss')
        self.maf_less_miss   = os.path.join(results_dir, output_name+'-chr1-22-mafless-missing.imiss')

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass'  : process_complete,
            'step'  : step,
            'output': outfiles_dict
        }

        return out_dict

    def execute_ibd(self)->dict:

        input_name = self.input_name
        output_name= self.output_name
        results_dir= self.results_dir

        step = "identity_by_descent"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # PLINK commands
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --genome --out {os.path.join(results_dir, output_name+'-ibd')} --threads {max_threads}"

        plink_cmd2 = f"plink --bfile {os.path.join(results_dir, input_name+'.LDpruned')} --allow-no-sex --missing --out {os.path.join(results_dir, output_name+'-ibd-missing')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self.ibd_miss = os.path.join(results_dir, output_name+'-ibd-missing.imiss')
        self.genome = os.path.join(results_dir, output_name+'-ibd.genome')

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass'  : process_complete,
            'step'  : step,
            'output': outfiles_dict
        }

        return out_dict

    def execute_kingship(self, kingship:float)->dict:
        
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

        self.kinship_miss = os.path.join(results_dir, output_name+'-kinship-pruned-duplicates.king.cutoff.out.id')

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass'  : process_complete,
            'step'  : step,
            'output': outfiles_dict
        }

        return out_dict
    
    def execute_duplicate_relatedness(self, kingship:float=0.354, use_king:bool=True)->dict:

        step = "duplicates_and_relatedness"
        
        if use_king:
            self.execute_kingship(kingship)
        else:
            self.execute_ibd()

        self.use_king = use_king

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': self.results_dir
        }

        out_dict = {
            'pass'  : process_complete,
            'step'  : step,
            'output': outfiles_dict
        }

        return out_dict

    def execute_recover_snp_names(self, rename:bool=True)->dict:
        
        """
        Executes the recovery of SNP names by renaming them based on a linkage file and 
        updating the corresponding .bim file. Optionally, it can skip the renaming process.

        Parameters:
        -----------
            rename (bool): Flag to determine whether to rename SNPs or not. Default is True.

        Returns:
        --------
            dict: A dictionary containing the status of the process, the step name, and the output file directory.

        Raises:
        -------
            FileNotFoundError: If the file linking old names and new names is not found.
        """

        input_path = self.input_path
        input_name = self.input_name
        output_name= self.output_name

        step = "recover_snp_names"

        if not rename:
            pass

        if not os.path.isfile(os.path.join(input_path, input_name + '.linkage')):
            raise FileNotFoundError("Linkage file not found.")

        df_bim_cleaned = pd.read_csv(
            os.path.join(self.clean_dir, output_name+'-clean-samples.bim'),
            sep   =r'\s+',
            engine='python',
            header=None,
        )

        df_linkage = pd.read_csv(
            os.path.join(input_path, input_name + '.linkage'),
            sep   =r'\s+',
            engine='python',
            header=None,
        )

        df_bim_cleaned[1] = df_linkage[0].copy()
        df_bim_cleaned.to_csv(
            os.path.join(self.clean_dir, output_name+'-clean-samples.bim'),
            sep   ='\t',
            header=False,
            index =False,
        )

        # PLINK command
        plink_cmd = f"plink --bfile {os.path.join(self.clean_dir, output_name+'-clean-samples')} --make-bed --out {os.path.join(self.clean_dir, output_name+'-clean-samples')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': self.clean_dir
        }

        out_dict = {
            'pass'  : process_complete,
            'step'  : step,
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
        summary_file= f"Summary-{os.path.basename(ped_file)}"
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
                    total= 0
                    hom  = 0
                    het  = 0

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

    def get_fail_samples(self, call_rate_thres:float, std_deviation_het:float, maf_het:float, ibd_threshold:float)->pd.DataFrame:
        
        """
        Identifies and reports samples that fail various quality control checks.

        Parameters:
        -----------
        call_rate_thres : float
            The threshold for the call rate check. Samples with call rates above this threshold will be flagged.
        std_deviation_het : float
            The standard deviation threshold for the heterozygosity rate check. Samples with heterozygosity rates threshold*std away from the mean will be flagged.
        maf_het : float
            The minor allele frequency threshold for the heterozygosity rate check. Used to split the heterozygosity rate check into two parts: MAF > threshold and MAF < threshold.

        Returns:
        --------
        pandas.DataFrame
            The function saves the results of the failed samples to a file and returns a summary DataFrame of the failure counts.
        """

        result_path = self.results_dir
        output_name = self.output_name
        plots_dir   = self.plots_dir

        # ==========================================================================================================
        #                                             CALL RATE CHECK
        # ==========================================================================================================

        # load samples who failed call rate check
        fail_call_rate = self.report_call_rate(
            directory    =result_path, 
            filename     =output_name+'-missing.imiss', 
            threshold    =call_rate_thres, 
            plots_dir    =plots_dir, 
            y_axis_cap   =10
        )

        print('Call rate check done')

        # ==========================================================================================================
        #                                             SEX CHECK
        # ==========================================================================================================

        fail_sexcheck = self.report_sex_check(
            directory          =result_path, 
            sex_check_filename =output_name+'-sexcheck.sexcheck', 
            xchr_imiss_filename=output_name+'-xchr-missing.imiss',
            plots_dir          =plots_dir
        )

        print('Sex check done')

        # ==========================================================================================================
        #                                       HETETROZYGOSITY RATE CHECK
        # ==========================================================================================================

        fail_het_greater = self.report_heterozygosity_rate(
            directory           = result_path, 
            summary_ped_filename= 'Summary-'+output_name+'-chr1-22-mafgreater-recode.ped', 
            autosomal_filename  = output_name+'-chr1-22-mafgreater-missing.imiss', 
            std_deviation_het   = std_deviation_het,
            maf                 = maf_het,
            split               = '>',
            plots_dir           = plots_dir
        )

        print('Heterozygosity rate check done for MAF > threshold')

        fail_het_less = self.report_heterozygosity_rate(
            directory           = result_path, 
            summary_ped_filename= 'Summary-'+output_name+'-chr1-22-mafless-recode.ped', 
            autosomal_filename  = output_name+'-chr1-22-mafless-missing.imiss', 
            std_deviation_het   = std_deviation_het,
            maf                 = maf_het,
            split               = '<',
            plots_dir           = plots_dir
        )

        print('Heterozygosity rate check done for MAF < threshold')

        # ==========================================================================================================
        #                                       DUPLICATES-RELATEDNESS CHECK
        # ==========================================================================================================

        if self.use_king:

            # load samples that failed duplicates and relatedness check
            duplicates_file = os.path.join(result_path, output_name+'-kinship-pruned-duplicates.king.cutoff.out.id')
            df_duplicates = pd.read_csv(
                duplicates_file,
                sep   =r'\s+',
                engine='python'
            )

            # filter samples that failed duplicates and relatedness check
            df_duplicates.columns = ['FID', 'IID']
            fail_duplicates = df_duplicates[['FID', 'IID']].reset_index(drop=True)
            fail_duplicates['Failure'] = 'Duplicates and relatedness (Kingship)'

            print('Duplicates and relatedness check done with kingship')

        else:

            fail_duplicates = self.report_ibd_analysis(ibd_threshold)

            print('Duplicates and relatedness check done with IBD')

        # ==========================================================================================================
        #                                       MERGE ALL FAILURES
        # ==========================================================================================================

        fails = [fail_call_rate, fail_sexcheck, fail_het_greater, fail_het_less, fail_duplicates] 

        df = pd.concat(fails, axis=0).reset_index(drop=True)

        summary = df['Failure'].value_counts().reset_index()
        num_dup = df.duplicated(subset=['FID', 'IID']).sum()

        df = df.drop_duplicates(subset=['FID', 'IID'])

        df.to_csv(os.path.join(self.fails_dir, 'fail_samples.txt'), index=False, sep='\t')

        totals = summary.select_dtypes(include="number").sum() - num_dup

        # Create the total row
        dups_row = pd.DataFrame({'Failure':['Duplicated Sample IDs'], 'count':[-num_dup]})
        total_row= pd.DataFrame({col: [totals[col] if col in totals.index else "Total"] for col in summary.columns})

        # Append the total row to the DataFrame
        summary = pd.concat([summary, dups_row, total_row], ignore_index=True)
        
        return summary
    
    def execute_drop_samples(self)->dict:
        
        """
        Executes the process of dropping samples using PLINK.

        This method constructs and runs a PLINK command to remove samples listed in a specified file.
        The resulting files are saved to the specified output directory with a modified name.

        Returns:
        --------
        dict: A dictionary containing the following keys:
                - 'pass' (bool): Indicates if the process completed successfully.
                - 'step' (str): The name of the step executed.
                - 'output' (dict): A dictionary with the key 'plink_out' pointing to the results directory.
        """

        input_path = self.input_path
        input_name = self.input_name
        clean_dir  = self.clean_dir
        output_name= self.output_name
        fails_dir  = self.fails_dir

        step = "drop_samples"

        # drop samples
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name+'-hh-missing')} --remove {os.path.join(fails_dir, 'fail_samples.txt')} --make-bed --out {os.path.join(clean_dir, output_name+'-clean-samples')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': clean_dir
        }

        out_dict = {
            'pass'  : process_complete,
            'step'  : step,
            'output': outfiles_dict
        }

        return out_dict
  
    def report_call_rate(self, directory:str, filename:str, threshold:float, plots_dir:str, y_axis_cap:int=10)->pd.DataFrame:
        
        """
        Generates a report on sample call rates, including histograms and scatter plots, and identifies samples that fail the call rate threshold.

        Parameters:
        -----------
        directory (str): 
            The directory where the input file is located.
        filename (str): 
            The name of the input file containing sample call rate data.
        threshold (float): 
            The threshold for the proportion of missing SNPs (F_MISS) above which samples are considered to have failed the call rate.
        plots_dir (str): 
            The directory where the output plots will be saved.
        y_axis_cap (int, optional): 
            The maximum value for the y-axis in the capped histogram. Default is 10.
        
        Returns:
        --------
        pandas.DataFrame: 
            A DataFrame containing the samples that fail the call rate threshold, with columns 'FID', 'IID', and 'Failure'.
        """

        # load samples that failed sex check
        df_call_rate = pd.read_csv(
            os.path.join(directory, filename),
            sep=r'\s+',
            engine='python'
        )

        # filter samples that fail call rate
        fail_call_rate = df_call_rate[df_call_rate['F_MISS'] > threshold][['FID', 'IID']].reset_index(drop=True)
        fail_call_rate['Failure'] = 'Call rate'

        # Create the figure and subplots
        fig1, axes1 = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

        # First subplot: Full histogram
        axes1[0] = sns.histplot(df_call_rate['F_MISS'], bins=30, color='blue', alpha=0.7, ax=axes1[0])
        axes1[0].set_title("Sample Call Rate Distribution")
        axes1[0].set_xlabel("Proportion of missing SNPs (F_MISS)")
        axes1[0].set_ylabel("Frequency")

        # Second subplot: Histogram with capped y-axis
        axes1[1] = sns.histplot(df_call_rate['F_MISS'], bins=30, color='blue', alpha=0.7, ax=axes1[1])
        axes1[1].set_ylim(0, y_axis_cap)  # Cap y-axis
        axes1[1].set_title("Sample Call Rate Distribution (Capped)")
        axes1[1].set_xlabel("Proportion of missing SNPs (F_MISS)")

        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, f"call_rate_{threshold}_histogram.jpeg"), dpi=400)
        plt.show(block=False)

        fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

        # First subplot: capped y-axis
        axes2[0] = sns.histplot(df_call_rate['F_MISS'], bins=50, color='blue', alpha=0.7, ax=axes2[0])
        axes2[0].set_ylim(0, y_axis_cap)  # Cap y-axis
        axes2[0].set_title("Sample Call Rate Distribution (Capped)")
        axes2[0].set_xlabel("Proportion of missing SNPs (F_MISS)")

        # Add a vertical line at the threshold
        axes2[0].axvline(threshold, linewidth=2, color='firebrick', linestyle='dashed')

        # Second subplot: Number of samples vs F_MISS
        df_call_rate_sorted = pd.DataFrame({
            'Index': range(len(df_call_rate['F_MISS'])),
            'F_MISS': sorted(df_call_rate['F_MISS'])
        })

        axes2[1] = sns.scatterplot(
            data  =df_call_rate_sorted,
            x     ='Index',
            y     ='F_MISS',
            marker='o',  
            color ='blue',
            ax    =axes2[1]
        ) 
        axes2[1].set_title("Sample Call Rate")
        axes2[1].set_xlabel(f"Number of samples")
        axes2[1].set_ylabel("F_MISS")

        # Add a vertical line at the threshold
        axes2[1].axhline(threshold, linewidth=2, color='firebrick', linestyle='dashed')

        # third subplot: Number of samples vs F_MISS
        axes2[2] = sns.scatterplot(
            x      =df_call_rate['F_MISS'],
            y      =np.random.normal(size=len(df_call_rate['F_MISS'])),
            markers='o',
            s      =20,
        )
        axes2[2].set_title("Sample Call Rate")
        axes2[2].set_xlabel("Proportion of missing SNPs (F_MISS)")
        axes2[2].set_ylabel(f"Samples")
    

        # Add a vertical line at the threshold
        axes2[2].axvline(threshold, linewidth=2, color='firebrick', linestyle='dashed')

        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, f"call_rate_{threshold}_scatterplot.jpeg"), dpi=400)
        plt.show(block=False)

        return fail_call_rate
    
    def report_sex_check(self, directory:str, sex_check_filename:str, xchr_imiss_filename:str, plots_dir:str)->pd.DataFrame:
        
        """
        Generates a report for sex check and creates a scatter plot visualizing the results.
        
        Parameters:
        -----------
        directory (str): 
            The directory where the input files are located.
        sex_check_filename (str): 
            The filename of the sex check data.
        xchr_imiss_filename (str): 
            The filename of the X chromosome missingness data.
        plots_dir (str): 
            The directory where the plot will be saved.
        
        Returns:
        --------
        pandas.DataFrame: 
            A DataFrame containing the FID and IID of samples that failed the sex check.
        """

        df_sexcheck = pd.read_csv(
            os.path.join(directory, sex_check_filename),
            sep   =r'\s+',
            engine='python'
        )

        df_xchr_imiss = pd.read_csv(
            os.path.join(directory, xchr_imiss_filename),
            sep   =r'\s+',
            engine='python'
        )

        df = pd.merge(df_sexcheck, df_xchr_imiss, on=['FID', 'IID'], how='inner')

        fail_sexcheck = df[df['STATUS'] == 'PROBLEM'][['FID', 'IID']].reset_index(drop=True)
        fail_sexcheck['Failure'] = 'Sex check'

        df['Category'] = 'General'
        df.loc[df['PEDSEX'] == 1, 'Category'] = 'Male PEDSEX'
        df.loc[df['PEDSEX'] == 2, 'Category'] = 'Female PEDSEX'

        df_problem = df[df['STATUS'] == 'PROBLEM'].reset_index(drop=True)
        df = df[df['STATUS'] != 'PROBLEM'].reset_index(drop=True)

        # Define the palette (color mapping)
        palette = {
            "Male PEDSEX"  : "blue",
            "Female PEDSEX": "green"
        }

        # Define the size mapping
        size_mapping = {
            "Male PEDSEX"  : 40,
            "Female PEDSEX": 40
        }

        # Create the Matplotlib scatter plot
        fig, ax = plt.subplots(figsize=(8, 6))

        # Iterate through categories to plot each group separately
        for category, group in df.groupby("Category"):
            ax.scatter(
                group["F"], 
                group["F_MISS"], 
                edgecolors=palette[category],     # Map color
                facecolors='none',                # Hollow circles
                s         =size_mapping[category],# Map size
                label     =category               # Add label for legend
            )

        ax.scatter(
            df_problem["F"], 
            df_problem["F_MISS"], 
            color     ='red',
            s         =25,
            marker    ='o',
            label     ='Problem Status',
            edgecolors=palette['Female PEDSEX'],
        )

        # Add vertical lines
        plt.axvline(x=0.8, color='red', linestyle='dotted')
        plt.axvline(x=0.2, color='red', linestyle='dotted')

        # Customize labels and legend
        plt.title("Sex Check")
        plt.xlabel("X chr inbreeding (homozygosity) estimate F")
        plt.ylabel("Proportion of missing SNPs for the X chr")
        plt.legend(title='', loc='best')
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'sex_check.jpeg'), dpi=400)

        return fail_sexcheck
    
    def report_heterozygosity_rate(self, directory:str, summary_ped_filename:str, autosomal_filename:str, std_deviation_het:float, maf:float, split:str, plots_dir:str, y_axis_cap:float=80)->pd.DataFrame:
        
        """
        Generates a report on heterozygosity rate and plots histograms and scatter plots for visualization.

        Parameters:
        -----------
        directory (str): 
            The directory where the input files are located.
        summary_ped_filename (str): 
            The filename of the summary PED file containing heterozygosity rates.
        autosomal_filename (str): 
            The filename of the autosomal file containing call rates.
        std_deviation_het (float): 
            The standard deviation threshold for heterozygosity rate exclusion.
        maf (float): 
            Minor Allele Frequency threshold.
        split (str): 
            A string identifier for the split (e.g., 'train', 'test').
        plots_dir (str): 
            The directory where the plots will be saved.
        y_axis_cap (float, optional): 
            The cap for the y-axis in the histogram plot. Default is 80.
        
        Returns:
        --------
        pandas.DataFrame: 
            A DataFrame containing samples that failed the heterozygosity rate check.
        """

        # load samples that failed heterozygosity rate check with MAF > threshold
        maf_file = os.path.join(directory, summary_ped_filename)
        df_maf = pd.read_csv(
            maf_file,
            sep   =r'\s+',
            engine='python'
        )

        # autosomal call rate per individual
        autosomal_file = os.path.join(directory, autosomal_filename)
        df_autosomal = pd.read_csv(
            autosomal_file,
            sep   =r'\s+',
            engine='python'
        )

        # merge both dataframes
        df_het = pd.merge(
            df_maf[['ID', 'Percent_het']],
            df_autosomal[['FID', 'IID', 'F_MISS']],
            left_on ='ID',
            right_on='IID',
            how     ='inner'
        )

        mean_percent= df_het['Percent_het'].mean()
        sd_percent  = df_het['Percent_het'].std()

        mask_plus = df_het['Percent_het'] > mean_percent + std_deviation_het*sd_percent
        mask_minus= df_het['Percent_het'] < mean_percent - std_deviation_het*sd_percent

        fail_het = df_het[mask_plus | mask_minus][['FID', 'IID']].reset_index(drop=True)

        if split == '>':
            fail_het['Failure'] = 'Heterozygosity rate greater'
        else:
            fail_het['Failure'] = 'Heterozygosity rate less'

        # plots

        fig1, axes1 = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

        axes1[0] = sns.histplot(df_het['Percent_het'], bins=30, color='green', alpha=0.7, ax=axes1[0])
        axes1[0].set_title("Autosomal heterozygosity")
        axes1[0].set_xlabel(f"% Heterozygosity MAF {split} {maf}")
        axes1[0].set_ylabel("Frequency")

        axes1[1] = sns.histplot(df_het['Percent_het'], bins=30, color='green', alpha=0.7, ax=axes1[1])
        axes1[1].set_title("Autosomal heterozygosity (capped)")
        axes1[1].set_xlabel(f"% Heterozygosity MAF {split} {maf}")
        axes1[1].set_ylim(0, y_axis_cap)  # Cap y-axis
        axes1[1].set_ylabel("Frequency")

        plt.tight_layout()
        if split == '>':
            plt.savefig(os.path.join(plots_dir, f"heterozygosity_rate_greater_{maf}_histogram.jpeg"), dpi=400)
        else:
            plt.savefig(os.path.join(plots_dir, f"heterozygosity_rate_less_{maf}_histogram.jpeg"), dpi=400)
        
        plt.show(block=False)

        df_het['Deviated'] = 'Not Excluded'
        df_het.loc[mask_plus, 'Deviated'] = f'{std_deviation_het}xSD Excluded'
        df_het.loc[mask_minus, 'Deviated']= f'{std_deviation_het}xSD Excluded'

        # Create the scatter plot
        plt.figure(figsize=(10, 6))
        sns.scatterplot(
            data   =df_het,
            x      ='Percent_het',
            y      ='F_MISS',
            hue    ='Deviated',
            palette={'Not Excluded': 'blue', f'{std_deviation_het}xSD Excluded': 'red'},
            markers={'Not Excluded': 'o', f'{std_deviation_het}xSD Excluded': 'o'},
            size   ='Deviated',
            sizes  ={'Not Excluded': 20, f'{std_deviation_het}xSD Excluded': 30}
        )
        plt.title("Autosomal heterozygosity and call rate")
        plt.xlabel(f"% Heterozygosity MAF {split} {maf}")
        plt.ylabel("Proportion of missing SNPs")
        plt.legend(title='Exclusion', loc='best')

        plt.tight_layout()
        if split == '>':
            plt.savefig(os.path.join(plots_dir, f"heterozygosity_rate_greater_{maf}_scatterplot.jpeg"), dpi=400)
        else:
            plt.savefig(os.path.join(plots_dir, f"heterozygosity_rate_less_{maf}_scatterplot.jpeg"), dpi=400)
        plt.show(block=False)

        return fail_het

    def report_ibd_analysis(self, ibd_threshold: float = 0.185, chunk_size: int = 100000) -> pd.DataFrame:
        """
        Identify individuals with high IBD (identity-by-descent) relatedness and report samples to remove.

        Parameters:
            ibd_threshold (float): Threshold for IBD (PI_HAT) above which samples are considered related.
            chunk_size (int): Number of rows to process at a time from the .genome file.

        Returns:
            pd.DataFrame: A dataframe containing FID, IID, and the reason for removal.
        """
        if not isinstance(ibd_threshold, float):
            raise TypeError("ibd_threshold should be a float")

        if self.use_king:
            return pd.DataFrame()

        # File paths
        results_dir = self.results_dir
        output_name = self.output_name

        imiss_path = os.path.join(results_dir, output_name + '-ibd-missing.imiss')
        genome_path= os.path.join(results_dir, output_name + '-ibd.genome')

        if not os.path.exists(imiss_path):
            raise FileNotFoundError(f"Missing file: {imiss_path}")
        if not os.path.exists(genome_path):
            raise FileNotFoundError(f"Missing file: {genome_path}")

        # Load .imiss file
        df_imiss = pd.read_csv(imiss_path, sep=r'\s+', engine='python')

        # Initialize dataframe for duplicates
        duplicates = []

        # Process the .genome file in chunks
        for chunk in pd.read_csv(
            genome_path,
            usecols  =['FID1', 'IID1', 'FID2', 'IID2', 'PI_HAT'],
            sep      =r'\s+',
            engine   ='python',
            chunksize=chunk_size,
        ):
            # Filter rows with PI_HAT > ibd_threshold
            filtered_chunk = chunk[chunk['PI_HAT'] > ibd_threshold]
            if not filtered_chunk.empty:
                duplicates.append(filtered_chunk)

        if not duplicates:
            return pd.DataFrame(columns=['FID', 'IID', 'Failure'])

        # Concatenate all filtered chunks
        df_dup = pd.concat(duplicates, ignore_index=True)

        # Merge with missingness data
        imiss_related1 = pd.merge(
            df_dup[['FID1', 'IID1']],
            df_imiss[['FID', 'IID', 'F_MISS']],
            left_on =['FID1', 'IID1'],
            right_on=['FID', 'IID'],
        ).rename(columns={'F_MISS': 'F_MISS_1'})

        imiss_related2 = pd.merge(
            df_dup[['FID2', 'IID2']],
            df_imiss[['FID', 'IID', 'F_MISS']],
            left_on =['FID2', 'IID2'],
            right_on=['FID', 'IID'],
        ).rename(columns={'F_MISS': 'F_MISS_2'})

        # Decide which samples to remove
        to_remove = pd.concat(
            [
                imiss_related1[['FID1', 'IID1', 'F_MISS_1']],
                imiss_related2[['FID2', 'IID2', 'F_MISS_2']],
            ],
            axis=1,
        )

        to_remove['FID'], to_remove['IID'] = np.where(
            to_remove['F_MISS_1'] > to_remove['F_MISS_2'],
            (to_remove['FID1'], to_remove['IID1']),
            (to_remove['FID2'], to_remove['IID2']),
        )

        to_remove = to_remove[['FID', 'IID']].drop_duplicates().reset_index(drop=True)
        to_remove['Failure'] = 'Duplicates and relatedness (IBD)'

        return to_remove

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
        norm  = Normalize(vmin=min(logFMISS), vmax=max(logFMISS))
        colors= colormaps['viridis']

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
            sep   =r"\s+",
            engine='python'
        )
        df_imiss = pd.read_csv(
            os.path.join(folder_path, file_name+'.imiss'),
            sep   =r"\s+",
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
            sep         ='\t', 
            index       =False, 
            header      =False
        )

        return df_imiss['logF_MISS'], df_het['meanHet']
