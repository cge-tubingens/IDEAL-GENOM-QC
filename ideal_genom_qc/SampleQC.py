"""
Module to perform sample quality control
"""

import os
import psutil
import warnings
import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import colormaps
import seaborn as sns

from ideal_genom_qc.Helpers import shell_do
from ideal_genom_qc.get_references import FetcherLDRegions

from pathlib import Path
from typing import Optional

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
            
            ld_fetcher = FetcherLDRegions(built=built)
            ld_fetcher.get_ld_regions()

            ld_regions = ld_fetcher.ld_regions
            if ld_regions is None:
                raise ValueError("Failed to fetch high LD regions file")
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
        
        """
        Executes the SNP ID renaming process using PLINK2.
        This method renames SNP IDs in the PLINK binary files to a standardized format of 'chr:pos:a1:a2'.
        The renaming is performed using PLINK2's --set-all-var-ids parameter.
        
        Parameter:
        ----------
        rename (bool, optional): Flag to control whether SNP renaming should be performed. 
            Defaults to True.

        Returns:
        --------
            None

        Raises:
        -------
            TypeError: If rename parameter is not a boolean.

        Notes:
        ------
            - The renamed files will be saved with '-renamed' suffix
            - Thread count is optimized based on available CPU cores
            - The new SNP ID format will be: chromosome:position:allele1:allele2
            - Sets self.renamed_snps to True if renaming is performed
        """

        if not isinstance(rename, bool):
            raise TypeError("rename must be a boolean")
        
        if not rename:
            logger.info(f"STEP: Rename SNPs. `rename` set to {rename}. Skipping renaming of SNPs in the study data")
            return
        else:
            logger.info(f"STEP: Rename SNPs. `rename` set to {rename}. Renaming SNPs in the study data to the format chr_pos_a1_a2")
            self.renamed_snps = True

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

        plink2_cmd = f"plink2 --bfile {self.input_path / self.input_name} --set-all-var-ids @:#:$r:$a --threads {max_threads} --make-bed --out {self.input_path / (self.input_name+ '-renamed')}"

        # Execute PLINK2 command
        shell_do(plink2_cmd, log=True)

        return
    
    def execute_haploid_to_missing(self, hh_to_missing: bool = True) -> None:

        """
        Convert haploid genotypes to missing values in PLINK binary files.
        This method uses PLINK's --set-hh-missing flag to convert haploid genotypes to 
        missing values in the genotype data. This is often useful for quality control 
        of genetic data, particularly for variants on sex chromosomes.
        
        Parameters
        ----------
        hh_to_missing : bool, default=True
            If True, converts haploid genotypes to missing values.
            If False, skips the conversion step.
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If hh_to_missing is not a boolean value.
        
        Notes
        -----
        The method uses PLINK to process the binary files (.bed, .bim, .fam) and creates
        new files with suffix '-hh-missing'. The input files are determined based on whether
        SNPs have been previously renamed (checks self.renamed_snps).
        """

        if not isinstance(hh_to_missing, bool):
            raise TypeError("hh_to_missing must be a boolean")
        
        if not hh_to_missing:
            logger.info(f"STEP: Convert haploid genotypes to missing values. `hh_to_missing` set to {hh_to_missing}. Skipping conversion of haploid genotypes to missing values")
            return
        else:
            logger.info(f"STEP: Convert haploid genotypes to missing values. `hh_to_missing` set to {hh_to_missing}. Converting haploid genotypes to missing values in the study data")
            self.hh_to_missing = True
        
        logger.info("STEP: Convert haploid genotypes to missing values")

        # Dynamically set the input file name based on whether SNPs are renamed
        input_file = self.input_name + '-renamed' if self.renamed_snps else self.input_name

        # PLINK command: convert haploid genotypes to missing
        plink_cmd = f"plink --bfile {self.input_path / input_file} --set-hh-missing --keep-allele-order --make-bed --out {self.input_path / (self.input_name+'-hh-missing')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return
    
    def execute_ld_pruning(self, ind_pair: list = [50, 5, 0.2]) -> None:
        """
        Execute LD (Linkage Disequilibrium) pruning on genetic data using PLINK.
        This method performs LD pruning in three steps:
        1. Excludes complex/high LD regions
        2. Identifies SNPs for pruning using indep-pairwise test
        3. Creates final pruned dataset

        Parameters
        ----------
        ind_pair : list, optional
            List of three elements for LD pruning parameters:
            - Window size (int): Number of SNPs to analyze in each window
            - Step size (int): Number of SNPs to shift window at each step
            - r² threshold (float): Correlation coefficient threshold for pruning
            Default is [50, 5, 0.2]
        
        Raises
        ------
        TypeError
            If ind_pair is not a list
            If first two elements of ind_pair are not integers
            If third element of ind_pair is not float
        ValueError
            If ind_pair does not contain exactly three elements
            If window size or step size is not positive
            If r² threshold is not between 0 and 1
        FileNotFoundError
            If required pruning input file is not found
        
        Notes
        -----
        - Uses available CPU cores (leaving 2 cores free) and 2/3 of available memory
        - Creates intermediate and final files with suffixes:
          * '-LDregionExcluded'
          * '-LDregionExcluded-prunning'
          * '-LDpruned'
        - Updates self.pruned_file with path to final pruned dataset
        """
        
        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        if len(ind_pair) != 3:
            raise ValueError("ind_pair must have exactly three elements")
        
        if not isinstance(ind_pair[0], int) or not isinstance(ind_pair[1], int):
            raise TypeError("The first two elements in ind_pair values should be integers (window size and step size)")
        if ind_pair[0] <= 0 or ind_pair[1] <= 0:
            raise ValueError("Window size and step size must be positive integers")
        
        if not isinstance(ind_pair[2], float):
            raise TypeError("The third element in ind_pair should be a float (r^2 threshold)")
        if not (0 < ind_pair[2] <= 1):
            raise ValueError("The r^2 threshold must be a float between 0 and 1")

        logger.info("STEP: LD pruning")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

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
        plink_cmd1 = f"plink --bfile {self.input_path / ld_input} --exclude {self.high_ld_file} --make-bed --out {self.results_dir / (self.input_name+'-LDregionExcluded')}"
        prune_in_file = (self.results_dir / (self.input_name+'-LDregionExcluded-prunning')).with_suffix('.prune.in')


        plink_cmd3 = f"plink --bfile {self.results_dir / (self.input_name+'-LDregionExcluded')} --extract {prune_in_file} --keep-allele-order --make-bed --out {self.results_dir / (self.input_name + '-LDpruned')} --memory {memory} --threads {max_threads}"
        # LD prune indep-pairwise test
        plink_cmd2 = f"plink --bfile {self.results_dir / (self.input_name+'-LDregionExcluded')} --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --keep-allele-order --make-bed --out {self.results_dir / (self.input_name+'-LDregionExcluded-prunning')} --memory {memory} --threads {max_threads}"

        plink_cmd3 = f"plink --bfile {self.results_dir / (self.input_name+'-LDregionExcluded')} --extract {(self.results_dir / (self.input_name+'-LDregionExcluded-prunning')).with_suffix('.prune.in')} --keep-allele-order --make-bed --out {self.results_dir / (self.input_name + '-LDpruned')} --memory {memory} --threads {max_threads}"

        self.pruned_file = self.results_dir / (self.input_name + '-LDpruned')

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2, plink_cmd3]
        for cmd in cmds:
            logger.info(f"Executing PLINK command: {cmd}")
            shell_do(cmd, log=True)

        return
    
    def execute_miss_genotype(self, mind: float = 0.2) -> None:
        """
        Execute missing genotype analysis using PLINK to identify and filter samples with high missingness rates.
        This method performs two main operations:
        1. Generates missingness statistics for all samples
        2. Filters samples based on the specified missingness threshold (mind)
        
        Parameters
        ----------
        mind : float, optional
            The missingness threshold for sample filtering (default is 0.2).
            Samples with missingness rates above this threshold will be removed.
            Recommended range is between 0.02 and 0.1.
        
        Returns
        -------
        dict
            Dictionary containing missingness analysis results
        
        Raises
        ------
        TypeError
            If mind parameter is not a float
        ValueError
            If mind parameter is not between 0 and 1
        FileNotFoundError
            If the output .imiss file is not generated
            If mind value is outside recommended range (0.02-0.1)
        
        Notes
        -----
        This function creates two files:
        - {input_name}-missing.imiss: Contains missingness statistics for all samples
        - {output_name}-mind.bed: New binary file with filtered samples
        """

        if not isinstance(mind, float):
            raise TypeError("mind should be a float")
        
        # Check if mind is in range
        if mind < 0 or mind > 1:
            raise ValueError("mind should be between 0 and 1")
        
        # Check if mind is around typical values
        if mind < 0.02 or mind > 0.1:
            warnings.warn(
                f"The 'mind' value {mind} is outside the recommended range of 0.02 to 0.1. "
                f"Values below 0.02 may exclude too few samples with high missingness, potentially affecting data quality, "
                f"while values above 0.1 may exclude too many samples, reducing the dataset size significantly.",
                UserWarning
            )

        logger.info(f"STEP: Missing genotype check. `mind` set to {mind}")

        # PLINK command: run mssingness across file genome-wide 
        plink_cmd1 = f"plink --bfile {self.pruned_file} --missing --out {self.results_dir / (self.input_name+'-missing')}"

        # PLINK command: produce a log file with samples excluded at CR 80% and generate plots
        plink_cmd2 = f"plink --bfile {self.pruned_file} --mind {mind} --keep-allele-order --make-bed --out {self.results_dir / (self.output_name+'-mind')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self.call_rate_miss = (self.results_dir / (self.input_name+'-missing')).with_suffix('.imiss')
        if not self.call_rate_miss.exists():
            raise FileNotFoundError(f"Missing file: {self.call_rate_miss}")

        return
    
    def execute_sex_check(self, sex_check: list = [0.2, 0.8]) -> None:
        """
        Execute sex check using PLINK to identify potential sex discrepancies in genetic data.
        This method performs sex check analysis by:
        1. Running PLINK's --check-sex command on pruned data
        2. Extracting X chromosome SNPs
        3. Calculating missingness rates for X chromosome SNPs
        
        Parameters
        ----------
        sex_check : list of float, default=[0.2, 0.8]
            List containing two float values that define the F-statistic boundaries for sex determination.
            The values must sum to 1.0. First value is the lower bound, second is the upper bound.
            Samples with F-statistics below the first value are called female, above the second value are called male.
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If sex_check is not a list or if its elements are not floats
        ValueError
            If sex_check doesn't contain exactly 2 elements or if they don't sum to 1
        
        Notes
        -----
        The method creates the following output files:
        - {output_name}-sexcheck.sexcheck : Contains sex check results
        - {output_name}-xchr.bed/bim/fam : X chromosome SNP data
        - {output_name}-xchr-missing.imiss : X chromosome missingness data
        The number of threads used is automatically determined based on available CPU cores,
        using max(available cores - 2, 1) or falling back to half of logical cores if CPU count
        cannot be determined.
        """

        if not isinstance(sex_check, list):
            raise TypeError("sex_check should be a list")
        if len(sex_check) != 2:
            raise ValueError("sex_check must have two elements")
        if not all(isinstance(i, float) for i in sex_check):
            raise TypeError("All elements in sex_check must be floats")
        if sum(sex_check) != 1:
            raise ValueError("The sum of sex_check elements must be equal to 1")
        
        logger.info(f"STEP: Check discordant sex information.")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

        plink_cmd1 = f"plink --bfile {self.pruned_file} --check-sex {sex_check[0]} {sex_check[1]} --threads {max_threads} --out {self.results_dir / (self.output_name+'-sexcheck')}"

        # extract xchr SNPs
        plink_cmd2 = f"plink --bfile {self.pruned_file} --chr 23 --keep-allele-order --make-bed --out {self.results_dir / (self.output_name+'-xchr')}"

        # run missingness on xchr SNPs
        plink_cmd3 = f"plink --bfile {self.results_dir / (self.output_name+'-xchr')} --missing --out {self.results_dir / (self.output_name+'-xchr-missing')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2, plink_cmd3]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self.sexcheck_miss = self.results_dir / (self.output_name + '-sexcheck.sexcheck')
        self.xchr_miss = self.results_dir / (self.output_name + '-xchr-missing.imiss')

        return

    def execute_heterozygosity_rate(self, maf: float = 0.01) -> None:
        """
        Executes heterozygosity rate analysis on genetic data using PLINK.

        This method performs a series of PLINK commands to analyze heterozygosity rates in genetic data,
        separating SNPs based on minor allele frequency (MAF) threshold and computing heterozygosity
        for both groups.

        Parameters
        ----------
        maf : float, optional
            Minor allele frequency threshold used to split SNPs into two groups.
            Must be between 0 and 0.5. Default is 0.01.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If maf is not a float
        ValueError
            If maf is not between 0 and 0.5
        FileNotFoundError
            If any of the expected output files are not created

        Notes
        -----
        The method:
        1. Extracts autosomal SNPs
        2. Splits SNPs based on MAF threshold
        3. Computes missingness
        4. Converts to PED/MAP format
        5. Computes heterozygosity for both MAF groups

        The computation uses optimized threading based on available CPU cores and memory.
        """

        if not isinstance(maf, float):
            raise TypeError("maf should be a float")
        if maf <= 0 or maf >= 0.5:
            raise ValueError("maf should be between 0 and 0.5")

        logger.info(f"STEP: Heterozygosity rate check. `maf` set to {maf}")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # extract autosomal SNPS
        plink_cmd1 = f"plink --bfile {self.pruned_file} --autosome --keep-allele-order --make-bed --out {self.results_dir / (self.output_name+'-chr1-22')}"

        # extract SNPs with minor allele frequency greater than threshold
        plink_cmd2 = f"plink --bfile {self.results_dir / (self.output_name+'-chr1-22')} --maf {maf} --keep-allele-order --make-bed --out {self.results_dir / (self.output_name+'-chr1-22-mafgreater')}"

        # extract SNPs with minor allele frequency less than threshold
        plink_cmd3 = f"plink --bfile {self.results_dir / (self.output_name+'-chr1-22')} --exclude {(self.results_dir / (self.output_name+'-chr1-22-mafgreater')).with_suffix('.bim')} --keep-allele-order --make-bed --out {self.results_dir / (self.output_name+'-chr1-22-mafless')}"

        # get missingness to plot against het
        plink_cmd4 = f"plink --bfile {self.results_dir / (self.output_name+'-chr1-22-mafgreater')} --missing --out {self.results_dir / (self.output_name+'-chr1-22-mafgreater-missing')}"
        plink_cmd5 = f"plink --bfile {self.results_dir / (self.output_name+'-chr1-22-mafless')} --missing --out {self.results_dir / (self.output_name+'-chr1-22-mafless-missing')}"

        # convert both to ped/map files for heterozigosity computation
        plink_cmd6 = f"plink --bfile {self.results_dir / (self.output_name+'-chr1-22-mafgreater')} --recode --out {self.results_dir / (self.output_name+'-chr1-22-mafgreater-recode')} --memory {memory} --threads {max_threads}"
        plink_cmd7 = f"plink --bfile {self.results_dir / (self.output_name+'-chr1-22-mafless')} --recode --out {self.results_dir / (self.output_name+'-chr1-22-mafless-recode')} --memory {memory} --threads {max_threads}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2, plink_cmd3, plink_cmd4, plink_cmd5, plink_cmd6, plink_cmd7]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self._compute_heterozigozity(
            ped_file=(self.results_dir / (self.output_name+'-chr1-22-mafgreater-recode')).with_suffix('.ped')
        )
        self._compute_heterozigozity(
            ped_file=(self.results_dir / (self.output_name+'-chr1-22-mafless-recode')).with_suffix('.ped')
        )

        self.summary_greater = self.results_dir / ('Summary-'+self.output_name+'-chr1-22-mafgreater-recode.ped')
        if not self.summary_greater.exists():
            raise FileNotFoundError(f"Missing file: {self.summary_greater}")
        self.summary_less    = self.results_dir / ('Summary-'+self.output_name+'-chr1-22-mafless-recode.ped')
        if not self.summary_less.exists():
            raise FileNotFoundError(f"Missing file: {self.summary_less}")
        self.maf_greater_miss= self.results_dir / (self.output_name+'-chr1-22-mafgreater-missing.imiss')
        if not self.maf_greater_miss.exists():
            raise FileNotFoundError(f"Missing file: {self.maf_greater_miss}")
        self.maf_less_miss   = self.results_dir / (self.output_name+'-chr1-22-mafless-missing.imiss')
        if not self.maf_less_miss.exists():
            raise FileNotFoundError(f"Missing file: {self.maf_less_miss}")

        return

    def execute_ibd(self) -> None:
        """
        Execute Identity by Descent (IBD) analysis using PLINK.

        This method performs duplicate and relatedness checks using IBD analysis. It runs two PLINK commands:
        1. Generates genome-wide IBD estimates
        2. Calculates missing genotype rates

        The method uses optimal thread count based on available CPU cores and validates input/output files.

        Returns:
        --------
            None

        Raises:
        -------
            FileNotFoundError: If required input pruned file is missing or if expected output files are not generated

        Required instance attributes:
            pruned_file: Path to pruned PLINK binary file
            results_dir: Directory path for output files
            output_name: Base name for output files
            ibd_miss: Path to missing genotype rate file (set by method)
            genome: Path to IBD estimates file (set by method)
        """

        logger.info("STEP: Duplicates and relatedness check with IBD")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

        if not self.pruned_file or not self.pruned_file.exists():
            raise FileNotFoundError(f"Missing file: {self.pruned_file}")

        # PLINK command
        plink_cmd1 = f"plink --bfile {self.pruned_file} --genome --out {self.results_dir / (self.output_name+'-ibd')} --threads {max_threads}"

        # PLINK command
        plink_cmd2 = f"plink --bfile {self.pruned_file} --allow-no-sex --missing --out {self.results_dir / (self.output_name+'-ibd-missing')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self.ibd_miss = self.results_dir / (self.output_name+'-ibd-missing.imiss')
        if not self.ibd_miss.exists():
            raise FileNotFoundError(f"Missing file: {self.ibd_miss}")
        self.genome = self.results_dir / (self.output_name+'-ibd.genome')
        if not self.genome.exists():
            raise FileNotFoundError(f"Missing file: {self.genome}")

        return

    def execute_kingship(self, kingship: float = 0.354) -> None:
        """
        Execute kinship analysis to identify and handle sample relatedness.
        This method performs kinship analysis using PLINK2 to identify duplicate samples and related individuals. 
        It first computes a kinship coefficient matrix for all samples and then prunes samples based on the 
        specified kingship threshold.
        
        Parameters
        ----------
        kingship : float, optional
            The kinship coefficient threshold used to identify related samples. Must be between 0 and 1.
            Samples with kinship coefficients above this threshold will be marked for removal.
            Default is 0.354 (equivalent to first-degree relatives).
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If kingship parameter is not a float.
        ValueError
            If kingship parameter is not between 0 and 1.
        FileNotFoundError
            If the expected output file from PLINK2 is not created.
        
        Notes
        -----
        - Uses PLINK2 to compute kinship coefficients and perform sample pruning
        - Automatically determines optimal thread count and memory usage based on system resources
        - Creates output files with kinship coefficient matrix and list of samples to be removed
        - Updates self.kinship_miss with path to file containing samples to be removed
        """

        if not isinstance(kingship, float):
            raise TypeError("kingship should be a float")
        if kingship < 0 or kingship >1:
            raise ValueError("kingship should be between 0 and 1")
        
        logger.info(f"STEP: Duplicates and relatedness check with Kingship. `kingship` set to {kingship}")
        
        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)
        
        # Compute kinship-coefficient matrix for all samples
        plink2_cmd1 = f"plink2 --bfile {self.pruned_file} --make-king triangle bin --out {self.results_dir / (self.output_name+'-kinship-coefficient-matrix')} --memory {memory} --threads {max_threads}"

        # Prune for Monozygotic Twins OR Duplicates
        plink2_cmd2 = f"plink2 --bfile {self.pruned_file} --king-cutoff {self.results_dir / (self.output_name+'-kinship-coefficient-matrix')} {kingship} --out {self.results_dir / (self.output_name+'-kinship-pruned-duplicates')} --memory {memory} --threads {max_threads}"

        # execute PLINK commands
        cmds = [plink2_cmd1, plink2_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        self.kinship_miss = (self.results_dir / (self.output_name+'-kinship-pruned-duplicates')).with_suffix('.king.cutoff.out.id')

        # Check if the file exists
        if not self.kinship_miss.exists():
            raise FileNotFoundError(f"Expected file {self.kinship_miss} was not created. Ensure the PLINK2 command executed successfully.")

        return
    
    def execute_duplicate_relatedness(self, kingship: float = 0.354, use_king: bool = True) -> None:
        """
        Execute duplicate and relatedness analysis on the genotype data.
        This method performs either IBD (Identity by Descent) or KING kinship coefficient
        analysis to identify duplicate samples and related individuals in the dataset.
        
        Parameters
        ----------
        kingship : float, optional
            The KING kinship coefficient threshold for identifying related samples.
            Default is 0.354, which corresponds to duplicates/MZ twins.
        use_king : bool, optional
            If True, uses KING algorithm for relatedness analysis.
            If False, uses traditional IBD analysis.
            Default is True.
        
        Returns
        -------
        None

        Raises
        ------
        TypeError
            If kingship is not a float or use_king is not a boolean.
        
        Notes
        -----
        The method will store the analysis type (KING or IBD) in the use_king attribute.
        """

        if not isinstance(use_king, bool):
            raise TypeError("use_king must be a boolean")
        if not isinstance(kingship, float):
            raise TypeError("kingship must be a float")

        logger.info("STEP: Duplicates and relatedness check")
        
        if use_king:
            self.execute_kingship(kingship)
        else:
            self.execute_ibd()

        self.use_king = use_king

        return

    def _compute_heterozigozity(self, ped_file: Path, map_file: Optional[Path] = None) -> None:
        """
        Computes heterozygosity statistics from a PED file and writes results to a summary file.
        This method analyzes a PED file to calculate homozygosity and heterozygosity rates
        for each individual. The results are written to a summary file in the same directory
        as the input PED file.
        
        Parameters
        ----------
        ped_file : Path
            Path to the input PED file containing genotype data
        map_file : Path, optional
            Path to the MAP file (not used in current implementation)
        
        Returns
        -------
        None
            Results are written to a summary file named "Summary-{ped_filename}"
        The output summary file contains the following columns:
        - ID: Individual identifier
        - total: Total number of valid genotypes
        - num_hom: Number of homozygous genotypes
        - num_het: Number of heterozygous genotypes  
        - Percent_hom: Percentage of homozygous genotypes
        - Percent_het: Percentage of heterozygous genotypes
        
        sNotes
        -----
        - Missing alleles (coded as '0' or 'N') are excluded from calculations
        - The method assumes PED file format with genotype data starting from column 7
        - Handles FileNotFound and IOError exceptions with appropriate error messages
        """
        
        # Define output file name
        summary_file= f"Summary-{ped_file.name}"
        output_path = ped_file.parent / summary_file

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
            logger.error(f"Error: {e}")

    def get_fail_samples(self, call_rate_thres: float, std_deviation_het: float, maf_het: float, ibd_threshold: float) -> pd.DataFrame:
        """
        Get samples that failed quality control checks and generate a summary DataFrame.
        This method performs multiple QC checks on samples:
        1. Call rate check
        2. Sex check
        3. Heterozygosity rate check (for MAF > and < threshold)
        4. Duplicates/Relatedness check (using either KING or IBD)
        
        Parameters
        ----------
        call_rate_thres : float
            Threshold for call rate filtering
        std_deviation_het : float 
            Number of standard deviations to use for heterozygosity filtering
        maf_het : float
            Minor allele frequency threshold for heterozygosity check
        ibd_threshold : float
            Threshold for IBD filtering (only used if use_king=False)
        
        Returns
        -------
        pd.DataFrame
            Summary DataFrame containing:
            - Counts of samples failing each QC check
            - Number of duplicated sample IDs
            - Total failures
        
        Raises
        ------
        FileNotFoundError
            If any required input files are missing
        TypeError
            If unexpected column types are found in summary DataFrame
        
        Notes
        -----
        The method saves a detailed fail_samples.txt file with all failed samples and their failure reasons.
        Samples failing multiple checks are only counted once in the final summary.
        """
        # Check if required files exist
        required_files = [
            self.call_rate_miss,
            self.results_dir / (self.output_name + '-sexcheck.sexcheck'),
            self.results_dir / (self.output_name + '-xchr-missing.imiss'),
            self.results_dir / ('Summary-' + self.output_name + '-chr1-22-mafgreater-recode.ped'),
            self.results_dir / (self.output_name + '-chr1-22-mafgreater-missing.imiss'),
            self.results_dir / ('Summary-' + self.output_name + '-chr1-22-mafless-recode.ped'),
            self.results_dir / (self.output_name + '-chr1-22-mafless-missing.imiss')
        ]

        if not self.use_king:
            required_files.append(self.results_dir / (self.output_name + '-ibd-missing.imiss'))
            required_files.append(self.results_dir / (self.output_name + '-ibd.genome'))

        for file in required_files:
            if not file.exists():
                raise FileNotFoundError(f"Required file not found: {file}")

        # ==========================================================================================================
        #                                             CALL RATE CHECK
        # ==========================================================================================================

        # load samples who failed call rate check
        fail_call_rate = self.report_call_rate(
            directory    =self.results_dir, 
            filename     =self.call_rate_miss,
            threshold    =call_rate_thres, 
            plots_dir    =self.plots_dir,
            y_axis_cap   =10
        )

        logger.info('Call rate check done')

        # ==========================================================================================================
        #                                             SEX CHECK
        # ==========================================================================================================

        fail_sexcheck = self.report_sex_check(
            directory          =self.results_dir, 
            sex_check_filename =self.output_name+'-sexcheck.sexcheck', 
            xchr_imiss_filename=self.output_name+'-xchr-missing.imiss',
            plots_dir          =self.plots_dir
        )

        logger.info('Sex check done')

        # ==========================================================================================================
        #                                       HETETROZYGOSITY RATE CHECK
        # ==========================================================================================================

        fail_het_greater = self.report_heterozygosity_rate(
            directory           = self.results_dir, 
            summary_ped_filename= 'Summary-'+self.output_name+'-chr1-22-mafgreater-recode.ped', 
            autosomal_filename  = self.output_name+'-chr1-22-mafgreater-missing.imiss', 
            std_deviation_het   = std_deviation_het,
            maf                 = maf_het,
            split               = '>',
            plots_dir           = self.plots_dir
        )

        logger.info(f'Heterozygosity rate check done for MAF > {maf_het}')

        fail_het_less = self.report_heterozygosity_rate(
            directory           = self.results_dir, 
            summary_ped_filename= 'Summary-'+self.output_name+'-chr1-22-mafless-recode.ped', 
            autosomal_filename  = self.output_name+'-chr1-22-mafless-missing.imiss', 
            std_deviation_het   = std_deviation_het,
            maf                 = maf_het,
            split               = '<',
            plots_dir           = self.plots_dir
        )

        logger.info(f'Heterozygosity rate check done for MAF < {maf_het}')

        # ==========================================================================================================
        #                                       DUPLICATES-RELATEDNESS CHECK
        # ==========================================================================================================

        if self.use_king:

            # load samples that failed duplicates and relatedness check
            duplicates_file = self.results_dir / (self.output_name+'-kinship-pruned-duplicates.king.cutoff.out.id')
            df_duplicates = pd.read_csv(
                duplicates_file,
                sep   =r'\s+',
                engine='python'
            )

            # filter samples that failed duplicates and relatedness check
            df_duplicates.columns = ['FID', 'IID']
            fail_duplicates = df_duplicates[['FID', 'IID']].reset_index(drop=True)
            fail_duplicates['Failure'] = 'Duplicates and relatedness (Kingship)'

            logger.info('Duplicates and relatedness check done with kingship')

        else:

            fail_duplicates = self.report_ibd_analysis(ibd_threshold)
            logger.info('Duplicates and relatedness check done with IBD')

        # ==========================================================================================================
        #                                       MERGE ALL FAILURES
        # ==========================================================================================================

        fails = [fail_call_rate, fail_sexcheck, fail_het_greater, fail_het_less, fail_duplicates] 

        df = pd.concat(fails, axis=0).reset_index(drop=True)

        summary = df['Failure'].value_counts().reset_index()
        num_dup = df.duplicated(subset=['FID', 'IID']).sum()

        df = df.drop_duplicates(subset=['FID', 'IID'])

        df.to_csv(self.fails_dir / 'fail_samples.txt', index=False, sep='\t')

        totals = summary.select_dtypes(include="number").sum() - num_dup
        dups_row = pd.DataFrame({summary.columns[0]: ['Duplicated Sample IDs'], summary.columns[1]: [-num_dup]})
        # Create the total row
        dups_row = pd.DataFrame({'Failure':['Duplicated Sample IDs'], 'count':[-num_dup]})
        # Validate column types in the summary DataFrame
        for col in summary.columns:
            if not (pd.api.types.is_numeric_dtype(summary[col]) or summary[col].dtype == 'object'):
                raise TypeError(f"Unexpected column type in summary DataFrame: {col} has type {summary[col].dtype}")

        # Construct the total_row DataFrame
        total_row = pd.DataFrame({col: [totals[col] if col in totals.index else "Total"] for col in summary.columns})

        # Append the total row to the DataFrame
        summary = pd.concat([summary, dups_row, total_row], ignore_index=True)
        
        return summary
    
    def execute_drop_samples(self) -> None:
        """
        Execute the removal of samples that failed quality control checks using PLINK.
        This method performs the following steps:
        1. Determines the appropriate binary file name based on previous processing steps
        2. Reads the fail_samples.txt file containing samples to be removed
        3. Executes PLINK command to create new binary files excluding failed samples
        
        Raises:
        -------
            FileNotFoundError: If the fail_samples.txt file is not found in the fails directory

        Returns:
        --------
            None

        Notes:
        ------
            - The output files will be created with suffix '-clean-samples'
            - The method preserves allele order during the operation
            - Input files must be in PLINK binary format (.bed, .bim, .fam)
        """
        
        logger.info("STEP: Drop samples that failed quality control checks")

        if self.hh_to_missing:
            binary_name = self.input_name+'-hh-missing'
        elif self.renamed_snps:
            binary_name = self.input_name+'-renamed'
        else:
            binary_name = self.input_name

        logger.info(f"Binary file name: {binary_name}")

        # drop samples
        fail_samples_file = self.fails_dir / 'fail_samples.txt'
        if not fail_samples_file.exists():
            raise FileNotFoundError(f"Required file {fail_samples_file} not found. Ensure the fail_samples.txt file is generated before executing this step.")

        plink_cmd = f"plink --bfile {self.input_path / binary_name} --remove {fail_samples_file} --keep-allele-order --make-bed --out {self.clean_dir / (self.output_name+'-clean-samples')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return
  
    def report_call_rate(self, directory: Path, filename: Path, threshold: float, plots_dir: Optional[Path] = None, y_axis_cap: int = 10, color: str = '#1B9E77', line_color: str = '#D95F02') -> pd.DataFrame:
        """
        Generate sample call rate analysis plots and identify samples failing the call rate threshold.
        This method reads a PLINK-format missing rate file, creates visualization plots, and identifies
        samples that fail the specified call rate threshold. It generates two sets of plots:
        1. Histograms showing the distribution of missing SNPs (F_MISS)
        2. Scatterplots showing different views of the call rate data
        
        Parameters
        ----------
        directory : Path
            Directory containing the input file
        filename : str
            Name of the PLINK format missing rate file
        threshold : float
            Call rate threshold for sample filtering (in terms of F_MISS)
        plots_dir : Path, optional
            Directory where plots will be saved. If None, uses default plots directory
        y_axis_cap : int, optional
            Maximum value for y-axis in capped histogram plots. Default is 10
        color : str, optional
            Color for the main plot elements. Default is '#1B9E77'
        line_color : str, optional
            Color for threshold lines in plots. Default is '#D95F02'
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing samples that failed the call rate threshold with columns:
            - FID: Family ID
            - IID: Individual ID
            - Failure: Always set to 'Call rate'
        
        Notes
        -----
        The method generates two JPEG files:
        - call_rate_{threshold}_histogram.jpeg: Contains histogram plots
        - call_rate_{threshold}_scatterplot.jpeg: Contains scatter plots
        """
        
        if not plots_dir:
            plots_dir = self.plots_dir

        # load samples that failed sex check
        df_call_rate = pd.read_csv(
            directory / filename,
            sep=r'\s+',
            engine='python'
        )

        # filter samples that fail call rate
        fail_call_rate = df_call_rate[df_call_rate['F_MISS'] > threshold][['FID', 'IID']].reset_index(drop=True)
        fail_call_rate['Failure'] = 'Call rate'

        # Create the figure and subplots
        fig1, axes1 = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

        # First subplot: Full histogram
        axes1[0] = sns.histplot(df_call_rate['F_MISS'], bins=30, color=color, alpha=0.7, ax=axes1[0]) # type: ignore
        axes1[0].set_title("Sample Call Rate Distribution")
        axes1[0].set_xlabel("Proportion of missing SNPs (F_MISS)")
        axes1[0].set_ylabel("Frequency")

        # Second subplot: Histogram with capped y-axis
        axes1[1] = sns.histplot(df_call_rate['F_MISS'], bins=30, color=color, alpha=0.7, ax=axes1[1]) # type: ignore
        axes1[1].set_ylim(0, y_axis_cap)  # Cap y-axis
        axes1[1].set_title("Sample Call Rate Distribution (Capped)")
        axes1[1].set_xlabel("Proportion of missing SNPs (F_MISS)")

        plt.tight_layout()
        plt.savefig(plots_dir / f"call_rate_{threshold}_histogram.jpeg", dpi=400)
        plt.show(block=False)

        fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

        # First subplot: capped y-axis
        axes2[0] = sns.histplot(df_call_rate['F_MISS'], bins=50, color=color, alpha=0.7, ax=axes2[0]) # type: ignore
        axes2[0].set_ylim(0, y_axis_cap)  # Cap y-axis
        axes2[0].set_title("Sample Call Rate Distribution (Capped)")
        axes2[0].set_xlabel("Proportion of missing SNPs (F_MISS)")

        # Add a vertical line at the threshold
        axes2[0].axvline(threshold, linewidth=2, color=line_color, linestyle='dashed')

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
            edgecolor='none',
            color =color,
            ax    =axes2[1]
        ) 
        axes2[1].set_title("Sample Call Rate")
        axes2[1].set_xlabel(f"Number of samples")
        axes2[1].set_ylabel("F_MISS")

        # Add a vertical line at the threshold
        axes2[1].axhline(threshold, linewidth=2, color=line_color, linestyle='dashed')

        # third subplot: Number of samples vs F_MISS
        axes2[2] = sns.scatterplot(
            x      =df_call_rate['F_MISS'],
            y      =np.random.normal(size=len(df_call_rate['F_MISS'])),
            markers='o',
            s      =20,
            color =color,
        )
        axes2[2].set_title("Sample Call Rate")
        axes2[2].set_xlabel("Proportion of missing SNPs (F_MISS)")
        axes2[2].set_ylabel(f"Samples")
        axes2[2].set_yticks([])
    

        # Add a vertical line at the threshold
        axes2[2].axvline(threshold, linewidth=2, color=line_color, linestyle='dashed')

        plt.tight_layout()
        plt.savefig(plots_dir / f"call_rate_{threshold}_scatterplot.jpeg", dpi=400)
        plt.show(block=False)

        return fail_call_rate
    
    def report_sex_check(self, directory: Path, sex_check_filename: str, xchr_imiss_filename: str, plots_dir: Optional[Path] = None) -> pd.DataFrame:
        """
        Creates a sex check report and visualization based on PLINK's sex check results.
        This function reads sex check data and X chromosome missingness data, merges them,
        and generates a scatter plot to visualize potential sex discrepancies. It also identifies
        samples that fail sex check quality control.
        
        Parameters
        ----------
        directory : Path
            Path to the directory containing input files
        sex_check_filename : str
            Filename of PLINK's sex check results (typically .sexcheck file)
        xchr_imiss_filename : str
            Filename of X chromosome missingness data
        plots_dir : Path, optional
            Directory where the plot should be saved. If None, uses default plots directory
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing samples that failed sex check QC with columns:
            - FID: Family ID
            - IID: Individual ID
            - Failure: Type of failure (always 'Sex check')
        
        Notes
        -----
        The function creates a scatter plot with:
        - Blue hollow circles for samples with Male PEDSEX
        - Green hollow circles for samples with Female PEDSEX
        - Red filled circles for problematic samples
        - Dotted red vertical lines at F=0.2 and F=0.8
        The plot is saved as 'sex_check.jpeg' in the specified plots directory.
        """
        
        if not plots_dir:
            plots_dir = self.plots_dir

        df_sexcheck = pd.read_csv(
            directory / sex_check_filename,
            sep   =r'\s+',
            engine='python'
        )

        df_xchr_imiss = pd.read_csv(
            directory / xchr_imiss_filename,
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
                edgecolors=palette[category],     # Map color # type: ignore
                facecolors='none',                # Hollow circles
                s         =size_mapping[category],# Map size # type: ignore
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
        plt.savefig(plots_dir / 'sex_check.jpeg', dpi=400)

        return fail_sexcheck
    
    def report_heterozygosity_rate(self, directory: Path, summary_ped_filename: str, autosomal_filename: str, std_deviation_het: float, maf: float, split: str, plots_dir: Path, y_axis_cap: float = 80) -> pd.DataFrame:
        """
        Analyze and report heterozygosity rates for samples, creating visualization plots and identifying samples that fail heterozygosity rate checks.
        This function loads heterozygosity and autosomal call rate data, merges them, identifies samples with deviant heterozygosity rates,
        and generates visualization plots to aid in quality control analysis.
        
        Parameters
        ----------
        directory : str
            Path to the directory containing input files
        summary_ped_filename : str
            Filename of the summary PED file containing heterozygosity information
        autosomal_filename : str
            Filename of the autosomal file containing call rate information
        std_deviation_het : float
            Number of standard deviations to use as threshold for identifying deviant samples
        maf : float
            Minor allele frequency threshold used in the analysis
        split : str
            Direction of MAF comparison ('>' or '<')
        plots_dir : str
            Directory where plot files will be saved
        y_axis_cap : float, optional
            Maximum value for y-axis in capped histogram plot (default: 80)
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing samples that failed heterozygosity rate check with columns:
            - FID: Family ID
            - IID: Individual ID
            - Failure: Description of failure type
        
        Notes
        -----
        The function generates two types of plots:
        1. Histograms of heterozygosity rates (both uncapped and capped)
        2. Scatter plot of heterozygosity rate vs missing SNP proportion
        Files are saved as JPEG images in the specified plots directory.
        """
        
        # load samples that failed heterozygosity rate check with MAF > threshold
        maf_file = directory / summary_ped_filename
        df_maf = pd.read_csv(
            maf_file,
            sep   =r'\s+',
            engine='python'
        )

        # autosomal call rate per individual
        autosomal_file = directory / autosomal_filename
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

        axes1[0] = sns.histplot(df_het['Percent_het'], bins=30, color='green', alpha=0.7, ax=axes1[0]) # type: ignore
        axes1[0].set_title("Autosomal heterozygosity")
        axes1[0].set_xlabel(f"% Heterozygosity MAF {split} {maf}")
        axes1[0].set_ylabel("Frequency")

        axes1[1] = sns.histplot(df_het['Percent_het'], bins=30, color='green', alpha=0.7, ax=axes1[1]) # type: ignore
        axes1[1].set_title("Autosomal heterozygosity (capped)")
        axes1[1].set_xlabel(f"% Heterozygosity MAF {split} {maf}")
        axes1[1].set_ylim(0, y_axis_cap)  # Cap y-axis
        axes1[1].set_ylabel("Frequency")

        plt.tight_layout()
        
        if split == '>':
            plt.savefig(plots_dir / f"heterozygosity_rate_greater_{maf}_histogram.jpeg", dpi=400)
        else:
            plt.savefig(plots_dir / f"heterozygosity_rate_less_{maf}_histogram.jpeg", dpi=400)
        
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
            plt.savefig(plots_dir / f"heterozygosity_rate_greater_{maf}_scatterplot.jpeg", dpi=400)
        else:
            plt.savefig(plots_dir / f"heterozygosity_rate_less_{maf}_scatterplot.jpeg", dpi=400)
        plt.show(block=False)

        return fail_het

    def report_ibd_analysis(self, ibd_threshold: float = 0.185, chunk_size: int = 100000) -> pd.DataFrame:
        """
        Analyze IBD (Identity By Descent) to identify duplicated or related samples.
        This method processes IBD analysis results to identify sample pairs with IBD scores
        above a specified threshold, indicating potential duplicates or related individuals.
        For identified pairs, it uses missingness data to determine which sample should be
        removed (keeping the sample with lower missingness rate).
        
        Parameters
        ----------
        ibd_threshold : float, default=0.185
            The PI_HAT threshold above which samples are considered related.
            Typical values: >0.98 for duplicates, >0.5 for first-degree relatives.
        chunk_size : int, default=100000
            Number of rows to process at a time when reading the genome file.
        
        Returns
        -------
        pd.DataFrame
            A DataFrame containing samples to be removed, with columns:
            - FID: Family ID
            - IID: Individual ID
            - Failure: Reason for removal ('Duplicates and relatedness (IBD)')
            Returns empty DataFrame if no related pairs are found or if KING is used.
        
        Raises
        ------
        TypeError
            If ibd_threshold is not a float.
        FileNotFoundError
            If required input files (*.imiss or *.genome) are not found.
        
        Notes
        -----
        The method requires two input files:
        - {output_name}-ibd-missing.imiss: Contains sample missingness information
        - {output_name}-ibd.genome: Contains pairwise IBD estimates
        For each related pair, the sample with higher missingness rate is marked for removal.
        """
        
        if not isinstance(ibd_threshold, float):
            raise TypeError("ibd_threshold should be a float")

        if self.use_king:
            return pd.DataFrame()

        # File paths

        imiss_path = self.results_dir / (self.output_name + '-ibd-missing.imiss')
        genome_path= self.results_dir / (self.output_name + '-ibd.genome')

        if not imiss_path.exists():
            raise FileNotFoundError(f"Missing file: {imiss_path}")
        if not genome_path.exists():
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
    
    def clean_input_folder(self) -> None:
        """Removes specific files from the input folder.

        This method cleans the input folder by removing files that contain 'missing' or
        'renamed' in their names, excluding log files. The cleaning process helps maintain
        organization by removing temporary or processed files.

        Returns:
            None
        """

        logger.info("STEP: Clean input folder")

        # remove all files from input folder
        for file in self.input_path.glob('*'):
            if file.is_file() and file.suffix != '.log' and 'missing' in file.name:
                file.unlink()
            elif file.is_file() and file.suffix != '.log' and 'renamed' in file.name:
                file.unlink()

        return
    
    def clean_result_folder(self) -> None:
        """
        Clean the results folder by removing specific file types.

        This method removes all .bed, .bim, and .fam files from the results directory,
        while preserving log files. The cleaning process helps maintain a tidy workspace
        by removing intermediate or temporary files from previous runs.

        Returns:
            None
        """

        logger.info("STEP: Clean results folder")

        # remove all files from output folder
        for file in self.results_dir.glob('*'):
            if file.is_file() and file.suffix != '.log' and '.bed' in file.name:
                file.unlink()
            elif file.is_file() and file.suffix != '.log' and '.bim' in file.name:
                file.unlink()
            elif file.is_file() and file.suffix != '.log' and '.fam' in file.name:
                file.unlink()

        return
