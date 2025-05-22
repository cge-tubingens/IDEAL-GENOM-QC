import os
import psutil
import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from typing import Union, Literal

from ideal_genom_qc.Helpers import shell_do, delete_temp_files
from ideal_genom_qc.get_references import Fetcher1000Genome, FetcherLDRegions

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class ReferenceGenomicMerger():

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name:str, high_ld_regions: Path, reference_files: dict, built: str = '38') -> None:
        """
        Initialize AncestryQC class.
        This class performs ancestry quality control on genetic data by comparing study samples against reference populations.
        
        Parameters:
        -----------
        input_path: Path 
            Path to directory containing input files
        input_name: str 
            Name of input file without extension
        output_path: Path 
            Path to directory for output files
        output_name: str 
            Name for output files without extension
        high_ld_regions: Path 
            Path to file containing high LD regions to exclude
        reference_files: dict 
            Dictionary containing paths to reference population files
        built: str (optional) 
            Genome build version ('37' or '38'). Defaults to '38'
        
        Raises:
        -------
        TypeError: 
            If input arguments are not of correct type
        ValueError: 
            If genome build version is not '37' or '38'
        FileNotFoundError: 
            If required input files/directories do not exist
        
        Attributes:
        -----------
            reference_AC_GT_filtered: Filtered reference allele counts and genotypes
            study_AC_GT_filtered: Filtered study allele counts and genotypes
            pruned_reference: LD-pruned reference data
            pruned_study: LD-pruned study data
            reference_fixed_chr: Reference data with fixed chromosomes
            reference_fixed_pos: Reference data with fixed positions
            reference_flipped: Reference data with flipped alleles
            reference_cleaned: Final cleaned reference data
        """

        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(high_ld_regions, Path):
            raise TypeError("high_ld_regions should be a Path object")
        if not isinstance(reference_files, dict):
            raise TypeError("reference_files should be a dictionary")
        if not isinstance(input_name, str):
            raise TypeError("input_name should be a string")
        if not isinstance(output_name, str):
            raise TypeError("output_name should be a string")
        if not isinstance(built, str):
            raise TypeError("built should be a string")
        if built not in ['37', '38']:
            raise ValueError("built should be either '37' or '38'")
        
        if not input_path.exists():
            raise FileNotFoundError(f"input_path does not exist: {input_path}")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")
        if not high_ld_regions.exists():
            raise FileNotFoundError("high_ld_regions does not exist")

        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.output_name= output_name
        self.high_ld_regions = high_ld_regions
        self.reference_files = reference_files

        self.reference_AC_GT_filtered= None
        self.study_AC_GT_filtered    = None
        self.pruned_reference        = None
        self.pruned_study            = None
        self.reference_fixed_chr     = None
        self.reference_fixed_pos     = None
        self.reference_flipped       = None
        self.reference_cleaned       = None

        pass

    def execute_rename_snpid(self) -> None:
        
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

        logger.info("STEP: Renaming SNP IDs in the study data using PLINK2")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

        plink2_cmd = f"plink2 --bfile {self.input_path / self.input_name} --set-all-var-ids @:#:$r:$a --threads {max_threads} --make-bed --out {self.output_path / (self.input_name+ '-renamed')}"

        # Execute PLINK2 command
        shell_do(plink2_cmd, log=True)

        return

    def execute_filter_prob_snps(self)->None:
        """
        Executes the filtering of problematic SNPs (A->T and C->G) from both study and reference data.
        This method performs the following operations:
        1. Identifies and filters A->T and C->G SNPs from study data
        2. Identifies and filters A->T and C->G SNPs from reference data
        3. Creates new PLINK binary files excluding the identified problematic SNPs
        4. Uses maximum available CPU threads (total cores - 2) and 2/3 of available memory
        The method handles both renamed and original SNP scenarios, determined by self.renamed_snps.
        
        Returns:
        --------
            None

        Side Effects:
        -------------
            - Creates filtered SNP list files in the output directory
            - Creates new PLINK binary files (.bed, .bim, .fam) in the output directory
            - Sets self.reference_AC_GT_filtered and self.study_AC_GT_filtered paths
            - Logs progress and statistics of filtering operations

        Requires:
        ---------
            - Valid PLINK binary files for both study and reference data
            - Proper initialization of input_path, output_path, and reference_files
        """

        logger.info("STEP: Filtering A->T and C->G SNPs from study and reference data.")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = cpu_count-2
        else:
            max_threads = 10
        
        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)


        # find A->T and C->G SNPs in study data
        filtered_study = self._filter_non_AT_or_GC_snps(target_bim=self.output_path / f"{self.input_name}-renamed.bim", output_filename=self.input_name)
        logger.info("STEP: Filtering problematic SNPs from the study data: filtered study data")

        # find A->T and C->G SNPs in reference data
        filtered_reference = self._filter_non_AT_or_GC_snps(target_bim=self.reference_files['bim'], output_filename=self.reference_files['bim'].stem)
        logger.info("STEP: Filtering problematic SNPs from the study data: filtered reference data")

        self.reference_AC_GT_filtered= self.output_path / f"{self.reference_files['bim'].stem}-no_ac_gt_snps"
        self.study_AC_GT_filtered    = self.output_path / f"{self.input_name}-no_ac_gt_snps"

        with open(filtered_study, 'r') as f:
            logger.info(f"STEP: Filtering problematic SNPs from the study data: {len(f.readlines())} SNPs filtered")
        with open(filtered_reference, 'r') as f:
            logger.info(f"STEP: Filtering problematic SNPs from the reference data: {len(f.readlines())} SNPs filtered")


        # PLINK command: generate cleaned study data files
        plink_cmd1 = f"plink --bfile  {self.output_path / (self.input_name+'-renamed')} --chr 1-22 --exclude {filtered_study} --keep-allele-order --threads {max_threads} --make-bed --out {self.study_AC_GT_filtered}"

        # Make sure the reference bim path is valid and extract the base filename
        if not self.reference_files.get('bim') or not isinstance(self.reference_files['bim'], Path):
            raise ValueError("reference_files dictionary must contain a valid 'bim' Path")
        
        reference_base = self.reference_files['bim'].with_suffix('')
        
        # PLINK command: generate cleaned reference data files
        plink_cmd2 = f"plink --bfile {reference_base} --biallelic-only strict --chr 1-22 --exclude {filtered_reference} --keep-allele-order --allow-extra-chr --memory {memory} --threads {max_threads} --make-bed --out {self.reference_AC_GT_filtered}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        return
    
    def execute_ld_pruning(self, ind_pair:list) -> None:
        """
        Execute linkage disequilibrium (LD) pruning on study and reference data.
        
        This method performs LD-based pruning using PLINK to remove highly correlated SNPs 
        from both study and reference datasets. The pruning is done using a sliding window 
        approach where SNPs are removed based on their pairwise correlation (r²).
        
        Parameters
        ----------
        ind_pair : list
            A list containing three elements:
            
            - ind_pair[0] (int): Window size in SNPs  
            - ind_pair[1] (int): Number of SNPs to shift the window at each step  
            - ind_pair[2] (float): r² threshold for pruning
        
        Raises
        ------
        TypeError
            If ind_pair is not a list.
        TypeError
            If first two elements of ind_pair are not integers.
        TypeError
            If third element of ind_pair is not a float.
        
        Returns
        -------
        None
        
        Notes
        -----
        - Uses PLINK's `--indep-pairwise` command for pruning.
        - Excludes high LD regions specified in `self.high_ld_regions`.
        - Creates pruned datasets for both study and reference data.
        - Updates `self.pruned_reference` and `self.pruned_study` with paths to pruned files.
        - Uses all available CPU threads except 2 for processing.
        """


        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        if not isinstance(ind_pair[0], int) or not isinstance(ind_pair[1], int):
            raise TypeError("The first two elements in ind_pair values should be integers (windows size and step size)")
        
        if not isinstance(ind_pair[2], float):
            raise TypeError("The third element in ind_pair should be a float (r^2 threshold)")
        
        logger.info("STEP: LD-based pruning of study and reference data")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = cpu_count-2
        else:
            max_threads = 10

        # PLINK command: generates prune.in and prune.out files from study data
        plink_cmd1 = f"plink --bfile {str(self.study_AC_GT_filtered)} --exclude range {self.high_ld_regions} --keep-allele-order --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --threads {max_threads} --out {str(self.output_path / self.input_name)}"

        # PLINK command: prune study data and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {str(self.study_AC_GT_filtered)} --extract {str((self.output_path / self.input_name).with_suffix('.prune.in'))} --keep-allele-order --threads {max_threads} --make-bed --out {str((self.output_path / (self.input_name+'-pruned')))}"

        # PLINK command: generates a pruned reference data files
        plink_cmd3 = f"plink --bfile {str(self.reference_AC_GT_filtered)} --extract {str((self.output_path / self.input_name).with_suffix('.prune.in'))} --keep-allele-order --make-bed --threads {max_threads} --out {str((self.output_path / (self.reference_files['bim'].stem+'-pruned')))}"

        self.pruned_reference = self.output_path / (self.reference_files['bim'].stem+'-pruned')
        self.pruned_study = self.output_path / (self.input_name+'-pruned')

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2, plink_cmd3]
        for cmd in cmds:
            shell_do(cmd, log=True)

        return
    
    def execute_fix_chromosome_mismatch(self) -> None:
        """
        Fix chromosome mismatch between study data and reference panel.

        This method executes PLINK commands to correct any chromosome mismatches between the study data
        and reference panel datasets. It identifies mismatches using internal methods and updates
        the chromosome assignments in the reference panel to match the study data.

        The method performs the following steps:
        1. Identifies chromosome mismatches between study and reference BIM files
        2. Creates an update file for chromosome reassignment
        3. Executes PLINK command to update chromosome assignments in reference panel

        Returns
        -------
        None

        Notes
        -----
        - Creates new PLINK binary files with updated chromosome assignments
        - The updated files are saved with '-updateChr' suffix
        
        Raises
        ------
        ValueError
            If pruned_study or pruned_reference is None, meaning execute_ld_pruning() was not called first
        """

        logger.info("STEP: Fixing chromosome mismatch between study data and reference panel")

        # Check if pruned_study and pruned_reference have been set
        if self.pruned_study is None:
            raise ValueError("pruned_study is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")
        if self.pruned_reference is None:
            raise ValueError("pruned_reference is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = cpu_count-2
        else:
            max_threads = 10

        # File paths - using with_suffix instead of with_name for more reliability
        study_bim = self.pruned_study.with_suffix('.bim')
        reference_bim = self.pruned_reference.with_suffix('.bim')

        to_update_chr_file = self._find_chromosome_mismatch(study_bim, reference_bim)

        self.reference_fixed_chr = self.output_path / f"{self.reference_files['bim'].stem}-updateChr"

        with open(to_update_chr_file, 'r') as f:
            logger.info(f"STEP: Fixing chromosome mismatch between study data and reference panel: {len(f.readlines())} SNPs to update")

        # PLINK command
        plink_cmd = f"plink --bfile {self.pruned_reference} --update-chr {to_update_chr_file} 1 2 --keep-allele-order --threads {max_threads} --make-bed --out {self.reference_fixed_chr}"

        # Execute PLINK command
        shell_do(plink_cmd, log=True)

        return
    
    def execute_fix_possition_mismatch(self) -> None:
        """
        Fixes position mismatches between study data and reference panel.

        This method executes PLINK commands to update the positions of SNPs in the reference panel
        to match those in the study data. It processes previously identified position mismatches
        and creates new binary PLINK files with corrected positions.

        The method:
        1. Determines optimal thread count for processing
        2. Identifies position mismatches between study and reference BIM files
        3. Updates reference panel positions using PLINK
        4. Creates new binary files with corrected positions

        Returns:
        --------
            None

        Side Effects:
        -------------
            - Creates new PLINK binary files (.bed, .bim, .fam) with updated positions
            - Logs the number of SNPs being updated
            - Modifies self.reference_fixed_pos with path to updated files

        Dependencies:
        -------------
            - Requires PLINK to be installed and accessible
            - Expects pruned study and reference files to exist
            - Requires previous chromosome fixing step to be completed
        """

        logger.info("STEP: Fixing position mismatch between study data and reference panel")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

        # Check if pruned_study and reference_fixed_chr have been properly set
        if self.pruned_study is None:
            raise ValueError("pruned_study is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")
        if self.reference_fixed_chr is None:
            raise ValueError("reference_fixed_chr is not set. Make sure execute_fix_chromosome_mismatch() is called before this method.")
            
        # File paths - using with_suffix instead of with_name for more reliability
        study_bim = self.pruned_study.with_suffix('.bim')
        reference_bim = self.reference_fixed_chr.with_suffix('.bim')

        to_update_pos_file = self._find_position_mismatch(study_bim, reference_bim)

        self.reference_fixed_pos = self.output_path / f"{self.reference_files['bim'].stem}-updatePos"

        with open(to_update_pos_file, 'r') as f:
            logger.info(f"STEP: Fixing position mismatch between study data and reference panel: {len(f.readlines())} SNPs to update")

        # PLINK command
        plink_cmd = f"plink --bfile {self.reference_fixed_chr} --update-map {to_update_pos_file} --keep-allele-order --threads {max_threads} --make-bed --out {self.reference_fixed_pos}"

        # Execute PLINK command
        shell_do(plink_cmd, log=True)

        return
    
    def execute_fix_allele_flip(self) -> None:
            """
            Executes the allele flipping process between study data and reference panel.
    
            This method performs the following steps:
            1. Identifies SNPs requiring allele flipping between study and reference data
            2. Creates a list of SNPs to flip
            3. Generates a new reference panel with flipped alleles using PLINK
    
            The method uses multi-threading capabilities based on available CPU cores,
            reserving 2 cores for system processes when possible.
    
            Returns:
            --------
                None
    
            Side Effects:
            -------------
                - Creates a .toFlip file containing SNPs requiring allele flipping
                - Generates new PLINK binary files (.bed, .bim, .fam) with flipped alleles
                - Logs the number of SNPs requiring flipping
                - Updates self.reference_flipped with the path to new flipped reference files
    
            Dependencies:
            -------------
                - PLINK must be installed and accessible in system PATH
                - Requires valid PLINK binary files for both study and reference data
                - Requires write permissions in output directory
            """
    
            logger.info("STEP: Allele flipping between study data and reference panel")
    
            cpu_count = os.cpu_count()
            if cpu_count is not None:
                max_threads = max(1, cpu_count - 2)
            else:
                max_threads = 10
                
            # Check if pruned_study and reference_fixed_pos have been properly set
            if self.pruned_study is None:
                raise ValueError("pruned_study is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")
            if self.reference_fixed_pos is None:
                raise ValueError("reference_fixed_pos is not set. Make sure execute_fix_possition_mismatch() is called before this method.")
                
            # File paths - using with_suffix for consistency and reliability
            study_bim = self.pruned_study.with_suffix('.bim')
            reference_bim = self.reference_fixed_pos.with_suffix('.bim')
    
            to_flip_file = self.output_path / f"{self.reference_files['bim'].stem}.toFlip"
            self._find_allele_flip(study_bim, reference_bim, to_flip_file)
    
            self.reference_flipped = self.output_path / f"{self.reference_files['bim'].stem}-flipped"
    
            with open(to_flip_file, 'r') as f:
                logger.info(f"STEP: Allele flipping between study data and reference panel: {len(f.readlines())} SNPs to flip")
    
            # plink command
            plink_cmd = f"plink --bfile {self.reference_fixed_pos} --flip {to_flip_file} --keep-allele-order --threads {max_threads} --make-bed --out {self.reference_flipped}"
    
            # execute PLINK command
            shell_do(plink_cmd, log=True)
    
            return

    def execute_remove_mismatches(self) -> None:
        """
        Removes mismatched SNPs from the reference data based on allele comparisons between study and reference datasets.

        This method performs the following steps:
        1. Determines optimal thread count for processing
        2. Identifies allele mismatches between study and reference BIM files
        3. Creates a list of SNPs to remove
        4. Generates a cleaned reference dataset excluding mismatched SNPs

        The method utilizes PLINK to perform the actual SNP removal while maintaining allele order.

        Returns:
        --------
            None

        Side Effects:
        -------------
            - Creates a file listing SNPs to be removed at {output_path}/{reference_bim_stem}.toRemove
            - Generates cleaned reference files at {output_path}/{reference_bim_stem}-cleaned.bed/bim/fam
            - Logs the number of SNPs being removed
        """

        logger.info("STEP: Removing mismatched SNPs from reference data")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

        if self.pruned_study is None:
            raise ValueError("pruned_study is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")
        if self.pruned_reference is None:
            raise ValueError("pruned_reference is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")

        # File paths - using with_suffix instead of with_name for more reliability
        study_bim = self.pruned_study.with_suffix('.bim')
        reference_bim = self.pruned_reference.with_suffix('.bim')

        mismatches_file = self.output_path / f"{self.reference_files['bim'].stem}.toRemove"
        self._find_allele_flip(study_bim, reference_bim, mismatches_file)

        self.reference_cleaned = self.output_path / f"{self.reference_files['bim'].stem}-cleaned"

        with open(mismatches_file, 'r') as f:
            logger.info(f"STEP: Removing mismatched SNPs from reference data: {len(f.readlines())} SNPs to remove")

        # plink command
        plink_cmd = f"plink --bfile {self.reference_flipped} --exclude {mismatches_file} --keep-allele-order --threads {max_threads} --make-bed --out {self.reference_cleaned}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return
    
    def execute_merge_data(self) -> None:
        """
        Merge study and reference data using PLINK.

        This method merges the pruned study data with the cleaned reference data using PLINK's
        --bmerge functionality. It automatically determines the optimal number of threads to use
        based on available CPU cores.

        The method:
        1. Calculates optimal thread count (CPU count - 2 or half of available cores)
        2. Constructs PLINK command for merging datasets
        3. Executes the merge operation via shell command

        Returns:
        --------
            None

        Side effects:
        -------------
            - Creates merged PLINK binary files (.bed, .bim, .fam) in the output directory
            - Logs the merge operation
        """

        logger.info("STEP: Merging study and reference data")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = cpu_count-2
        else:
            max_threads = 10

        if self.reference_cleaned is None:
            raise ValueError("reference_cleaned is not set. Make sure execute_remove_mismatches() is called before this method and completed successfully.")

        # plink command
        plink_cmd = f"plink --bfile {self.pruned_study} --bmerge {str(self.reference_cleaned.with_suffix('.bed'))} {str(self.reference_cleaned.with_suffix('.bim'))} {str(self.reference_cleaned.with_suffix('.fam'))} --keep-allele-order --threads {max_threads} --make-bed --out {self.output_path / (self.output_name+'-merged')}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        return

    def _filter_non_AT_or_GC_snps(self, target_bim: Path, output_filename: str) -> Path:
        """
        Filter SNPs that are not A/T or G/C variants from a PLINK BIM file.
        This method reads a BIM file and identifies SNPs that are either A/T or G/C variants.
        These variants are known as strand-ambiguous SNPs because their complementary alleles 
        are the same as their original alleles, making it impossible to determine the correct 
        strand without additional information.
        
        Parameters
        ----------
        target_bim : Path
            Path to the input BIM file containing SNP information
        output_filename : str
            Base name for the output file (without extension)
        
        Returns
        -------
        Path
            Path to the output file containing filtered SNP IDs with .ac_get_snps extension
        
        Notes
        -----
        The input BIM file should be tab-delimited with standard PLINK BIM format.
        Only columns containing SNP ID (column 2) and alleles (columns 5 and 6) are used.
        """

        df = pd.read_csv(
            target_bim, sep="\t", header=None, usecols=[1, 4, 5], names=["SNP", "A1", "A2"]
        )

        output_file = self.output_path / f"{output_filename}.ac_get_snps"

        filtered_snps = df[df[['A1', 'A2']].apply(lambda x: ''.join(sorted(x)) in {"AT", "TA", "GC", "CG"}, axis=1)]
        
        filtered_snps[["SNP"]].to_csv(output_file, index=False, header=False)

        return output_file
    
    def _find_chromosome_mismatch(self, study_bim: Path, reference_bim: Path) -> Path:
        """
        Find chromosome mismatches between study and reference BIM files.

        This function identifies SNPs where the chromosome assignment differs between
        the study dataset and the reference panel, despite having the same rsID.
        Sex chromosomes (X, Y) are excluded from the update list.

        Parameters
        ----------
        study_bim : Path
            Path to the study BIM file to check for mismatches
        reference_bim : Path
            Path to the reference BIM file to compare against

        Returns
        -------
        Path
            Path to output file containing SNPs that need chromosome updates.
            File format is tab-separated with columns: chromosome, rsID
        """

        col_names = ["chr", "rsid", "pos_cm", "pos_bp", "allele1", "allele2"]
        study_df = pd.read_csv(study_bim, sep='\t', names=col_names)
        reference_df = pd.read_csv(reference_bim, sep='\t', names=col_names)

        # Find mismatches where rsID is the same but chromosome differs
        mismatch_df = reference_df.merge(study_df[["chr", "rsid"]], on="rsid", suffixes=("_ref", "_study"))
        chromosome_mismatch_df = mismatch_df[mismatch_df["chr_ref"] != mismatch_df["chr_study"]]

        # Exclude chromosomes X and Y from updates
        mismatch_df = mismatch_df[~mismatch_df["chr_study"].astype(str).isin(["X", "Y"])]

        to_update_chr_file = self.output_path / "all_phase3.toUpdateChr"

        # Save the mismatch data to a file
        chromosome_mismatch_df[["chr_study", "rsid"]].to_csv(to_update_chr_file, sep="\t", header=False, index=False)

        return to_update_chr_file
    
    def _find_position_mismatch(self, study_bim: Path, reference_bim: Path) -> Path:
        """
        Find SNPs with mismatched positions between study and reference datasets.

        This method compares the base pair positions of SNPs between a study dataset and a 
        reference dataset to identify SNPs that have different positions despite having the 
        same rsID.

        Parameters
        ----------
        study_bim : Path
            Path to the PLINK .bim file of the study dataset.
        reference_bim : Path
            Path to the PLINK .bim file of the reference dataset.

        Returns
        -------
        Path
            Path to the output file containing SNPs that need position updates.
            The output file contains two columns (rsID and new position) without headers.

        Notes
        -----
        The output file format is compatible with PLINK's --update-map command for updating
        SNP positions in the study dataset.
        """

        col_names = ["chr", "rsid", "pos_cm", "pos_bp", "allele1", "allele2"]
        study_df = pd.read_csv(study_bim, sep='\t', names=col_names)
        reference_df = pd.read_csv(reference_bim, sep='\t', names=col_names)

        # Create a dictionary from file1 with column 2 as key and column 4 as value
        a = dict(zip(study_df['rsid'], study_df['pos_bp']))

        # Filter rows in reference_df where column 2 exists in 'a' and the values isn column 4 differ
        filtered = reference_df[reference_df['rsid'].map(a).notna() & (reference_df['pos_bp'] != reference_df['rsid'].map(a))]

        # Print the result to a file
        to_update_pos_file = self.output_path / f"{self.reference_files['bim'].stem}.toUpdatePos"
        filtered[['rsid', 'pos_bp']].to_csv(to_update_pos_file, sep="\t", header=False, index=False)

        return to_update_pos_file
    
    def _find_allele_flip(self, study_bim: Path, reference_bim: Path, output_filename: Path) -> None:
        """
        Find SNPs with allele flips between study and reference datasets.

        This method identifies Single Nucleotide Polymorphisms (SNPs) where the alleles are
        flipped between the study and reference datasets. A flip occurs when the allele
        pairs don't match in either order.

        Parameters
        ----------
        study_bim : Path
            Path to the study .bim file containing SNP information
        reference_bim : Path
            Path to the reference .bim file containing SNP information
        output_filename : Path
            Path where the list of flipped SNPs will be saved

        Returns
        -------
        None
            Writes rsids of flipped SNPs to the specified output file

        Notes
        -----
        The .bim files should be tab-separated with columns:
        chromosome, rsid, genetic_distance, base_pair_position, allele1, allele2

        The output file will contain one rsid per line for SNPs where alleles don't match
        between study and reference in either order (A1/A2 or A2/A1).
        """

        col_names = ["chr", "rsid", "pos_cm", "pos_bp", "allele1", "allele2"]
        study_df = pd.read_csv(study_bim, sep='\t', names=col_names)
        reference_df = pd.read_csv(reference_bim, sep='\t', names=col_names)

        # Create a dictionary with the composite key from file1
        a = {f"{row['chr']}{row['rsid']}{row['pos_bp']}": f"{row['allele1']}{row['allele2']}" for _, row in study_df.iterrows()}

        # Filtering the rows in file2 based on the conditions
        filtered = reference_df[
            reference_df.apply(
                lambda row: (
                    f"{row['chr']}{row['rsid']}{row['pos_bp']}" in a and 
                    a[f"{row['chr']}{row['rsid']}{row['pos_bp']}"] not in {f"{row['allele1']}{row['allele2']}", f"{row['allele2']}{row['allele1']}"}
                ), axis=1
            )
        ]

        # Save the second column of filtered rows to a file
        filtered['rsid'].to_csv(output_filename, sep="\t", header=False, index=False)

        return
    
class GenomicOutlierAnalyzer:

    def __init__(self, input_path: Path, input_name: str, merged_file: Path, reference_tags: Path, output_path: Path, output_name: str) -> None:
        """
        Initialize GenomicOutlierAnalyzer object with input and output parameters.

        Parameters
        ----------
        input_path : Path
            Path to input directory containing files to process
        input_name : str
            Name of input file 
        merged_file : Path
            Path to merged genotype file
        reference_tags : Path
            Path to file containing reference population tags
        output_path : Path
            Path to output directory
        output_name : str
            Name for output files

        Attributes
        ----------
        einvectors : numpy.ndarray, None
            Principal component eigenvectors, initialized as None
        eigenvalues : numpy.ndarray, None
            Principal component eigenvalues, initialized as None
        ancestry_fails : list, None
            List of samples failing ancestry QC, initialized as None
        population_tags : pandas.DataFrame, None
            DataFrame containing population reference tags, initialized as None
        """

        self.merged_file = merged_file
        self.reference_tags = reference_tags
        self.output_path= output_path
        self.output_name= output_name
        self.input_path = input_path
        self.input_name = input_name

        self.einvectors = None
        self.eigenvalues = None
        self.ancestry_fails = None
        self.population_tags = None

        pass

    def execute_pca(self, pca: int = 10, maf: float = 0.01) -> None:
        """
        Perform Principal Component Analysis (PCA) on the genetic data using PLINK.

        This method executes PCA on the merged genetic data file, calculating the specified
        number of principal components. It automatically determines the optimal number of
        threads and memory allocation based on system resources.

        Parameters
        ----------
        pca : int, default=10
            Number of principal components to calculate.
            Must be a positive integer.
        maf : float, default=0.01
            Minor allele frequency threshold for filtering variants.
            Must be between 0 and 0.5.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If pca is not an integer or maf is not a float
        ValueError
            If pca is not positive or maf is not between 0 and 0.5

        Notes
        -----
        The method creates two output files:
        - {output_name}-pca.eigenvec: Contains the eigenvectors (PC loadings)
        - {output_name}-pca.eigenval: Contains the eigenvalues

        The results are stored in self.einvectors and self.eigenvalues attributes.
        """

        if not isinstance(pca, int):
            raise TypeError("pca should be an integer")
        if pca <= 0:
            raise ValueError("pca should be a positive integer")
        if not isinstance(maf, float):
            raise TypeError("maf should be a float")
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be a float between 0 and 0.5")

        logger.info("STEP: Performing principal component decomposition")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = cpu_count-2
        else:
            max_threads = 10

        # Get the virtual memory details
        memory_info = psutil.virtual_memory()
        available_memory_mb = memory_info.available / (1024 * 1024)
        memory = round(2*available_memory_mb/3,0)

        # PLINK command: generate PCA for reference data
        plink_cmd = f"plink --bfile {str(self.merged_file)} --keep-allele-order --maf {maf} --out {str(self.output_path / (self.output_name+'-pca'))} --pca {pca} --memory {memory} --threads {max_threads}"

        # execute PLINK command
        shell_do(plink_cmd, log=True)

        self.einvectors = self.output_path / (self.output_name+'-pca.eigenvec')
        self.eigenvalues = self.output_path / (self.output_name+'-pca.eigenval')

        return
    
    def find_ancestry_outliers(self, ref_threshold: float, stu_threshold: float, reference_pop: str, num_pcs: int = 2, fails_dir: Path = Path()) -> None:
        """
        Identifies ancestry outliers in the dataset based on PCA analysis.
        This method analyzes population structure using principal component analysis (PCA) and identifies
        samples that are potential ancestry outliers based on their distance from reference populations.
        
        Parameters
        ----------
        ref_threshold : float
            Distance threshold for reference population samples
        stu_threshold : float
            Distance threshold for study population samples
        reference_pop : str
            Name of the reference population to compare against
        num_pcs : int, optional
            Number of principal components to use in the analysis (default is 2)
        fails_dir : Path, optional
            Directory path to save failed samples information (default is empty Path)
        
        Returns
        -------
        None
            Results are stored in the ancestry_fails attribute

        Raises
        ------
        TypeError
            If parameters are not of the expected type
        ValueError
            If num_pcs is not a positive integer
        
        Notes
        -----
        The method requires:
        - A reference tags file with population information
        - An eigenvectors file from PCA analysis
        - Both files should be previously set in the class instance
        The results are saved in:
        - population_tags: CSV file with population assignments
        - ancestry_fails: List of samples identified as ancestry outliers
        """

        if not isinstance(ref_threshold, (float, int)):
            raise TypeError("ref_threshold should be a float")
        if not isinstance(stu_threshold, (float, int)):
            raise TypeError("stu_threshold should be a float")
        if not isinstance(reference_pop, str):
            raise TypeError("reference_pop should be a string")
        if not isinstance(num_pcs, int):
            raise TypeError("num_pcs should be an integer")
        if num_pcs <= 0:
            raise ValueError("num_pcs should be a positive integer")
        if not isinstance(fails_dir, Path):
            raise TypeError("fails_dir should be a Path object")
        
        if not fails_dir.exists():
            logger.info("STEP: Identifying ancestry outliers: `fails_dir` does not exist.")
            logger.info(f"STEP: Identifying ancestry outliers: ancestry outliers will be saved in {self.output_path}")
            fails_dir = self.output_path
        
        logger.info("STEP: Identifying ancestry outliers")

        df_tags = pd.read_csv(self.reference_tags, sep="\t", usecols=['#IID', 'SuperPop'])
        df_tags['ID'] = '0'
        df_tags = df_tags[['ID', '#IID', 'SuperPop']]
        df_tags = df_tags.rename(columns={'ID': 'ID1', '#IID': 'ID2', 'SuperPop': 'SuperPop'})

        if self.einvectors is None:
            raise ValueError("einvectors is not set. Make sure execute_pca() is called before this method and completed successfully.")

        df = pd.read_csv(self.einvectors, sep=r"\s+",engine='python', header=None)
        logger.info("STEP: Identifying ancestry outliers: read eigenvec file")

        df = df[[0, 1]]
        df = df.rename(columns = {0: 'ID1', 1:'ID2'})

        df['ID2'] = df['ID2'].astype(str)
        df['ID1'] = df['ID1'].astype(str)

        df = pd.merge(df, df_tags, on=['ID1', 'ID2'], how='left')
        df['SuperPop'] = df['SuperPop'].fillna('StPop', inplace=False)

        df.to_csv((self.output_path / (self.output_name + 'pop_tags.csv')), sep='\t', index=False)

        self.population_tags = self.output_path / (self.output_name + 'pop_tags.csv')

        # filter samples who are ethnicity outliers
        ancestry_fails = self._find_pca_fails(
            output_path  = fails_dir,
            df_tags      = df,
            ref_threshold= ref_threshold,
            stu_threshold= stu_threshold,
            reference_pop= reference_pop,
            num_pcs      = num_pcs
        )

        self.ancestry_fails = ancestry_fails

        return
    
    def execute_drop_ancestry_outliers(self, output_dir: Path = Path()) -> None:
        """
        Drop ancestry outliers from the study data by removing samples identified as ancestry outliers
        using PLINK command line tool.
        This method reads a file containing samples identified as ancestry outliers and creates new
        binary PLINK files excluding these samples.

        Parameters
        ----------
        output_dir : Path, optional
            Directory where the cleaned files will be saved. If not provided or doesn't exist,
            files will be saved in self.output_path.
        
        Returns
        -------
        None

        Raises
        ------
        TypeError
            If output_dir is not a Path object.
        
        Notes
        -----
        The method creates new PLINK binary files (.bed, .bim, .fam) with the suffix '-ancestry-cleaned'
        excluding the samples listed in self.ancestry_fails file.
        """

        logger.info("STEP: Dropping ancestry outliers from the study data")

        if not isinstance(output_dir, Path):
            raise TypeError("output_dir should be a Path object")
        
        if not output_dir.exists():
            logger.info("STEP: Dropping ancestry outliers from the study data: `output_dir` does not exist.")
            logger.info(f"STEP: Dropping ancestry outliers from the study data: ancestry outliers will be saved in {self.output_path}")
            output_dir = self.output_path

        if self.ancestry_fails is None:
            raise ValueError("ancestry_fails is not set. Make sure find_ancestry_outliers() is called before this method and completed successfully.")

        with open(self.ancestry_fails, 'r') as f:
            logger.info(f"STEP: Dropping ancestry outliers from the study data: {len(f.readlines())} samples identified as ancestry outliers")

        # create cleaned binary files
        plink_cmd2 = f"plink --bfile {str(self.input_path / self.input_name)} --allow-no-sex --remove {str(self.ancestry_fails)} --make-bed --out {str(output_dir / (self.output_name+'-ancestry-cleaned'))}"

        # execute PLINK command
        shell_do(plink_cmd2, log=True)

        return
    
    def draw_pca_plot(self, reference_pop: str, aspect_ratio: Union[Literal['auto', 'equal'], float], exclude_outliers: bool = False, plot_dir: Path = Path(), plot_name: str = 'pca_plot.pdf') -> None:
        """
        Generate 2D and 3D PCA plots from eigenvector data and population tags.
        This method creates two PCA visualization plots:
        - A 2D scatter plot showing PC1 vs PC2 colored by super-population
        - A 3D scatter plot showing PC1 vs PC2 vs PC3 colored by super-population
        
        Parameters
        ----------
        plot_dir : Path, optional
            Directory path where plots will be saved. Defaults to current directory.
            If directory doesn't exist, plots will be saved in self.output_path
        plot_name : str, optional
            Base name for the plot files. Defaults to 'pca_plot.jpeg'.
            Final filenames will be prefixed with '2D-' and '3D-'
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If plot_dir is not a Path object
            If plot_name is not a string
        
        Notes
        -----
        Requires the following class attributes to be set:
        - self.population_tags : Path to population tags file (tab-separated)
        - self.einvectors : Path to eigenvectors file (space-separated)
        - self.output_path : Path to output directory (used if plot_dir doesn't exist)
        The population tags file should contain columns 'ID1', 'ID2', and 'SuperPop'
        The eigenvectors file should contain the principal components data
        """

        logger.info("STEP: Generating PCA plots")

        if not isinstance(plot_dir, Path):
            raise TypeError("plot_dir should be a Path object")
        if not isinstance(plot_name, str):
            raise TypeError("plot_name should be a string")
        
        if not plot_dir.exists():
            logger.info('STEP: Generating PCA plots: `plot_dir` does not exist.')
            logger.info(f'STEP: Generating PCA plots: pca plots will be saved in {self.output_path}')
            plot_dir = self.output_path

        if self.population_tags is None:
            raise ValueError("population_tags is not set. Make sure find_ancestry_outliers() is called before this method and completed successfully.")

        # add population tags to pca output
        df_tags = pd.read_csv(self.population_tags, sep='\t')
        df_tags['ID1'] = df_tags['ID1'].astype(str)

        if self.einvectors is None:
            raise ValueError("einvectors is not set. Make sure execute_pca() is called before this method and completed successfully.")
        if self.eigenvalues is None:
            raise ValueError("eigenvalues is not set. Make sure execute_pca() is called before this method and completed successfully.")
        
        # load .eigenval file and calculate variance explained by the first two PCs
        df_eigenval = pd.read_csv(
            self.eigenvalues,
            header=None,
            sep   =r"\s+",
            engine='python'
        )

        total_variance = df_eigenval[0].sum()
        pc1_var = df_eigenval[0][0]
        pc2_var = df_eigenval[0][1]

        pc1_var_perc = round((pc1_var / total_variance) * 100, 2)
        pc2_var_perc = round((pc2_var / total_variance) * 100, 2)

        # load .eigenvec file and keep the first three principal components
        df_eigenvec = pd.read_csv(
            self.einvectors,
            header=None,
            sep   =r"\s+",
            engine='python'
        )
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:5]].copy()
        df_eigenvec.columns = ['ID1', 'ID2', 'pc_1', 'pc_2', 'pc_3']
        df_eigenvec['ID1'] = df_eigenvec['ID1'].astype(str)

        if exclude_outliers:
            # load ancestry outliers
            if self.ancestry_fails is None:
                raise ValueError("ancestry_fails is not set. Make sure find_ancestry_outliers() is called before this method and completed successfully.")
            logger.info("STEP: Generating PCA plots: excluding ancestry outliers")

            df_outliers = pd.read_csv(self.ancestry_fails, sep=r'\s+', header=None, engine='python')
            df_outliers.columns = ['ID1', 'ID2']
            df_outliers['ID1'] = df_outliers['ID1'].astype(str)
            df_outliers['ID2'] = df_outliers['ID2'].astype(str)

            df_eigenvec = df_eigenvec.merge(df_outliers, on=['ID1', 'ID2'], how='left', indicator=True)
            df_eigenvec = df_eigenvec[df_eigenvec['_merge'] == 'left_only'].drop(columns=['_merge'])

            plot_name = f'no-outliers-{plot_name}'

        # merge to get data with tagged populations
        df = pd.merge(df_eigenvec, df_tags, on=['ID1', 'ID2'])

        # generates a 2D scatter plot
        fig, ax = plt.subplots(figsize=(10,10))
        sns.scatterplot(data=df, x='pc_1', y='pc_2', hue='SuperPop', ax=ax, marker='.', s=70)
        ax.set_aspect(aspect_ratio, adjustable='datalim')
        plt.xlabel(f'PC_1 ({pc1_var_perc}%)')
        plt.ylabel(f'PC_2 ({pc2_var_perc}%)')
        fig.savefig(plot_dir / f'2D-aspect-{aspect_ratio}-{plot_name}', dpi=400)

        fig.clf()
        plt.close()

        fig3, ax3 = plt.subplots(figsize=(10,10))
        df_zoom = df[(df['SuperPop'] == 'StPop') | (df['SuperPop'] == reference_pop)].reset_index(drop=True)
        sns.scatterplot(data=df_zoom, x='pc_1', y='pc_2', hue='SuperPop', ax=ax3, marker='.', s=70)
        ax.set_aspect(aspect_ratio, adjustable='datalim')
        plt.xlabel(f'PC_1 ({pc1_var_perc}%)')
        plt.ylabel(f'PC_2 ({pc2_var_perc}%)')
        fig3.savefig(plot_dir / f'2D-zoom-aspect-{aspect_ratio}-{plot_name}', dpi=400)

        # generates a 3D scatter plot
        fig2= plt.figure()
        ax  = fig2.add_subplot(111, projection='3d')

        grouped = df.groupby('SuperPop')
        for s, group in grouped:
            ax.scatter(
                group['pc_1'],
                group['pc_2'],
                group['pc_3'],
                label=s
            )

        ax.legend()
        plt.savefig(plot_dir / f'3D-{plot_name}', dpi=400)
        plt.close()

        return
    
    def _set_population_tags(self, psam_path: Path, study_fam_path: Path) -> pd.DataFrame:
        """
        Sets population tags for genetic data by combining information from a PSAM file and a study FAM file.

        This method processes population information from reference data (PSAM file) and study data (FAM file), 
        combining them into a single DataFrame with consistent column naming and structure.

        Parameters
        ----------
        psam_path : Path
            Path to the PSAM file containing reference population information.
        study_fam_path : Path
            Path to the FAM file containing study individual IDs.

        Returns
        -------
        pd.DataFrame
            Combined DataFrame containing:
                - ID1: Family or group identifier (0 for reference data)
                - ID2: Individual identifier
                - SuperPop: Population tag ('StPop' for study individuals, actual population for reference data)

        Notes
        -----
        The PSAM file should contain at least '#IID' and 'SuperPop' columns.
        The FAM file should be space-separated with no header.
        """

        # Read population information from the .psam file
        df_psam = pd.read_csv(
            psam_path,
            sep='\t',
            usecols=['#IID', 'SuperPop']
        )

        # Set an ID column and rename columns for consistency
        df_psam['ID'] = 0
        df_psam = df_psam[['ID', '#IID', 'SuperPop']]
        df_psam.columns = ['ID1', 'ID2', 'SuperPop']

        # read individual IDs from the study .fam file
        df_fam = pd.read_csv(
            study_fam_path,
            sep=' ',
            header=None,
            index_col=False
        )

        # select relevant columns, assign a placeholder population tag, and rename columns
        df_fam = df_fam[df_fam.columns[:2]].copy()
        df_fam['SuperPop'] = 'StPop'
        df_fam.columns = ['ID1', 'ID2', 'SuperPop']

        # concatenate the two DataFrames to merge the information
        return pd.concat([df_fam, df_psam], axis=0)
    
    def _find_pca_fails(self, output_path: Path, df_tags: pd.DataFrame, ref_threshold: float, stu_threshold: float, reference_pop: str, num_pcs: int = 2) -> Path:
        """
        Identifies ancestry outliers based on PCA results using two thresholds:
        one for reference population and another for study population.

        Parameters
        ----------
        output_path : Path
            Path where the output file will be saved
        df_tags : pd.DataFrame
            DataFrame containing subject IDs and population tags
        ref_threshold : int
            Number of standard deviations from reference population mean to consider a subject as outlier
        stu_threshold : int
            Number of standard deviations from study population mean to consider a subject as outlier
        reference_pop : str
            Reference population name as it appears in df_tags
        num_pcs : int, optional
            Number of principal components to use in the analysis (default is 2)

        Returns
        -------
        str
            Path to the output file containing the IDs of subjects identified as ancestry outliers

        Raises
        ------
        TypeError
            If ref_threshold, stu_threshold are not numeric
            If reference_pop is not a string
            If num_pcs is not an integer
        ValueError
            If ref_threshold, stu_threshold are not positive
            If num_pcs is less than 1
            If num_pcs is greater than available PCs in eigenvec file

        Notes
        -----
        The method identifies outliers that deviate significantly from both:
        1. The reference population mean (by ref_threshold standard deviations)
        2. The study population mean (by stu_threshold standard deviations)
        Only subjects that are outliers in both criteria are included in the final output.
        """

        if not isinstance(ref_threshold, (float, int)):
            raise TypeError("ref_threshold should be an integer or float value")
        if not isinstance(stu_threshold, (float, int)):
            raise TypeError("stu_threshold should be an integer or float value")
        if stu_threshold<=0:
            raise ValueError("stu_threshold should be a positive value")
        if ref_threshold<=0:
            raise ValueError("ref_threshold should be a positive value")
        if not isinstance(reference_pop, str):
            raise TypeError("reference_pop should be a string")
        if not isinstance(num_pcs, int):
            raise TypeError("num_pcs should be an integer value")
        if num_pcs<1:
            raise ValueError("num_pcs should be a positive integer")

        # filters reference subjects
        mask1 = (df_tags['SuperPop']==reference_pop)
        # filters subjects from study data
        mask2 = (df_tags['SuperPop']=='StPop')

        # generates two data frames with filtered subjects
        df_ref = df_tags[mask1].reset_index(drop=True)
        df_stu = df_tags[mask2].reset_index(drop=True)

        if self.einvectors is None:
            raise ValueError("einvectors is not set. Make sure execute_pca() is called before this method and completed successfully.")

        # read .eigenvec file
        df_eigenvec = pd.read_csv(
            self.einvectors,
            header=None,
            sep   =r"\s+",
            engine='python'
        )

        if num_pcs>df_eigenvec.shape[1]-2:
            raise ValueError("num_pcs should be less than or equal to the number of principal components in the .eigenvec file")
        
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:2+num_pcs]].copy()

        # renames columns for consistency
        new_col_names = []
        for k in range(2+num_pcs):
            if k<2:
                new_col_names.append(f"ID{k+1}")
            else:
                new_col_names.append(f"pc_{k-1}")
        df_eigenvec.columns = new_col_names

        df_eigenvec['ID1'] = df_eigenvec['ID1'].astype(str)
        df_eigenvec['ID2'] = df_eigenvec['ID2'].astype(str)

        # merge filtered subjects with its principal components
        df_ref = df_ref.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)
        df_stu = df_stu.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)

        # computes mean and standard deviation by columns in reference data
        mean_ref= df_ref[df_ref.columns[2:]].mean()
        std_ref = df_ref[df_ref.columns[2:]].std()

        # creates empty data frame
        outliers_1 = pd.DataFrame(columns=df_ref.columns)
        outliers_1[df_stu.columns[:2]] = df_stu[df_stu.columns[:2]]

        # identifies subjects with more than `ref_threshold` std deviations from the reference mean
        for col in outliers_1.columns[2:]:
            outliers_1[col] = (np.abs(df_stu[col] - mean_ref[col]) > ref_threshold*std_ref[col])

        outliers_1['is_out'] = (np.sum(outliers_1.iloc[:,2:], axis=1) >0)

        df_1 = outliers_1[outliers_1['is_out']].reset_index(drop=True)[['ID1', 'ID2']].copy()

        # computes mean and standard deviation by columns in study data
        mean_stu= df_stu[df_stu.columns[2:]].mean()
        std_stu = df_stu[df_stu.columns[2:]].std()

        # creates empty data frame
        outliers_2 = pd.DataFrame(columns=df_ref.columns)
        outliers_2[df_stu.columns[:2]] = df_stu[df_stu.columns[:2]]

        # identifies subjects with more than `stu_threshold` std deviation from the study mean
        for col in outliers_2.columns[2:]:
            outliers_2[col] = (np.abs(df_stu[col] - mean_stu[col]) > stu_threshold*std_stu[col])

        outliers_2['is_out'] = (np.sum(outliers_2.iloc[:,2:], axis=1) >0)

        df_2 = outliers_2[outliers_2['is_out']].reset_index(drop=True)[['ID1', 'ID2']].copy()

        df = pd.merge(df_1, df_2, on=['ID1', 'ID2'])

        ancestry_fails = output_path / (self.output_name + '_fail-ancestry-qc.txt')

        logger.info(f"STEP: Identifying ancestry outliers: {df.shape[0]} samples identified as ancestry outliers")

        # save samples considered as ethnicity outliers
        df.to_csv(
            ancestry_fails,
            header=False,
            index =False,
            sep   ='\t'
        )

        return ancestry_fails

class AncestryQC:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name: str, high_ld_file: Path, reference_files: dict = dict(), recompute_merge: bool = True, built: str = '38', rename_snps: bool = False) -> None:
        """
        Initialize AncestryQC class.
        This class performs ancestry quality control analysis on genetic data by merging it with 1000 Genomes reference data
        and running principal component analysis.

        Parameters:
        -----------
        input_path: Path 
            Path to directory containing input files
        input_name: str 
            Base name of input files (without extension) 
        output_path: Path 
            Path to directory where output files will be saved
        output_name: str 
            Base name for output files
        high_ld_file: Path 
            Path to file containing high LD regions to exclude
        reference_files: dict (optional) 
            Dictionary with paths to reference files. Must contain 'bim', 'bed', 'fam' and 'psam' keys. 
            If not provided, will download 1000 Genomes reference files. Defaults to empty dict.
        recompute_merge: bool (optional): 
            Whether to recompute merge with reference even if merged files exist. Defaults to True.
        built: str (optional) 
            Genome build version, either '37' or '38'. Defaults to '38'.
        rename_snps: bool (optional): 
            Whether to rename SNPs to avoid duplicates during merge. Defaults to False.
        
        Raises:
        -------
            TypeError: If input arguments are not of expected types
            ValueError: If built is not '37' or '38'
            FileNotFoundError: If input_path or output_path do not exist
        
        Note:
        -----
            Creates the following directory structure under output_path:
            - ancestry_qc_results/
                - merging/
                - plots/ 
                - fail_samples/
                - clean_files/
        """

        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(high_ld_file, Path):
            raise TypeError("high_ld_regions should be a Path object")
        if not isinstance(reference_files, dict):
            raise TypeError("reference_files should be a dictionary")
        if not isinstance(input_name, str): 
            raise TypeError("input_name should be a string")
        if not isinstance(output_name, str):
            raise TypeError("output_name should be a string")
        if not isinstance(recompute_merge, bool):
            raise TypeError("recompute_merge should be a boolean")
        if not isinstance(built, str):
            raise TypeError("built should be a string")
        if built not in ['37', '38']:
            raise ValueError("built should be either '37' or '38'")
        if not isinstance(rename_snps, bool):
            raise TypeError("rename_snps should be a boolean")
        
        if not input_path.exists():
            raise FileNotFoundError("input_path does not exist")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")
        if not high_ld_file.is_file():
            logger.info(f"High LD file not found at {high_ld_file}")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions(built=built)
            ld_fetcher.get_ld_regions()

            ld_regions = ld_fetcher.ld_regions
            if ld_regions is None:
                raise ValueError("Failed to fetch high LD regions file")
            high_ld_file = ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_file}")
        
        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.output_name= output_name
        self.reference_files = reference_files
        self.high_ld_regions = high_ld_file
        self.recompute_merge = recompute_merge
        self.built = built
        self.rename_snps = rename_snps

        if not reference_files:

            logger.info(f"No reference files provided. Fetching 1000 Genomes reference data for built {self.built}")

            fetcher = Fetcher1000Genome(built=self.built)
            fetcher.get_1000genomes()
            fetcher.get_1000genomes_binaries()

            self.reference_files = {
                'bim': fetcher.bim_file,
                'bed': fetcher.bed_file,
                'fam': fetcher.fam_file,
                'psam': fetcher.psam_file
            }

        self.results_dir = self.output_path / 'ancestry_qc_results' 
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.merging_dir = self.results_dir / 'merging'
        self.merging_dir.mkdir(parents=True, exist_ok=True)

        self.plots_dir = self.results_dir / 'plots'
        self.plots_dir.mkdir(parents=True, exist_ok=True)

        self.fail_samples_dir = self.results_dir / 'fail_samples'
        self.fail_samples_dir.mkdir(parents=True, exist_ok=True)

        self.clean_files = self.results_dir / 'clean_files'
        self.clean_files.mkdir(parents=True, exist_ok=True)

        pass

    def merge_reference_study(self, ind_pair: list = [50, 5, 0.2]) -> None:
        """
        Merge reference and study data by applying quality control filters and merging steps.
        This method performs a series of quality control steps to merge study data with reference data:
        1. Filters problematic SNPs
        2. Performs LD pruning
        3. Fixes chromosome mismatches
        4. Fixes position mismatches  
        5. Fixes allele flips
        6. Removes remaining mismatches
        7. Merges the datasets
        
        Parameters
        ----------
        ind_pair : list, default [50, 5, 0.2]
            Parameters for LD pruning: [window size, step size, r2 threshold]
        
        Returns
        -------
        None
        
        Notes
        -----
        If recompute_merge is False, the method will skip the merging process and expect
        merged data to already exist in the merging directory.
        
        Raises
        ------
        TypeError
            If ind_pair is not a list
        """

        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        if not self.recompute_merge:
            logger.info("STEP: Merging study and reference data: recompute_merge is set to False. Skipping merging step")
            logger.info(f"STEP: Merging study and reference data: merged data is expected to be in {self.merging_dir}")
            return

        rgm = ReferenceGenomicMerger(
            input_path= self.input_path,
            input_name= self.input_name,
            output_path= self.merging_dir, 
            output_name= self.output_name,
            high_ld_regions =self.high_ld_regions, 
            reference_files = self.reference_files,
        )

        rgm.execute_rename_snpid()
        rgm.execute_filter_prob_snps()
        rgm.execute_ld_pruning(ind_pair=ind_pair)
        rgm.execute_fix_chromosome_mismatch()
        rgm.execute_fix_possition_mismatch()
        rgm.execute_fix_allele_flip()
        rgm.execute_remove_mismatches()
        rgm.execute_merge_data()

        return
    
    def _clean_merging_dir(self) -> None:
        """
        Cleans up the merging directory by removing unnecessary files.
        This method removes all files in the merging directory except:
        - Files containing '-merged' in their name
        - Log files with '.log' extension
        The cleanup helps manage disk space and removes intermediate files that are no longer needed
        after the merging process is complete.
        
        Returns:
        --------

            None
        """
        

        for file in self.merging_dir.iterdir():
            if file.is_file() and '-merged' not in file.name and file.suffix != '.log':
                file.unlink()

        return
    
    def run_pca(self, ref_population: str, pca: int = 10, maf: float = 0.01, num_pca: int = 10, ref_threshold: float = 4, stu_threshold: float = 4, aspect_ratio: Union[Literal['auto', 'equal'], float]='equal') -> None:
        """
        Performs Principal Component Analysis (PCA) on genetic data and identifies ancestry outliers.

        This method executes a complete PCA workflow including:
        1. Running the PCA analysis
        2. Identifying ancestry outliers
        3. Removing identified outliers
        4. Generating PCA plots

        Parameters
        ----------
        ref_population : str
            Reference population identifier for ancestry comparison
        pca : int, optional
            Number of principal components to calculate (default=10)
        maf : float, optional
            Minor allele frequency threshold for filtering (default=0.01)
        num_pca : int, optional
            Number of principal components to use in outlier detection (default=10)
        ref_threshold : float, optional
            Threshold for reference population outlier detection (default=4)
        stu_threshold : float, optional
            Threshold for study population outlier detection (default=4)

        Returns
        -------
        None
            Results are saved to specified output directories

        Notes
        -----
        The method uses the GenomicOutlierAnalyzer class to perform the analysis and 
        saves results in the directories specified during class initialization.
        """

        # Make sure the reference tag path is valid before creating the analyzer
        if 'psam' not in self.reference_files or not isinstance(self.reference_files['psam'], Path):
            raise ValueError("Reference files dictionary must contain a valid 'psam' Path")
        
        goa = GenomicOutlierAnalyzer(
            input_path= self.input_path, 
            input_name= self.input_name,
            merged_file= self.merging_dir / (self.output_name + '-merged'),
            reference_tags= self.reference_files['psam'],
            output_path= self.results_dir, 
            output_name= self.output_name
        )

        logger.info(f"STEP: Running PCA analysis: `ref_population` = {ref_population}")
        logger.info(f"STEP: Running PCA analysis: `pca` = {pca}")
        logger.info(f"STEP: Running PCA analysis: `maf` = {maf}")
        logger.info(f"STEP: Running PCA analysis: `num_pca` = {num_pca}")
        logger.info(f"STEP: Running PCA analysis: `ref_threshold` = {ref_threshold}")
        logger.info(f"STEP: Running PCA analysis: `stu_threshold` = {stu_threshold}")
        logger.info(f"STEP: Running PCA analysis: `psam_file` = {self.reference_files['psam']}")

        goa.execute_pca(pca=pca, maf=maf)
        goa.find_ancestry_outliers(
            ref_threshold=ref_threshold, 
            stu_threshold=stu_threshold, 
            reference_pop=ref_population, 
            num_pcs      =num_pca, 
            fails_dir    =self.fail_samples_dir
        )
        goa.execute_drop_ancestry_outliers(output_dir=self.clean_files)
        goa.draw_pca_plot(plot_dir=self.plots_dir, reference_pop=ref_population, aspect_ratio=aspect_ratio)
        goa.draw_pca_plot(plot_dir=self.plots_dir, reference_pop=ref_population, aspect_ratio=aspect_ratio, exclude_outliers=True)

        return
