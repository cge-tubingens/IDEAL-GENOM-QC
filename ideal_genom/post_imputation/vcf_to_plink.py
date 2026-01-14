
"""Module for converting VCF files to PLINK binary format.

This module provides the GetPLINK class for converting post-imputation VCF files
to PLINK binary format for downstream genomic analysis. The conversion process
handles resource management, file validation, and supports both initial conversion
and family information updates.
"""

import logging

from pathlib import Path
from typing import Optional

import pandas as pd

from ..core.utils import validate_input_file, get_optimal_threads, get_available_memory
from ..core.executor import run_plink2

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class GetPLINK:

    """A class for converting concatenated VCF files to PLINK binary format.

    This class handles conversion of a concatenated VCF file to a PLINK binary ready for further analysis.

    Attributes
    ----------
    input_path : Path
        Directory path where the input VCF file is located.
    output_path : Path
        Directory path where the output files will be saved.
    input_name : str
        Name of the input VCF file (must end with .vcf or .vcf.gz).
    output_name : str, optional
        Name for the output file. If not provided, it will be derived from input_name.
        
    Raises
    ------
    TypeError
        If input_path or output_path is not a Path object, or if input_name or output_name is not a string.
    FileNotFoundError
        If input_path or output_path does not exist.
    NotADirectoryError
        If input_path or output_path is not a directory.
    ValueError
        If input_name is not provided or if it doesn't end with .vcf or .vcf.gz.
    """

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name: str) -> None:
        

        if not isinstance(input_path, Path):
            raise TypeError(f"input_path should be of type Path, got {type(input_path)}")
        if not isinstance(output_path, Path):
            raise TypeError(f"output_path should be of type Path, got {type(output_path)}")
        if not input_path.exists():
            raise FileNotFoundError(f"Input path {input_path} does not exist.")
        if not output_path.exists():
            raise FileNotFoundError(f"Output path {output_path} does not exist.")
        if not input_path.is_dir():
            raise NotADirectoryError(f"Output path {input_path} is not a directory.")
        if not output_path.is_dir():
            raise NotADirectoryError(f"Output path {output_path} is not a directory.")
        
        if input_name:
            if not isinstance(input_name, str):
                raise TypeError(f"input_name should be of type str, got {type(input_name)}")
            # Validate the actual file exists and has correct extension
            input_file_path = input_path / input_name
            validate_input_file(input_file_path, ['.vcf', '.vcf.gz'])
        else:
            raise ValueError("input_name must be provided")

        if not output_name:
            output_name = input_name.split('.vcf')[0]
        else:
            if not isinstance(output_name, str):
                raise TypeError(f"output_name should be of type str, got {type(output_name)}")
        
        self.input_path = input_path
        self.output_path= output_path
        self.input_name = input_name
        self.output_name = output_name

        self.analysis_ready = self.input_path / 'analysis_ready'
        self.analysis_ready.mkdir(parents=True, exist_ok=True)

        self.updated_fam = None
        self.renamed_snps = None
        
        pass

    def convert_vcf_to_plink(self, double_id: bool = True, threads: Optional[int] = None, memory: Optional[int] = None) -> None:
        """Convert a VCF file to PLINK binary format (.bed, .bim, .fam).
        
        This method runs the plink2 command-line tool to convert the input VCF file to PLINK
        binary format, filtering for SNPs with standard ACGT alleles only.
        
        Parameters
        ----------
        double_id : bool, optional 
            Whether to use the --double-id flag in plink2 command, which sets both FID and IID to the sample ID. Defaults to True.
        threads : int, optional 
            Number of CPU threads to use. If None, defaults to (available CPU cores - 2) or 10 if CPU count can't be determined.
        memory : int, optional
            Memory allocation in MB for plink2. If None, defaults to approximately 2/3 of available system memory.

        Returns
        -------
        None
        
        Side Effects
        ------------
        Creates PLINK binary files (.bed, .bim, .fam) in the self.analysis_ready directory with the prefix self.output_name + "-nosex".
        
        Raises
        ------
        subprocess.CalledProcessError 
            If the plink2 command execution fails.
        """

        threads = threads or get_optimal_threads()
        memory = memory or get_available_memory()

        if double_id:
            # plink2 command
            plink2_args = [
                "--vcf", (self.input_path / self.input_name).as_posix(),
                "--snps-only", "just-acgt", "--double-id",
                "--make-bed",
                "--out", (self.analysis_ready / (self.output_name + "-nosex")).as_posix(),
                "--threads", str(threads),
                "--memory", str(memory)
            ]
        else:
            # plink2 command
            plink2_args = [
                "--vcf", (self.input_path / self.input_name).as_posix(),
                "--snps-only", "just-acgt",
                "--make-bed",
                "--out", (self.analysis_ready / (self.output_name + "-nosex")).as_posix(),
                "--threads", str(threads),
                "--memory", str(memory)
            ]

        # execute plink2 command
        try:
            run_plink2(plink2_args)
            print(f"PLINK2 command executed successfully. Output files saved with prefix: {self.output_name}-nosex")
        except Exception as e:
            print(f"Error running PLINK2: {e}")

        pass

    def update_fam(self, for_fam_update_file: Path, threads: Optional[int] = None, memory: Optional[int] = None) -> None:
        """Add family information to the PLINK .fam file.

        This method reads a family information file and updates the PLINK .fam file
        using the provided family information, via PLINK2.

        Parameters
        ----------
        for_fam_update_file : Path
            Path to the family information file (.fam or without suffix).
        threads : int, optional
            Number of threads to use for PLINK2 (defaults to available CPUs - 2).

        Returns
        -------
        None
        """
        if not for_fam_update_file:
            raise ValueError("for_fam_update_file must be provided")
        
        if isinstance(for_fam_update_file, str):
            for_fam_update_file = Path(for_fam_update_file)

        if not isinstance(for_fam_update_file, Path):
            raise TypeError(f"for_fam_update_file should be of type Path, got {type(for_fam_update_file)}")

        fam_file = for_fam_update_file

        # Ensure the path points to a .fam file
        if not fam_file.as_posix().endswith('.fam'):
            fam_file = Path(str(fam_file) + '.fam')

        # Validate the fam file exists and is a file
        validate_input_file(fam_file, ['.fam'])

        logger.info(f"Updating family information in {self.output_name}-nosex.fam with {for_fam_update_file}")
   
        threads = threads or get_optimal_threads()
        memory = memory or get_available_memory()
        
        # PLINK2 command
        plink2_args = [
            "--bfile", (self.analysis_ready / (self.output_name + "-nosex")).as_posix(),
            "--make-bed",
            "--out", (self.analysis_ready / (self.output_name + '-updated')).as_posix(),
            "--fam", fam_file.as_posix(),
            "--threads", str(threads),
            "--memory", str(memory)
        ]
        try:
            run_plink2(plink2_args)
            logger.info(f"Family information updated successfully. Output files saved with prefix: {self.output_name}")
        except Exception as e:
            logger.error(f"Error updating family information with PLINK2: {e}")

        self.updated_fam = self.analysis_ready / (self.output_name + '-updated')

        pass

    def execute_intermediate_cleanup(self) -> None:
        """Cleans up intermediate files generated during the PLINK conversion process.

        This method removes temporary files such as the initial .bed, .bim, and .fam files
        created during the conversion from VCF to PLINK binary format.

        Returns
        -------
        None
        """
        import os

        intermediate_files = [
            self.analysis_ready / (self.output_name + "-nosex.bed"),
            self.analysis_ready / (self.output_name + "-nosex.bim"),
            self.analysis_ready / (self.output_name + "-nosex.fam"),
            self.analysis_ready / (self.output_name + "-updated.bed"),
            self.analysis_ready / (self.output_name + "-updated.bim"),
            self.analysis_ready / (self.output_name + "-updated.fam"),
        ]

        for file_path in intermediate_files:
            try:
                if file_path.exists():
                    os.remove(file_path)
                    logger.info(f"Removed intermediate file: {file_path}")
            except Exception as e:
                logger.error(f"Error removing intermediate file {file_path}: {e}")

        return
    
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

        max_threads = get_optimal_threads()
        memory = get_available_memory()

        plink2_args = [
            '--bfile', str(self.updated_fam),
            '--set-all-var-ids', '@:#:$r:$a',
            '--threads', str(max_threads),
            '--memory', str(memory),
            '--make-bed',
            '--out', str(self.analysis_ready / (self.output_name ))
        ]

        self._get_id_links()

        # Execute PLINK2 command
        run_plink2(plink2_args)

        self.renamed_snps = str(self.input_path / (self.output_name ))

        return


    def execute_plink_conversion_pipeline(self, plink_params: dict) -> None:
        """Execute the full PLINK conversion pipeline: VCF to PLINK binary and optional family info update.

        This method orchestrates the conversion of a VCF file to PLINK binary format
        and optionally updates the .fam file with family information.

        Parameters
        ----------
        double_id : bool, optional
            Whether to use the --double-id flag in plink2 command. Defaults to True.
        for_fam_update_file : Path, optional
            Path to the family information file for updating the .fam file. If None, no update is performed.
        threads : int, optional
            Number of CPU threads to use. If None, defaults to (available CPU cores - 2) or 10 if CPU count can't be determined.
        memory : int, optional
            Memory allocation in MB for plink2. If None, defaults to approximately 2/3 of available system memory.

        Returns
        -------
        None
        """

        double_id           = plink_params.get('double_id', True)
        for_fam_update_file = plink_params.get('for_fam_update_file', None)
        threads             = plink_params.get('threads', None)
        memory              = plink_params.get('memory', None)

        self.convert_vcf_to_plink(double_id=double_id, threads=threads, memory=memory)

        if for_fam_update_file:
            self.update_fam(for_fam_update_file=for_fam_update_file, threads=threads, memory=memory)
        else:
            logger.info("No family update file provided; skipping .fam update step.")
            logger.info("Retaining initial PLINK files with '-nosex' suffix.")

        self.execute_rename_snpid()
        self.execute_intermediate_cleanup()

        return
    
    def _get_id_links(self)->None:

        """Generates a mapping of original SNP IDs to renamed SNP IDs.

        This private method reads the BIM file of the input PLINK dataset and
        creates a mapping of original SNP IDs to a new format based on chromosome,
        position, and alleles. The new SNP ID format is 'chr:pos:a1:a2'.

        Returns
        -------
        None

        Notes
        -----
        The mapping is stored in the `snp_id_map` attribute as a dictionary.
        """
        if self.updated_fam is None:
            raise ValueError("Updated fam file not found. Please run update_fam() before generating SNP ID links.")
        
        bim_file = self.updated_fam.with_suffix('.bim')
        chunk_size = 1000000  # Number of lines to read at a time

        df_links = pd.DataFrame()

        for chunk in pd.read_csv(bim_file, sep='\t', header=None, chunksize=chunk_size):
            # BIM columns: chrom, snp, cm, pos, a1, a2
            chunk['new_snp_id'] = chunk[0].astype(str) + ':' + chunk[3].astype(str) + ':' + chunk[4].astype(str) + ':' + chunk[5].astype(str)
            df_links = pd.concat([df_links, chunk[[1, 'new_snp_id']]], ignore_index=True)

        df_links.columns = ['original_snp_id', 'renamed_snp_id']
        df_links.to_csv(self.analysis_ready / f"{self.output_name}_snp_id_links.tsv", sep='\t', index=False)
    