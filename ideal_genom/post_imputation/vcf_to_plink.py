
import logging
import os
import subprocess
import psutil

from pathlib import Path
from typing import Optional
from typing import Optional

from core.utils import validate_input_file
from core.executor import run_plink2

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

        if threads is None:
            cpu_count = os.cpu_count()
            if cpu_count is not None:
                threads = max(1, cpu_count - 2)  # Ensure at least 1 thread
            else:
                threads = 10

        if memory is None:
            # get virtual memory details
            memory_info = psutil.virtual_memory()
            available_memory_mb = memory_info.available / (1024 * 1024)
            memory = int(round(2*available_memory_mb/3,0))

        if double_id:
            # plink2 command
            plink2_cmd = [
                "plink2",
                "--vcf", (self.input_path / self.input_name).as_posix(),
                "--snps-only", "just-acgt", "--double-id",
                "--make-bed",
                "--out", (self.analysis_ready / (self.output_name + "-nosex")).as_posix(),
                "--threads", str(threads),
                "--memory", str(memory)
            ]
        else:
            # plink2 command
            plink2_cmd = [
                "plink2",
                "--vcf", (self.input_path / self.input_name).as_posix(),
                "--snps-only", "just-acgt",
                "--make-bed",
                "--out", (self.analysis_ready / (self.output_name + "-nosex")).as_posix(),
                "--threads", str(threads),
                "--memory", str(memory)
            ]

        # execute plink2 command
        try:
            run_plink2(plink2_cmd)
            print(f"PLINK2 command executed successfully. Output files saved with prefix: {self.output_name}-nosex")
        except subprocess.CalledProcessError as e:
            print(f"Error running PLINK2: {e}")

        pass

    def update_fam(self, for_fam_update_file: Path, threads: Optional[int] = None) -> None:
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

        if not isinstance(for_fam_update_file, Path):
            raise TypeError(f"for_fam_update_file should be of type Path, got {type(for_fam_update_file)}")

        fam_file = for_fam_update_file

        # Ensure the path points to a .fam file
        if not fam_file.as_posix().endswith('.fam'):
            fam_file = Path(str(fam_file) + '.fam')

        # Validate the fam file exists and is a file
        validate_input_file(fam_file, ['.fam'])

        logger.info(f"Updating family information in {self.output_name}-nosex.fam with {for_fam_update_file}")
   
        if threads is None:
            cpu_count = os.cpu_count()
            if cpu_count is not None:
                threads = max(1, cpu_count - 2)  # Ensure at least 1 thread
            else:
                threads = 10
        
        # PLINK2 command
        run_plink2([
            "plink2",
            "--bfile", (self.analysis_ready / (self.output_name + "-nosex")).as_posix(),
            "--make-bed",
            "--out", (self.analysis_ready / self.output_name).as_posix(),
            "--fam", fam_file.as_posix(),
            "--threads", str(threads)
        ])

        pass
    