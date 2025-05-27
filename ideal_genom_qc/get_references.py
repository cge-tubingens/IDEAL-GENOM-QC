import requests
import logging

import pandas as pd

from pathlib import Path
from typing import Optional

from ideal_genom_qc.Helpers import shell_do

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class Fetcher1000Genome:

    def __init__(self, destination: Optional[Path] = None, built: str = '38'):
        """
        Initialize a reference data handler.
        This class manages reference data files from 1000 Genomes Project.
        
        Paramters:
        ----------
        destination: Path (optional) 
            Path where reference files will be stored. If not provided, defaults to '../data/1000genomes_built_{built}'.
        built: str (optional): 
            Human genome build version. Defaults to '38'.
        
        Attributes:
        -----------
        destination: Path 
            Directory path where reference files are stored
        built: str 
            Human genome build version being used
        pgen_file: Path 
            Path to PGEN format file
        pvar_file: Path 
            Path to PVAR format file 
        psam_file: Path 
            Path to PSAM format file
        bed_file: Path 
            Path to BED format file
        bim_file: Path 
            Path to BIM format file 
        fam_file: Path 
            Path to FAM format file
        """

        if not destination:
            destination = Path(__file__).resolve().parent.parent / "data" / f"1000genomes_built_{built}"

        logger.info(f"Destination folder: {destination}")
        
        self.destination = destination
        self.built = built

        self.pgen_file = None
        self.pvar_file = None
        self.psam_file = None

        self.bed_file = None
        self.bim_file = None
        self.fam_file = None
        
        pass

    def get_1000genomes(self, url_pgen: Optional[str] = None, url_pvar: Optional[str] = None, url_psam: Optional[str] = None)-> Path:
        """
        Download and decompress 1000 Genomes reference data.
        This method downloads the PLINK2 binary files (.pgen, .pvar, .psam) for the 1000 Genomes 
        reference dataset, corresponding to the specified genome build (37 or 38). If the files 
        already exist in the destination directory, the download is skipped.
        
        Parameters:
        -----------
        url_pgen (str, optional): Custom URL for downloading the .pgen file. 
            If None, uses default URL based on genome build.
        url_pvar (str, optional): Custom URL for downloading the .pvar file.
            If None, uses default URL based on genome build.
        url_psam (str, optional): Custom URL for downloading the .psam file.
            If None, uses default URL based on genome build.
        
        Returns:
        --------
            Path: Path object pointing to the decompressed .pgen file location.
        
        Note:
        -----
            The method requires plink2 to be installed and accessible in the system path
            for decompressing the .pgen file.
        """

        self.destination.mkdir(parents=True, exist_ok=True)

        if self.built == '38':
            if url_pgen is None:
                url_pgen = r"https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1"
            if url_pvar is None:
                url_pvar = r"https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1"
            if url_psam is None:
                url_psam = r"https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1"
        
        elif self.built == '37':
            if url_pgen is None:
                url_pgen = r"https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1"
            if url_pvar is None:
                url_pvar = r"https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1"
            if url_psam is None:
                url_psam = r"https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1"

        if self._check_if_binaries_exist():

            logger.info("1000 Genomes binaries already exist. Skipping download.")

            self.pvar_file = self.destination / "all_phase3.pvar.zst"
            self.psam_file = self.destination / "all_phase3.psam"
            self.pgen_decompressed = self.destination / "all_phase3.pgen"

            return self.pgen_decompressed
        
        logger.info("Downloading 1000 Genomes data...")
        
        if url_pgen is not None:
            self.pgen_file = self._download_file(url_pgen, self.destination / "all_phase3.pgen.zst")
        if url_pvar is not None:
            self.pvar_file = self._download_file(url_pvar, self.destination / "all_phase3.pvar.zst")
        if url_psam is not None:
            self.psam_file = self._download_file(url_psam, self.destination / "all_phase3.psam")

        pgen_file = self.destination / "all_phase3.pgen.zst"
        pgen_decompressed = self.destination / "all_phase3.pgen"

        logger.info("Decompressing pgen file from 1000 Genomes data...")

        # plink2 command
        plink2_cmd = f"plink2 --zst-decompress {pgen_file} {pgen_decompressed}"

        # execute plink2 command
        shell_do(plink2_cmd)

        self.pgen_file = pgen_decompressed

        return pgen_decompressed

    def _download_file(self, url: str, destination: Path) -> Path:
        """
        Downloads a file from a given URL and saves it to the specified destination.

        Parameters:
        -----------
        url: str 
            The URL of the file to download.
        destination: Path 
            The path where the downloaded file will be saved.

        Returns:
        --------
        Path: 
            The path to the downloaded file if successful, None if download fails.

        Raises:
        -------
            None: Exceptions are caught and logged internally.
        """
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            with open(destination, 'wb') as out_file:
                for chunk in response.iter_content(chunk_size=1024):
                    out_file.write(chunk)
            return destination
        except requests.RequestException as e:
            logging.error(f"Failed to download {url}: {e}")
            return Path()

    def get_1000genomes_binaries(self) -> Path:
        """
        Convert downloaded 1000 Genomes data into PLINK binary files (.bed, .bim, .fam).
        This method processes the downloaded 1000 Genomes data files and converts them into PLINK binary format.
        If the binary files already exist, it skips the conversion process. The method handles file cleanup
        and proper renaming of output files.
        The conversion is done in two steps:
        1. Convert pfile to binary format including only SNPs from chromosomes 1-22,X,Y,MT
        2. Update variant IDs and create final binary files
        
        Returns
        -------
        Path
            Path object pointing to the generated binary files (without extension)
            The actual files created will be .bed, .bim, .fam and .psam with the same prefix
        """

        if self._check_if_binaries_exist():

            logger.info("1000 Genomes binaries already exist. Skipping conversion into bfiles...")

            (self.destination / "all_phase3.pgen").unlink(missing_ok=True)
            (self.destination / "all_phase3.pgen.zst").unlink(missing_ok=True)
            (self.destination / "all_phase3.pvar.zst").unlink(missing_ok=True)

            self.bed_file = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.bed')
            self.bim_file = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.bim')
            self.fam_file = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.fam')
            self.psam_file = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.psam')

            return self.destination / "all_phase3"
        
        logger.info("Converting 1000 Genomes data into bfiles...")

        # plink2 command
        plink2_cmd = f"plink2 --pfile {self.destination / 'all_phase3'} vzs --chr 1-22,X,Y,MT --snps-only --max-alleles 2 --make-bed --out {self.destination / 'all_phase3'}"
        
        # execute plink2 command
        shell_do(plink2_cmd)

        (self.destination / "all_phase3.pgen").unlink(missing_ok=True)
        (self.destination / "all_phase3.pgen.zst").unlink(missing_ok=True)
        (self.destination / "all_phase3.pvar.zst").unlink(missing_ok=True)

        logger.info("Downloaded 1000 Genomes data deleted.")

        # PLINK2 command
        plink2_cmd = f"plink2 --bfile {self.destination / 'all_phase3'} --set-all-var-ids @:#:$r:$a --make-bed --out {self.destination / f'1kG_phase3_GRCh{self.built}'}"

        shell_do(plink2_cmd, log=True)

        (self.destination / "all_phase3.bed").unlink(missing_ok=True)
        (self.destination / "all_phase3.bim").unlink(missing_ok=True)
        (self.destination / "all_phase3.fam").unlink(missing_ok=True)

        self.bed_file = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.bed')
        self.bim_file = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.bim')
        self.fam_file = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.fam')

        psam_renamed = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.psam')
        self.psam_file= (self.destination / "all_phase3.psam").rename(psam_renamed)

        return self.destination / f'1kG_phase3_GRCh{self.built}'
    
    def _check_if_binaries_exist(self) -> bool:
        """
        Checks if all required binary files exist in the destination directory.

        This method verifies the existence of .bed, .bim, .fam, and .psam files
        for the 1000 Genomes Phase 3 reference panel in the specified genome build.

        Returns:
        --------
            bool: True if all required files exist, False otherwise.
        """

        check_bed = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.bed').exists()
        check_bim = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.bim').exists()
        check_fam = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.fam').exists()
        check_psam = (self.destination / f'1kG_phase3_GRCh{self.built}').with_suffix('.psam').exists()

        return check_bed and check_bim and check_fam and check_psam
    
class FetcherLDRegions:

    def __init__(self, destination: Optional[Path] = None, built: str = '38'):
        """
        Initialize LDRegions object.
        This initializer sets up the destination path for LD regions files and the genome build version.
        If no destination is provided, it defaults to a 'data/ld_regions_files' directory relative to
        the parent directory of the current file.
        
        Parameters
        ----------
        destination : Path, optional
            Path where LD region files will be stored. If None, uses default path.
        built : str, optional
            Genome build version, defaults to '38'.
        
        Attributes
        ----------
        destination : Path
            Directory path where LD region files are stored
        built : str
            Genome build version being used
        ld_regions : None
            Placeholder for LD regions data, initially set to None
        """

        if not destination:
            destination = Path(__file__).resolve().parent.parent / "data" / "ld_regions_files"

        self.destination = destination
        self.built = built

        self.ld_regions = None
        
        pass

    def get_ld_regions(self)-> Path:
        """
        Downloads or creates high LD regions file based on genome build version.
        This method handles the retrieval of high Linkage Disequilibrium (LD) regions for
        different genome builds (37 or 38). For build 37, it downloads the regions from a
        GitHub repository. For build 38, it creates the file from predefined coordinates.
        
        Returns:
        --------
        Path: Path to the created/downloaded LD regions file. Returns empty Path if
              download fails for build 37.
        
        Raises:
        -------
            None explicitly, but may raise standard I/O related exceptions.
        
        Notes:
        -----
            - For build 37: Downloads from genepi-freiburg/gwas repository
            - For build 38: Creates file from hardcoded coordinates from GWAS-pipeline
            - Files are named as 'high-LD-regions_GRCh{build}.txt'
            - Creates destination directory if it doesn't exist
        """

        self.destination.mkdir(parents=True, exist_ok=True)

        out_dir = self.destination

        if self.built == '37':
            url_ld_regions = r"https://raw.githubusercontent.com/genepi-freiburg/gwas/refs/heads/master/single-pca/high-LD-regions.txt"
        
            ld = requests.get(url_ld_regions)


            if ld.status_code == 200:
                with open((out_dir / f"high-LD-regions_GRCh{self.built}.txt"), "wb") as f:
                    f.write(ld.content)
                logger.info(f"LD regions file for built {self.built} downloaded successfully to {out_dir}")

                self.ld_regions = out_dir / f"high-LD-regions_GRCh{self.built}.txt"
                return out_dir / f"high-LD-regions_GRCh{self.built}.txt"
            else:
                logger.info(f"Failed to download .bim file: {ld.status_code}")

                return Path()

        elif self.built == '38':
            # extracted from
            # https://github.com/neurogenetics/GWAS-pipeline
            data = [
                (1, 47534328, 51534328, "r1"),
                (2, 133742429, 137242430, "r2"),
                (2, 182135273, 189135274, "r3"),
                (3, 47458510, 49962567, "r4"),
                (3, 83450849, 86950850, "r5"),
                (5, 98664296, 101164296, "r6"),
                (5, 129664307, 132664308, "r7"),
                (5, 136164311, 139164311, "r8"),
                (6, 24999772, 35032223, "r9"),
                (6, 139678863, 142178863, "r10"),
                (8, 7142478, 13142491, "r11"),
                (8, 110987771, 113987771, "r12"),
                (11, 87789108, 90766832, "r13"),
                (12, 109062195, 111562196, "r14"),
                (20, 33412194, 35912078, "r15")
            ]

            with open(out_dir / f'high-LD-regions_GRCH{self.built}.txt', 'w') as file:
                for line in data:
                    file.write(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\n")
            self.ld_regions = out_dir / f"high-LD-regions_GRCH{self.built}.txt"
            return out_dir / f'high-LD-regions_GRCH{self.built}.txt'
        else:
            logger.error(f"Unsupported genome build: {self.built}. Only '37' and '38' are supported.")
            return Path()
