import requests
import logging

import pandas as pd

from pathlib import Path

from ideal_genom_qc.Helpers import shell_do

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class Fetcher1000Genome:

    def __init__(self, destination: Path = None, built: str = '38'):

        if not destination:
            destination = Path(__file__).resolve().parent.parent / "data" / f"1000genomes_built_{built}"
        
        self.destination = destination
        self.built = built

        self.pgen_file = None
        self.pvar_file = None
        self.psam_file = None

        self.bed_file = None
        self.bim_file = None
        self.fam_file = None
        
        pass

    def get_1000genomes(self, url_pgen: str = None, url_pvar: str = None, url_psam: str = None)-> Path:

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
        
        self._download_file(url_pgen, self.destination / "all_phase3.pgen.zst")
        self.pvar_file = self._download_file(url_pvar, self.destination / "all_phase3.pvar.zst")
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
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            with open(destination, 'wb') as out_file:
                for chunk in response.iter_content(chunk_size=1024):
                    out_file.write(chunk)
            return destination
        except requests.RequestException as e:
            logging.error(f"Failed to download {url}: {e}")
            return None

    def get_1000genomes_binaries(self) -> Path:

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

        (self.destination / "all_phase3.bim").unlink(missing_ok=True)
        (self.destination / "all_phase3.bed").unlink(missing_ok=True)
        (self.destination / "all_phase3.fam").unlink(missing_ok=True)

        bim = (self.destination / "all_phase3_renamed.bim").with_name('all_phase3.bim')
        bed = (self.destination / "all_phase3_renamed.bed").with_name('all_phase3.bed')
        fam = (self.destination / "all_phase3_renamed.fam").with_name('all_phase3.fam')

        (self.destination / "all_phase3_renamed.bim").rename(bim)
        (self.destination / "all_phase3_renamed.bed").rename(bed)
        (self.destination / "all_phase3_renamed.fam").rename(fam)

        self.bed_file = bed
        self.bim_file = bim
        self.fam_file = fam

        return self.destination / "all_phase3"
    
    def _check_if_binaries_exist(self) -> bool:

        check_bed = (self.destination / "all_phase3.bed").exists()
        check_bim = (self.destination / "all_phase3.bim").exists()
        check_fam = (self.destination / "all_phase3.fam").exists()
        check_psam = (self.destination / "all_phase3.psam").exists()


        return check_bed and check_bim and check_fam and check_psam
    
class FetcherLDRegions:

    def __init__(self, destination: Path = Path(), built: str = '38'):

        if not destination:
            destination = Path(__file__).resolve().parent.parent / "data" / "ld_regions_files"

        self.destination = destination
        self.built = built

        self.ld_regions = None
        
        pass

    def get_ld_regions(self)-> Path:

        self.destination.mkdir(parents=True, exist_ok=True)

        if self.built == '37':
            url_ld_regions = r"https://raw.githubusercontent.com/genepi-freiburg/gwas/refs/heads/master/single-pca/high-LD-regions.txt"
        
            ld = requests.get(url_ld_regions)

            out_dir = self.destination

            if ld.status_code == 200:
                with open((out_dir / f"high-LD-regions_GRCh{self.built}.txt"), "wb") as f:
                    f.write(ld.content)
                logger.info(f"LD regions file for built {self.built} downloaded successfully to {out_dir}")
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
            return out_dir / f'high-LD-regions_GRCH{self.built}.txt'
