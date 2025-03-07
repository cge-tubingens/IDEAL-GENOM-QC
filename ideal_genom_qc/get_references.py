import requests
import logging

from pathlib import Path

from ideal_genom_qc.Helpers import shell_do

class Fetcher1000Genome:

    def __init__(self, destination: Path = Path(__file__).resolve().parent.parent / "data" / "1000genomes"):

        self.destination = destination

        self.pgen_file = None
        self.pvar_file = None
        self.psam_file = None
        
        pass

    def get_1000genomes(self, url_pgen: str = None, url_pvar: str = None, url_psam: str = None)-> Path:

        self.destination.mkdir(parents=True, exist_ok=True)

        if url_pgen is None:
            url_pgen = r"https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1"
        if url_pvar is None:
            url_pvar = r"https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1"
        if url_psam is None:
            url_psam = r"https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1"

        if self._check_if_binaries_exist():

            self.pvar_file = self.destination / "all_phase3.pvar.zst"
            self.psam_file = self.destination / "all_phase3.psam"
            self.pgen_decompressed = self.destination / "all_phase3.pgen"

            return self.pgen_decompressed
        
        self._download_file(url_pgen, self.destination / "all_phase3.pgen.zst")
        self.pvar_file = self._download_file(url_pvar, self.destination / "all_phase3.pvar.zst")
        self.psam_file = self._download_file(url_psam, self.destination / "all_phase3.psam")

        pgen_file = self.destination / "all_phase3.pgen.zst"
        pgen_decompressed = self.destination / "all_phase3.pgen"

        # plink2 command
        plink2_cmd = f"plink2 --zst-decompress {str(pgen_file)} {str(pgen_decompressed)}"

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

            (self.destination / "all_phase3.pgen").unlink(missing_ok=True)
            (self.destination / "all_phase3.pgen.zst").unlink(missing_ok=True)
            (self.destination / "all_phase3.pvar.zst").unlink(missing_ok=True)

            return self.destination / "all_phase3"

        # plink2 command
        plink2_cmd = f"plink2 --pfile {str(self.destination / "all_phase3")} vzs --max-alleles 2 --make-bed --out {str(self.destination / "all_phase3")}"

        # execute plink2 command
        shell_do(plink2_cmd)

        (self.destination / "all_phase3.pgen").unlink()
        (self.destination / "all_phase3.pgen.zst").unlink()
        (self.destination / "all_phase3.pvar.zst").unlink()

        return self.destination / "all_phase3"
    
    def _check_if_binaries_exist(self) -> bool:

        check_bed = (self.destination / "all_phase3.bed").exists()
        check_bim = (self.destination / "all_phase3.bim").exists()
        check_fam = (self.destination / "all_phase3.fam").exists()
        check_psam = (self.destination / "all_phase3.psam").exists()


        return check_bed and check_bim and check_fam and check_psam

