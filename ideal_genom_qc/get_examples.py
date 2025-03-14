import requests
import logging

from pathlib import Path

from ideal_genom_qc.get_references import FetcherLDRegions

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class ExampleDataFetcher:

    def __init__(self, destination: Path = Path(__file__).resolve().parent.parent / "data" / "example_data"):

        self.destination = destination

        self.input_file = self.destination / "inputData"
        self.input_file.mkdir(parents=True, exist_ok=True)

        self.bim = None
        self.bed = None
        self.fam = None

        pass

    def get_sample_binary(self) -> None:

        url_bim = r"https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/refs/heads/master/inst/extdata/data.bim"
        url_bed = r"https://github.com/meyer-lab-cshl/plinkQC/raw/refs/heads/master/inst/extdata/data.bed"
        url_fam = r"https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/refs/heads/master/inst/extdata/data.fam"

        bim = requests.get(url_bim)
        bed = requests.get(url_bed)
        fam = requests.get(url_fam)

        out_dir = self.input_file

        if bim.status_code == 200:
            with open((out_dir / "example.bim"), "wb") as f:
                f.write(bim.content)
            logger.info(f"example.bim file downloaded successfully to {out_dir}")
            self.bim = out_dir / "example.bim"
        else:
            logger.info(f"Failed to download example.bim file: {bim.status_code}")

        if fam.status_code == 200:
            with open((out_dir / "example.fam"), "wb") as f:
                f.write(fam.content)
            logger.info(f"example.fam file downloaded successfully to {out_dir}")
            self.fam = out_dir / "example.fam"
        else:
            logger.info(f"Failed to download example.fam file: {fam.status_code}")

        if bed.status_code == 200:
            with open((out_dir / "example.bed"), "wb") as f:
                f.write(bed.content)
            logger.info(f"example.bed file downloaded successfully to {out_dir}")
            self.bed = out_dir / "example.bed"
        else:
            logger.info(f"Failed to download example.bed file: {bed.status_code}")

        return
    
class PrepareExampleProject():

    def __init__(self, destination: Path = Path(__file__).resolve().parent.parent / "data" / "example_data", built: str = '38'):

        self.destination = destination
        self.built = built

        self.input_path = self.destination / "inputData"
        self.input_path.mkdir(parents=True, exist_ok=True)

        pass

    def prepare_folders(self) -> None:

        if not self._check_binary_files():
            logger.info("Downloading example binary files.")
            ExampleDataFetcher(destination=self.destination).get_sample_binary()
        else:
            logger.info("Example binary files already exist.")
            self.bim = self.input_path / "example.bim"
            self.bed = self.input_path / "example.bed"
            self.fam = self.input_path / "example.fam"
        
        self.output_path = self.destination / "outputData"
        self.output_path.mkdir(parents=True, exist_ok=True)

        self.dependables = self.destination / "dependables"
        self.dependables.mkdir(parents=True, exist_ok=True)

        logger.info(f'Fetching high LD regions file for built {self.built}.')
        self.ld_file  = self._get_ld_region_file()

        return
    
    def _get_ld_region_file(self) -> Path:

        ld_file = self.dependables / f"high-LD-regions_GRCh{self.built}.txt"

        if not ld_file.exists():
            logger.info("LD regions file not found.")
            logger.info("Downloading LD regions file.")
            fetcher = FetcherLDRegions(destination=self.dependables, built=self.built)
            ld_file = fetcher.get_ld_regions()
            return ld_file

        return ld_file
        
    def _check_binary_files(self) -> bool:

        bim = self.input_path / "example.bim"
        bed = self.input_path / "example.bed"
        fam = self.input_path / "example.fam"

        if bim.exists() and bed.exists() and fam.exists():
            return True
        else:
            return False
        