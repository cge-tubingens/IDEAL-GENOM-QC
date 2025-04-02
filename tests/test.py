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
            logger.info(f".bim file downloaded successfully to {out_dir}")
        else:
            logger.info(f"Failed to download .bim file: {bim.status_code}")

        if fam.status_code == 200:
            with open((out_dir / "example.fam"), "wb") as f:
                f.write(fam.content)
            logger.info(f".fam file downloaded successfully to {out_dir}")
        else:
            logger.info(f"Failed to download .fam file: {fam.status_code}")

        if bed.status_code == 200:
            with open((out_dir / "example.bed"), "wb") as f:
                f.write(bed.content)
            logger.info(f".bed file downloaded successfully to {out_dir}")
        else:
            logger.info(f"Failed to download .bed file: {bed.status_code}")

        return
    
class PrepareExampleProject():

    def __init__(self, destination: Path = Path(__file__).resolve().parent.parent / "data" / "example_data", built: str = '38'):

        self.destination = destination
        self.built = built

        self.input_file = self.destination / "inputData"
        self.input_file.mkdir(parents=True, exist_ok=True)

        pass

    def prepare_folders(self) -> None:

        if not self._check_binary_files():
            logger.info("Downloading example binary files.")
            ExampleDataFetcher(destination=self.destination).get_sample_binary()
        
        self.output_file = self.destination / "outputData"
        self.output_file.mkdir(parents=True, exist_ok=True)

        self.dependables = self.destination / "dependables"
        self.dependables.mkdir(parents=True, exist_ok=True)

        return
    
    def _get_ld_region_file(self) -> None:

        ld_file = self.dependables / f"high-LD-regions_GRCh{self.built}.txt"

        if not ld_file.exists():
            logger.info("LD regions file not found.")
            logger.info("Downloading LD regions file.")
            fetcher = FetcherLDRegions(built=self.built)
            fetcher.get_ld_regions()
            return

        return
        
    def _check_binary_files(self) -> bool:

        bim = self.input_file / "example.bim"
        bed = self.input_file / "example.bed"
        fam = self.input_file / "example.fam"

        if bim.exists() and bed.exists() and fam.exists():
            return True
        else:
            return False
        