"""This module provides a class for performing Genome-Wide Association Studies (GWAS) using a Generalized Linear Model (GLM) with PLINK2. 

It includes methods for association analysis, obtaining top hits, and annotating SNPs with gene information.
"""
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from ..core.executor import run_plink2, run_gcta
from ..core.utils import get_available_memory, get_optimal_threads
from ideal_genom.utilities.annotations import annotate_snp

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

class GWAS_GLM:

    """Class for performing Genome-Wide Association Studies (GWAS) using a Generalized Linear Model (GLM) with PLINK2.

    This class provides methods to perform association analysis, obtain top hits, and annotate SNPs with gene information.

    Attributes
    ----------
    input_path : Path
        Path to the input directory.
    output_path : Path
        Path to the output directory.
    input_name : str 
        Base name of the input PLINK files.
    output_name : str 
        Base name for the output files.
    recompute : bool 
        Flag indicating whether to recompute the analysis.
    results_dir : Path
        Directory where the results will be saved.
        
    Raises
    ------
    ValueError
        If input_path, output_path, input_name, or output_name are not provided.
    FileNotFoundError
        If the specified input_path or output_path does not exist.
    FileNotFoundError
        If the required PLINK files (.bed, .bim, .fam) are not found in the input_path.
    TypeError
        If input_name or output_name are not strings, or if recompute is not a boolean.
    """

    def __init__(self, input_path: str | Path, input_name: str, output_path: str | Path, output_name: str, recompute: bool = True) -> None:
    
        # check if paths are set
        if input_path is None or output_path is None:
            raise ValueError("Values for input_path and output_path must be set upon initialization.")
        
        if not isinstance(input_path, (str, Path)) or not isinstance(output_path, (str, Path)):
            raise TypeError("input_path and output_path should be of type str or Path.")
        
        input_path = Path(input_path)
        output_path = Path(output_path)
        
        if not input_path.exists() or not input_path.is_dir():
            raise FileNotFoundError(f"Input path does not exist: {input_path}")
        if not output_path.exists() or not output_path.is_dir():
            raise FileNotFoundError(f"Output path does not exist: {output_path}")
        
        # check if input_name and output_name are set
        if input_name is None or output_name is None:
            raise ValueError("Values for input_name and output_name must be set upon initialization.")
        if not isinstance(input_name, str) or not isinstance(output_name, str):
            raise TypeError("input_name and output_name should be of type str.")
        
        # check existence of PLINK files
        if not (input_path / f'{input_name}.bed').exists():
            raise FileNotFoundError(f"PLINK bed file was not found: {input_path / f'{input_name}.bed'}")
        if not (input_path / f'{input_name}.bim').exists():
            raise FileNotFoundError(f"PLINK bim file was not found: {input_path / f'{input_name}.bim'}")
        if not (input_path / f'{input_name}.fam').exists():
            raise FileNotFoundError(f"PLINK fam file was not found: {input_path / f'{input_name}.fam'}")
        
        if not isinstance(recompute, bool):
            raise TypeError("recompute should be of type bool.")

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.recompute   = recompute

        # create results folder
        self.results_dir = output_path / 'gwas_glm'
        self.results_dir.mkdir(parents=True, exist_ok=True)

        logger.info("\033[1;32mGWAS fixed effects (GLM) analysis initialized.\033[0m")

    def glm_association_analysis(self, maf: float = 0.01, mind: float = 0.1, hwe: float = 5e-6, ci: float = 0.95) -> None:
        """Perform fixed model association analysis using PLINK2.
        
        This method performs a fixed model association analysis on genomic data using PLINK2. 
        It checks the validity of the input parameters, ensures necessary files exist, 
        and executes the PLINK2 command to perform the analysis.

        Parameters
        ----------
        maf : float
            Minor allele frequency threshold. Must be between 0 and 0.5.
        mind : float 
            Individual missingness threshold. Must be between 0 and 1.
        hwe : float 
            Hardy-Weinberg equilibrium threshold. Must be between 0 and 1.
        ci : float 
            Confidence interval threshold. Must be between 0 and 1.
        
        Returns
        -------
        dict
            A dictionary containing the status of the process, the step name, and the output directory.
        
        Raises
        ------
        TypeError
            If any of the input parameters are not of type float.
        ValueError
            If any of the input parameters are out of their respective valid ranges.
        FileNotFoundError
            If the required PCA file is not found.
        """

        output_name= self.output_name
        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir
        recompute  = self.recompute

        logger.info("Starting association analysis with GLM.")

        # Check type of maf
        if not isinstance(maf, float):
             raise TypeError("maf should be of type float.")

        # Check type of mind
        if not isinstance(mind, float):
            raise TypeError("mind should be of type float.")
        
        # Check type of hwe
        if not isinstance(hwe, float):
            raise TypeError("hwe should be of type float.")
        
        # Check type of ci
        if not isinstance(ci, float):
            raise TypeError("ci should be of type float.")
        
        # Check if maf is in range
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be between 0 and 0.5")
        
        # Check if mind is in range
        if mind < 0 or mind > 1:
            raise ValueError("mind should be between 0 and 1")
        
        # Check if hwe is in range
        if hwe < 0 or hwe > 1:
            raise ValueError("hwe should be between 0 and 1")
        
        # Check if ci is in range
        if ci <= 0 or ci >= 1:
            raise ValueError("ci should be between 0 and 1")
        
        # check if the PCA file exists
        if not (input_path / f'{input_name}.eigenvec').exists():
            raise FileNotFoundError(f"PCA file was not found: {input_path / f'{input_name}.eigenvec'}")

        threads = get_optimal_threads()
        memory = get_available_memory()

        if recompute:

            # plink2 command to perform association analysis
            run_plink2([
                '--bfile', str(input_path / input_name),
                '--adjust',
                '--ci', str(ci),
                '--maf', str(maf),
                '--mind', str(mind),
                '--hwe', str(hwe),
                '--covar', str(input_path / f'{input_name}.eigenvec'),
                '--glm', 'hide-covar', 'omit-ref', 'sex', 'cols=+a1freq,+beta',
                '--out', str(results_dir / f'{output_name}_glm'),
                '--threads', str(threads),
                '--memory', str(memory)
            ])

        logger.info("Association analysis with GLM completed.")

        return

    def get_top_hits(self, maf: float = 0.01) -> None:
        
        """Get the top hits from the GWAS results.

        Parameters
        ----------
        maf : float
            Minor allele frequency threshold. Must be a float between 0 and 0.5.

        Returns
        -------
        dict
            A dictionary containing the process status, step name, and output directory.

        Raises
        ------
        TypeError
            If maf is not of type float.
        ValueError
            If maf is not between 0 and 0.5.

        Notes
        -----
        The function performs the following steps:
            1. Validates the type and range of the maf parameter.
            2. Computes the number of threads to use based on the available CPU cores.
            3. Loads the results of the association analysis and renames columns according to GCTA requirements.
            4. Prepares a .ma file with the necessary columns.
            5. If recompute is True, constructs and executes a GCTA command to perform conditional and joint analysis.
            6. Returns a dictionary with the process status, step name, and output directory.
        """

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name
        recompute   = self.recompute

        # check type and range of maf
        if not isinstance(maf, float):
            raise TypeError("maf should be of type float.")
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be between 0 and 0.5")

        logger.info("Getting top hits from GWAS results.")

        threads = get_optimal_threads()

        # load results of association analysis and rename columns according GCTA requirements
        df = pd.read_csv(results_dir / f'{output_name}_glm.PHENO1.glm.logistic.hybrid', sep="\t")
        rename = {
            '#CHROM'          : 'CHR',	
            'POS'             : 'POS',	
            'ID'              : 'SNP',
            'REF'             : 'A2',	
            'ALT'             : 'ALT',	
            'PROVISIONAL_REF?': 'PROVISIONAL_REF',	
            'A1'              : 'A1',	
            'OMITTED'         : 'OMITTED',	
            'A1_FREQ'         : 'freq',	
            'FIRTH?'          : 'FIRTH',	
            'TEST'            : 'TEST',	
            'OBS_CT'          : 'N',	
            'BETA'            : 'b',	
            'SE'              : 'se',	
            'L95'             : 'L95',	
            'U95'             : 'U95',	
            'Z_STAT'          : 'Z_STAT',	
            'P'               : 'p',	
            'ERRCODE'         : 'ERRCODE'
        }
        df = df.rename(columns=rename)

        # prepare .ma file
        df = df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].copy()

        df.to_csv(results_dir / 'cojo_file.ma', sep="\t", index=False)

        del df

        if recompute:
            # gcta command for conditional analysis
            run_gcta([
                '--bfile', str(input_path / input_name),
                '--maf', str(maf),
                '--cojo-slct',
                '--cojo-file', str(results_dir / 'cojo_file.ma'),
                '--out', str(results_dir / f'{output_name}-cojo'),
                '--thread-num', str(threads)
            ])

        logger.info("Top hits extraction completed.")

        return
    
    def annotate_top_hits(self, gtf_path: Optional[str] = None, build: str = '38', anno_source: str = "ensembl") -> None:
        """Annotate top SNP hits from COJO analysis with gene information.
        
        This method reads the COJO joint analysis results, extracts the top SNPs, 
        and annotates them with gene information using the specified genome build 
        and annotation source. The annotated results are saved to a TSV file.
        
        Parameters
        ----------
        gtf_path : Optional[str], default=None
            Path to the GTF (Gene Transfer Format) file for custom annotation.
            If None, the annotation will use default resources.
        build : str, default='38'
            Genome build version to use for annotation ('38' for GRCh38, etc.).
        anno_source : str, default="ensembl"
            Source of annotations to use (e.g., "ensembl", "refseq").
        
        Returns
        -------
        dict
            A dictionary containing:
            - 'pass': Boolean indicating if the process completed successfully
            - 'step': The name of the step ('annotate_hits')
            - 'output': Dictionary with output file paths
        
        Raises
        ------
        FileExistsError
            If the COJO results file is not found in the results directory.
        
        Notes
        -----
        The annotated results are saved to 'top_hits_annotated.tsv' in the results directory.
        """

        results_dir = self.results_dir

        logger.info("Starting annotation of top hits.")

        # load the data
        cojo_file_path = results_dir / (self.output_name + '-cojo.jma.cojo')
        if cojo_file_path.exists():
            df_hits = pd.read_csv(cojo_file_path, sep="\t")
        else:
            raise FileNotFoundError(f"File {self.output_name}-cojo.jma.cojo not found in the results directory: {results_dir}")
        df_hits = df_hits[['Chr', 'SNP', 'bp']].copy()

        if not df_hits.empty:
            df_hits = annotate_snp(
                df_hits,
                chrom  ='Chr',
                pos    ='bp',
                build  =build,
                source =anno_source,
                gtf_path=gtf_path # type: ignore
            ).rename(columns={"GENE":"GENENAME"})

        df_hits.to_csv(results_dir / 'top_hits_annotated.tsv', sep="\t", index=False)
        
        logger.info("Top hits annotation completed and saved.")
        
        return

    def execute_gwas_glm_pipeline(self, glm_params: dict) -> None:
        """Execute the complete GWAS fixed effects pipeline.

        This method orchestrates the full GWAS analysis workflow using a generalized linear
        model (GLM). It sequentially executes all necessary steps: performing the association
        analysis, extracting top hits, and annotating them with gene information.

        Parameters
        ----------
        maf : float, optional (default=0.01)
            Minor allele frequency threshold for filtering SNPs. Must be between 0 and 0.5.
        mind : float, optional (default=0.1)
            Individual missingness threshold. Must be between 0 and 1.
        hwe : float, optional (default=5e-6)
            Hardy-Weinberg equilibrium threshold. Must be between 0 and 1.
        ci : float, optional (default=0.95)
            Confidence interval threshold. Must be between 0 and 1.
        gtf_path : Optional[str], optional (default=None)
            Path to the GTF file for custom annotation. If None, uses default annotation resources.
        build : str, optional (default='38')
            Genome build version to use for annotation ('38' for GRCh38, '37' for GRCh37).
        anno_source : str, optional (default='ensembl')
            Source of annotations to use (e.g., 'ensembl', 'refseq').

        Returns
        -------
        None
            Results are saved to the results directory specified during initialization.

        Raises
        ------
        TypeError
            If any of the numeric parameters are not of the correct type.
        ValueError
            If any of the parameters are out of their respective valid ranges.
        FileNotFoundError
            If required input files (e.g., PCA file) are not found.

        Notes
        -----
        This method sequentially calls:
        1. glm_association_analysis() - Performs the GLM-based GWAS analysis
        2. get_top_hits() - Extracts top significant hits using conditional analysis
        3. annotate_top_hits() - Annotates hits with gene information

        The pipeline expects that preparatory steps (LD pruning and PCA) have already been
        performed, as it requires the existence of .eigenvec files.

        Examples
        --------
        >>> gwas = GWASfixed(input_path='data/', input_name='genotypes',
        ...                  output_path='results/', output_name='gwas_results')
        >>> gwas.execute_gwasfixed_pipeline(maf=0.05, hwe=1e-6, build='38')
        """

        maf         = glm_params.get('maf', 0.01)
        mind        = glm_params.get('mind', 0.1)
        hwe         = glm_params.get('hwe', 5e-6)
        ci          = glm_params.get('ci', 0.95)
        gtf_path    = glm_params.get('gtf_path', None)
        build       = glm_params.get('build', '38')
        anno_source = glm_params.get('anno_source', 'ensembl')

        logger.info("\033[1;34m" + "="*80 + "\033[0m")
        logger.info("\033[1;34mStarting GWAS Fixed Effects (GLM) Pipeline\033[0m")
        logger.info("\033[1;34m" + "="*80 + "\033[0m")

        # Step 1: Run association analysis with GLM
        logger.info("\033[1;33m[1/3] Running GLM association analysis...\033[0m")
        self.glm_association_analysis(maf=maf, mind=mind, hwe=hwe, ci=ci)
        logger.info("\033[1;32m✓ Association analysis completed\033[0m\n")

        # Step 2: Extract top hits
        logger.info("\033[1;33m[2/3] Extracting top hits from GWAS results...\033[0m")
        self.get_top_hits(maf=maf)
        logger.info("\033[1;32m✓ Top hits extracted\033[0m\n")

        # Step 3: Annotate top hits
        logger.info("\033[1;33m[3/3] Annotating top hits with gene information...\033[0m")
        self.annotate_top_hits(gtf_path=gtf_path, build=build, anno_source=anno_source)
        logger.info("\033[1;32m✓ Annotation completed\033[0m\n")

        logger.info("\033[1;34m" + "="*80 + "\033[0m")
        logger.info("\033[1;32m✓ GWAS Fixed Effects Pipeline Completed Successfully!\033[0m")
        logger.info("\033[1;34m" + "="*80 + "\033[0m")
        logger.info(f"\033[1mResults saved to: {self.results_dir}\033[0m\n")

        return