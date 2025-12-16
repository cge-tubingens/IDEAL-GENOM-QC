"""This module provides a class for performing Genome-Wide Association Studies (GWAS) using a Generalized Linear Mixed Model (GLMM) with PLINK2. 

It includes methods for association analysis, obtaining top hits, and annotating SNPs with gene information.
"""
from pathlib import Path
import pandas as pd

from ..core.executor import run_gcta
from ..core.utils import get_optimal_threads
from ideal_genom.utilities.annotations import annotate_snp

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class GWASrandom:

    """Class for performing Genome-Wide Association Studies (GWAS) using a Generalized Linear Mixed Model (GLM) with PLINK2.

    This class provides methods to perform association analysis, obtain top hits, and annotate SNPs with gene information.
        
    Parameters
    ----------
    input_path : str 
        Path to the input directory containing PLINK files.
    input_name : str 
        Base name of the input PLINK files (without extensions).
    output_path : str 
        Path to the output directory where results will be saved.
    output_name : str 
        Base name for the output files.
    recompute : bool 
        Flag indicating whether to recompute the analysis if results already exist. Default is True.
    
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

    def __init__(self, input_path:Path, input_name:str, output_path:Path, output_name:str, recompute:bool = True) -> None:

       # check if paths are set
        if input_path is None or output_path is None:
            raise ValueError("Values for input_path, output_path and dependables_path must be set upon initialization.")
        
        # check if input_name and output_name are set
        if input_name is None or output_name is None:
            raise ValueError("Values for input_name and output_name must be set upon initialization.")
        if not isinstance(input_name, str) or not isinstance(output_name, str):
            raise TypeError("input_name and output_name should be of type str.")
        
        # check existence of PLINK files
        if not (input_path / f'{input_name}.bed').exists():
            raise FileNotFoundError(f"PLINK bed file was not found: {(input_path / f'{input_name}.bed')}")
        if not (input_path / f'{input_name}.bim').exists():
            raise FileNotFoundError(f"PLINK bim file was not found: {(input_path / f'{input_name}.bim')}")
        if not (input_path / f'{input_name}.fam').exists():
            raise FileNotFoundError(f"PLINK fam file was not found: {(input_path / f'{input_name}.fam')}")
        
        if not isinstance(recompute, bool):
            raise TypeError("recompute should be of type bool.")
        
        if not input_path.exists():
            raise FileNotFoundError(f"Input path does not exist: {input_path}")
        if not output_path.exists():
            raise FileNotFoundError(f"Output path does not exist: {output_path}")

        self.input_path  = input_path
        self.output_path = output_path
        self.input_name  = input_name
        self.output_name = output_name
        self.recompute   = recompute

        # create results folder
        self.results_dir = output_path / 'gwas_random'
        if not self.results_dir.exists():
            self.results_dir.mkdir(parents=True)

        print("\033[1;32mAnalysis of GWAS data using a random effect model initialized.\033[0m")

        pass

    def prepare_aux_files(self) -> None:
        
        """Prepares auxiliary files for GWAS analysis by processing phenotype and sex data.
        
        This function reads a .fam file, extracts and recodes phenotype and sex information, 
        and writes the processed data to new files in the specified results directory.

        Returns
        -------
        dict
            A dictionary containing the status of the process, the step name, and the output directory.
        """

        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir

        step = "prepare_aux_files"

        df_fam = pd.read_csv(
            input_path / f'{input_name}.fam',
            sep   =r'\s+',
            engine='python',
            header=None
        )

        df_pheno = df_fam[[df_fam.columns[0], df_fam.columns[1], df_fam.columns[5]]].copy()

        # recode phenotype
        df_pheno[5] = df_pheno[5]-1

        df_pheno.to_csv(
            results_dir / f'{input_name}_pheno.phen',
            sep   ='\t',
            header=False,
            index =False
        )

        # recode sex
        df_sex = df_fam[[0,1,4]].copy()

        df_sex.to_csv(
            results_dir / f'{input_name}_sex.covar',
            sep   ='\t',
            header=False,
            index =False
        )

        return
    
    def compute_grm(self, max_threads: Optional[int] = None) -> None:
        
        """Compute the Genetic Relationship Matrix (GRM) using GCTA software.

        This method computes the GRM for the given input data using the GCTA software. 
        It allows for multi-threaded execution and can optionally recompute the GRM if specified.

        Parameters
        ----------
        max_threads : int, optional
            The maximum number of threads to use for computation. If not specified, 
            it defaults to the number of available CPU cores minus two. 
            If the number of CPU cores cannot be determined, it defaults to 10.

        Returns
        -------
        dict 
            A dictionary containing the following keys:
                - 'pass' (bool): Indicates whether the process completed successfully.
                - 'step' (str): The name of the step performed ('compute_grm').
                - 'output' (dict): A dictionary containing the output file paths with the key 'gcta_out'.
        """

        results_dir= self.results_dir
        input_path = self.input_path
        input_name = self.input_name
        recompute  = self.recompute

        step = "compute_grm"

        # compute the number of threads to use
        threads = max_threads or get_optimal_threads()

        if recompute:
                # gcta commands as lists
                gcta_args1 = [
                    '--bfile', str(input_path / f'{input_name}-pruned'),
                    '--make-grm',
                    '--thread-num', str(threads),
                    '--out', str(results_dir / f'{input_name}_grm')
                ]
                gcta_args2 = [
                    '--grm', str(results_dir / f'{input_name}_grm'),
                    '--make-bK-sparse', '0.05',
                    '--out', str(results_dir / f'{input_name}_sparse')
                ]

                # run gcta commands
                run_gcta(gcta_args1)
                run_gcta(gcta_args2)
        return
    
    def run_gwas_glmm(self, maf: float = 0.01) -> None:
        
        """Runs a Genome-Wide Association Study (GWAS) using a generalized linear mixed model (GLMMM).

        Parameters
        ----------
        maf : float 
            Minor allele frequency threshold for filtering SNPs. Default is 0.01.

        Returns
        -------
        dict
            A dictionary containing the status of the process, the step name, and the output directory.

        Raises
        ------
        TypeError 
            If `maf` is not of type float.
        ValueError 
            If `maf` is not between 0 and 1.
        FileExistsError 
            If required input files are not found in the results directory.
        """

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name

        if not isinstance(maf, float):
            raise TypeError("maf should be of type float.")
        if maf < 0 or maf > 1:
            raise ValueError("maf should be between 0 and 1.")

        

        # compute the number of threads to use
        threads = get_optimal_threads()

        if not (results_dir / f'{input_name}_sparse.grm.id').exists():
            raise FileExistsError(f"File {input_name+'_sparse.grm.id'} is not in the results directory.")
        if not (results_dir / f'{input_name}_sparse.grm.sp').exists():
            raise FileExistsError(f"File {input_name+'_sparse.grm.id'} is not in the results directory.")
        if not (results_dir / f'{input_name}_sex.covar').exists():
            raise FileExistsError(f"File {input_name+'sex.covar'} is not in the results directory.")
        if not (results_dir / f'{input_name}_pheno.phen').exists():
            raise FileExistsError(f"File {input_name+'_pheno.phen'} is not in the results directory.")

        # gcta command
        gcta_args = [
            '--bfile', str(input_path / input_name),
            '--fastGWA-mlm-binary',
            '--maf', str(maf),
            '--grm-sparse', str(results_dir / f'{input_name}_sparse'),
            '--qcovar', str(input_path / f'{input_name}.eigenvec'),
            '--covar', str(results_dir / f'{input_name}_sex.covar'),
            '--pheno', str(results_dir / f'{input_name}_pheno.phen'),
            '--out', str(results_dir / f'{output_name}_assocSparseCovar_pca_sex-mlm-binary'),
            '--thread-num', str(threads)
        ]

        # run gcta command
        run_gcta(gcta_args)

        return
    
    def get_top_hits(self, maf: float = 0.01) -> None:
        
        """Get the top hits from the GWAS results.

        This function processes the results of a genome-wide association study (GWAS) 
        to identify the top hits based. It prepares the necessary files and 
        optionally recomputes the results using GCTA.

        Parameters
        ----------
        maf : float, optional 
            Minor allele frequency threshold. Default is 0.01. Must be between 0 and 1.

        Returns
        -------
        dict
            A dictionary containing the status of the process, the step name, and the output directory.

        Raises
        ------
        TypeError 
            If `maf` is not of type float.
        ValueError 
            If `maf` is not between 0 and 0.5.
        """

        results_dir = self.results_dir
        input_name  = self.input_name
        input_path  = self.input_path
        output_name = self.output_name
        recompute   = self.recompute

        if not isinstance(maf, float):
            raise TypeError("maf should be of type float.")
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be between 0 and 1.")

        step = "get_top_hits"

        # compute the number of threads to use
        threads = get_optimal_threads()

        # load results of association analysis and rename columns according to GCTA requirements
        df = pd.read_csv(results_dir / f'{output_name}_assocSparseCovar_pca_sex-mlm-binary.fastGWA', sep="\t")
        rename = {
            'CHR'     :'CHR',	
            'SNP'     :'SNP',
            'POS'     :'POS',	
            'A1'      :'A1', 
            'A2'      :'A2', 
            'N'       :'N', 
            'AF1'     :'freq', 
            'T'       :'T', 
            'SE_T'    :'SE_T', 
            'P_noSPA' :'P_noSPA',
            'BETA'    :'b', 
            'SE'      :'se', 
            'P'       :'p', 
            'CONVERGE':'CONVERGE'
        }
        df = df.rename(columns=rename)        

        # prepare .ma file
        df = df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].copy()

        df.to_csv(results_dir / 'cojo_file.ma', sep="\t", index=False)

        del df

        if recompute:
            # gcta command as list
            gcta_args = [
                '--bfile', str(input_path / input_name),
                '--maf', str(maf),
                '--cojo-slct',
                '--cojo-file', str(results_dir / 'cojo_file.ma'),
                '--out', str(results_dir / f'{output_name}_assocSparseCovar_pca_sex-mlm-binary-cojo'),
                '--thread-num', str(threads)
            ]

            # execute gcta command
            run_gcta(gcta_args)

        return

    def annotate_top_hits(self, gtf_path: Optional[str] = None, build: str = '38', anno_source: str = 'ensembl') -> None:
        """Annotate top genetic hits from GWAS analysis with gene information.
        
        This method loads top hits from COJO analysis results, annotates them with gene information
        using the specified genome build and annotation source, and saves the annotated results 
        to a TSV file.
        
        Parameters
        ----------
        gtf_path : Optional[str], default=None
            Path to a GTF file for custom annotation. If None, will use built-in annotation resources.
        build : str, default='38'
            Genome build version to use for annotation (e.g., '38', '37').
        anno_source : str, default='ensembl'
            Source of the annotation data (e.g., 'ensembl').
        
        Returns
        -------
        dict
            A dictionary containing:
            - 'pass': bool - Whether the process completed successfully
            - 'step': str - The name of the processing step
            - 'output': dict - Dictionary of output file paths
        
        Raises
        ------
        FileExistsError
            If the COJO file is not found in the results directory.
        
        Notes
        -----
        The annotated results are saved to 'top_hits_annotated.tsv' in the results directory.
        """

        results_dir = self.results_dir
        output_name = self.output_name

        step = "annotate_hits"

        # load the data
        cojo_file_path = results_dir / f'{output_name}_assocSparseCovar_pca_sex-mlm-binary-cojo.jma.cojo'

        # check if .jma file exists
        if cojo_file_path.exists():
            df_hits = pd.read_csv(cojo_file_path, sep="\t")
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
        else:
            raise FileExistsError("File cojo_file.jma not found in the results directory.")

        
        return
