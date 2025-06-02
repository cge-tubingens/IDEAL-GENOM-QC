"""
Module to draw plots based on UMAP dimension reduction
"""

import os
import umap
import warnings
import logging
import psutil
from typing import Optional

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path

from ideal_genom_qc.Helpers import shell_do
from ideal_genom_qc.get_references import FetcherLDRegions, Fetcher1000Genome
from ideal_genom_qc.AncestryQC import ReferenceGenomicMerger
from sklearn.model_selection import ParameterGrid

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class UMAPplot:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, high_ld_file: Path=Path(), built: str = '38', recompute_pca: bool = True) -> None:
        """
        Initialize UMAPplot object for population structure analysis.
        This class handles the creation of UMAP plots for genetic data, managing input/output paths
        and configuration for the analysis.

        Parameters
        ----------
        input_path : Path
            Path to the directory containing input files
        input_name : str
            Name of the input file
        output_path : Path
            Path to the directory where results will be saved
        output_name : str
            Name for the output files
        high_ld_file : Path
            Path to the file containing high LD regions
        built : str, optional
            Genome build version, either '37' or '38' (default is '38')
        recompute_pca : bool, optional
            Whether to recompute PCA analysis (default is True)

        Raises
        ------
        TypeError
            If input types are incorrect for any parameter
        ValueError
            If built is not '37' or '38'
        FileNotFoundError
            If input_path or output_path do not exist
        
        Notes
        -----
        If high_ld_file is not found, it will be automatically fetched from the package.
        Creates 'umap_results' and 'plots' directories in the output path.
        """

        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(high_ld_file, Path):
            raise TypeError("high_ld_regions should be a Path object")
        if not isinstance(input_name, str): 
            raise TypeError("input_name should be a string")
        if not isinstance(recompute_pca, bool):
            raise TypeError("recompute_merge should be a boolean")
        if not isinstance(built, str):
            raise TypeError("built should be a string")
        if built not in ['37', '38']:
            raise ValueError("built should be either '37' or '38'")        
        if not input_path.exists():
            raise FileNotFoundError("input_path does not exist")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")
        if not high_ld_file.is_file():
            logger.info(f"High LD file not found at {high_ld_file}")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions()
            ld_fetcher.get_ld_regions()

            if ld_fetcher.ld_regions is None:
                raise FileNotFoundError("Could not fetch LD regions file.")
                
            high_ld_file = ld_fetcher.ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_file}")

        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.high_ld_regions = high_ld_file
        self.recompute_pca = recompute_pca
        self.built = built

        self.files_to_keep= []

        self.results_dir = self.output_path / 'umap_results' 
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.plots_dir = self.results_dir / 'plots'
        self.plots_dir.mkdir(parents=True, exist_ok=True)

        pass

    def ld_pruning(self, maf: float = 0.001, geno: float = 0.1, mind: float = 0.2, hwe: float = 5e-8, ind_pair: list = [50, 5, 0.2]) -> None:
        """
        Perform Linkage Disequilibrium (LD) pruning on genetic data using PLINK.
        This method filters SNPs based on specified thresholds for various quality control metrics
        and performs LD-based pruning to remove highly correlated variants.
        
        Parameters
        ----------
        maf : float, default=0.001
            Minor allele frequency threshold. Variants with MAF below this value are removed.
            Must be between 0 and 0.5.
        geno : float, default=0.1
            Maximum per-SNP missing rate. Variants with missing rate above this are removed.
            Must be between 0 and 1.
        mind : float, default=0.2
            Maximum per-individual missing rate. Samples with missing rate above this are removed.
            Must be between 0 and 1. Recommended range is 0.02 to 0.1.
        hwe : float, default=5e-8
            Hardy-Weinberg equilibrium exact test p-value threshold.
            Variants with p-value below this are removed. Must be between 0 and 1.
        ind_pair : list, default=[50, 5, 0.2]
            Parameters for pairwise LD pruning: [window size, step size, r² threshold].
        
        Returns
        -------
        None
            Creates pruned PLINK binary files in the results directory.
        
        Notes
        -----
        - Skips processing if recompute_pca is False
        - Uses multithreading with optimal thread count based on system CPU
        - Generates intermediate files: .prune.in and .prune.out
        - Creates final LD-pruned dataset with '-LDpruned' suffix
        
        Raises
        ------
        TypeError
            If input parameters are not of type float
        ValueError
            If input parameters are outside their valid ranges
        UserWarning
            If mind parameter is outside recommended range
        """

        if not self.recompute_pca:
            logger.info(f"`recompuite_pca` is set to {self.recompute_pca}. LD pruning will be skipped.")
            logger.info("LD pruning already performed. Skipping this step.")
            return

        # Check type of maf
        if not isinstance(maf, float):
             raise TypeError("maf should be of type float.")

        # Check type of geno
        if not isinstance(geno, float):
            raise TypeError("geno should be of type float.")

        # Check type of mind
        if not isinstance(mind, float):
            raise TypeError("mind should be of type float.")
        
        # Check type of hwe
        if not isinstance(hwe, float):
            raise TypeError("hwe should be of type float.")
        
        # Check if maf is in range
        if maf <= 0.0 or maf >= 0.5:
            raise ValueError("maf should be between 0 and 0.5")
        
        # Check if geno is in range
        if geno <= 0 or geno >= 1:
            raise ValueError("geno should be between 0 and 1")
        
        # Check if mind is in range
        if mind < 0 or mind > 1:
            raise ValueError("mind should be between 0 and 1")
        
        # Check if mind is around typical values
        if mind <= 0.02 and mind >= 0.1:
            warnings.warn(f"The 'mind' value {mind} is outside the recommended range of 0.02 to 0.1.", UserWarning)

        # Check if hwe is in range
        if hwe < 0 or hwe > 1:
            raise ValueError("hwe should be between 0 and 1")

        logger.info("Executing LD pruning with the following parameters:")
        logger.info(f"LD pruning parameters: maf={maf}, geno={geno}, mind={mind}, hwe={hwe}, ind_pair={ind_pair}")

        cpu_count = os.cpu_count()
        if cpu_count is not None:
            max_threads = max(1, cpu_count - 2)
        else:
            # Dynamically calculate fallback as half of available cores or default to 2
            max_threads = max(1, (psutil.cpu_count(logical=True) or 2) // 2)

        # generates prune.in and prune.out files
        plink_cmd1 = f"plink --bfile {self.input_path / self.input_name} --maf {maf} --geno {geno} --mind {mind} --hwe {hwe} --exclude {self.high_ld_regions} --range --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --threads {max_threads} --out {self.results_dir / self.input_name}"

        # prune and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {self.input_path / self.input_name} --keep-allele-order --extract {self.results_dir / (self.input_name+'.prune.in')} --make-bed --threads {max_threads} --out {self.results_dir / (self.input_name+'-LDpruned')}"

        # execute plink command
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        return

    def compute_pcas(self, pca: int = 10) -> None:
        """
        Computes Principal Component Analysis (PCA) using PLINK.

        This method performs PCA on the LD-pruned dataset using PLINK's --pca command.
        The analysis generates eigenvectors and eigenvalues that can be used for
        population structure analysis and visualization.

        Parameters
        ----------
        pca : int, default=10
            Number of principal components to compute. Should be a positive integer.
            Values below 3 will trigger a warning as they may be insufficient for 
            meaningful analysis.

        Returns
        -------
        None
            Results are written to disk in the results directory with the input_name prefix.

        Raises
        ------
        TypeError
            If pca parameter is not an integer.
        ValueError
            If pca parameter is not positive.

        Notes
        -----
        - If recompute_pca is False, the method will skip PCA computation
        - Uses PLINK's --pca command on the LD-pruned dataset
        - Output files are saved in the results directory specified during initialization
        """

        if not self.recompute_pca:
            logger.info(f"`recompuite_pca` is set to {self.recompute_pca}. PCA will be skipped.")
            logger.info("PCA already performed. Skipping this step.")
            return

        # Check type of pca
        if not isinstance(pca, int):
            raise TypeError("pca should be of type int.")
        if pca <= 0:
            raise ValueError("pca should be a positive integer.")
        if pca <= 3:
            warnings.warn(f"The 'pca' value {pca} is low. Consider increasing it for better results.", UserWarning)

        logger.info("Executing PCA with the following parameters:")
        logger.info(f"PCA parameters: pca={pca}")

        # runs pca analysis
        plink_cmd1 = f"plink --bfile {self.results_dir / (self.input_name+'-LDpruned')} --keep-allele-order --out {self.results_dir / self.input_name} --pca {pca}"

        shell_do(plink_cmd1, log=True)

        return
    
    def generate_plots(self, color_hue_file: Optional[Path] = None, case_control_markers: bool = True, n_neighbors: list = [5], min_dist: list = [0.5], metric: list = ['euclidean'], random_state: Optional[int] = None, umap_kwargs: dict = dict()) -> None:
        """
        Generate UMAP plots with different parameter combinations.
        This method generates UMAP (Uniform Manifold Approximation and Projection) plots using various 
        combinations of parameters. It can incorporate color coding based on metadata and case-control markers.
        
        Parameters
        ----------
        color_hue_file : Path, optional
            Path to a tab-separated file containing color hue information. The file should have at least 3 columns,
            where the first two are ID1 and ID2, and the third column contains the values for color coding.
            Default is None.
        case_control_markers : bool, optional
            Whether to include case-control markers in the plots. If True, reads from the .fam file.
            Default is True. If color_hue_file is not provided, the difference between case and control will be used as hue.
        n_neighbors : list of int, optional
            List of values for the n_neighbors parameter in UMAP. Each value must be positive.
            Default is [5].
        min_dist : list of float, optional
            List of values for the min_dist parameter in UMAP. Each value must be non-negative.
            Default is [0.5].
        metric : list of str, optional
            List of distance metrics to use in UMAP.
            Default is ['euclidean'].
        random_state : int, optional
            Random seed for reproducibility. Must be non-negative.
            Default is None.
        umap_kwargs : dict, optional
            Additional keyword arguments to pass to the UMAP constructor.
            Default is an empty dictionary.
        
        Returns
        -------
        None
            Saves UMAP plots as JPEG files and parameters as a CSV file in the results directory.
        
        Raises
        ------
        TypeError
            If input parameters are not of the correct type.
        ValueError
            If input parameters have invalid values.
        FileNotFoundError
            If color_hue_file is specified but not found.
        
        Notes
        -----
        The method creates a grid of all possible parameter combinations and generates a UMAP plot for each.
        Parameters for each plot are saved in 'plots_parameters.csv'.
        """


        # Check type of n_neighbors
        if not isinstance(n_neighbors, list):
            raise TypeError("n_neighbors should be of type list.")
        if not all(isinstance(i, int) for i in n_neighbors):
            raise TypeError("n_neighbors should be a list of integers.")
        if not all(i > 0 for i in n_neighbors):
            raise ValueError("n_neighbors should be a list of positive integers.")
        if len(n_neighbors) == 0:
            raise ValueError("n_neighbors should not be an empty list.")
        
        # Check type of min_dist
        if not isinstance(min_dist, list):
            raise TypeError("min_dist should be of type list.")
        if not all(isinstance(i, float) for i in min_dist):
            raise TypeError("min_dist should be a list of floats.")
        if not all(i >= 0 for i in min_dist):
            raise ValueError("min_dist should be a list of non-negative floats.")
        if len(min_dist) == 0:
            raise ValueError("min_dist should not be an empty list.")
        
        # Check type of metric
        if not isinstance(metric, list):
            raise TypeError("metric should be of type list.")
        if not all(isinstance(i, str) for i in metric):
            raise TypeError("metric should be a list of strings.")
        if len(metric) == 0:
            raise ValueError("metric should not be an empty list.")
        
        # Check type of random_state
        if random_state is not None:
            if not isinstance(random_state, int):
                raise TypeError("random_state should be of type int.")
            if random_state < 0:
                raise ValueError("random_state should be a non-negative integer.")
            
        # Check if color_hue_file is a file
        if color_hue_file is not None:
            if not isinstance(color_hue_file, Path):
                raise TypeError("color_hue_file should be a Path object.")
            if not color_hue_file.is_file():
                raise FileNotFoundError(f"color_hue_file not found at {color_hue_file}")
        
        # Check if case_control_markers is a boolean
        if not isinstance(case_control_markers, bool):
            raise TypeError("case_control_markers should be of type bool.")

        logger.info("Generating UMAP plots with the following parameters:")
        logger.info(f"UMAP parameters: n_neighbors={n_neighbors}, min_dist={min_dist}, metric={metric}")
        logger.info(f"Random state: {random_state}")
        logger.info(f"Color hue file: {color_hue_file}")
        logger.info(f"Case control markers: {case_control_markers}")

        # generate a parameter grid
        params_dict = {
            'n_neighbors': n_neighbors,
            'min_dist'   : min_dist,
            'metric'     : metric
        }
        param_grid = ParameterGrid(params_dict)

        if color_hue_file is not None:

            if color_hue_file.is_file():
            # load color hue file
                df_color_hue = pd.read_csv(
                    color_hue_file,
                    sep='\t',
                )
                logger.info(f"Color hue file loaded from {color_hue_file}")
                logger.info(f"Column {df_color_hue.columns[2]} will be used for color hue")
                df_color_hue.columns = ["ID1", "ID2", df_color_hue.columns[2]]
                logger.info(f"Color hue file has {df_color_hue.shape[0]} rows and {df_color_hue.shape[1]} columns")
                hue_col = df_color_hue.columns[2]
            else:
                raise FileNotFoundError(f"color_hue_file not found at {color_hue_file}")
        else:
            hue_col = None

        if case_control_markers:
            # load case control markers
            df_fam = pd.read_csv(
                self.input_path / (self.input_name+'.fam'),
                sep=r'\s+',
                engine='python'
            )
            logger.info(f"Case-control labels loaded from {self.input_path / (self.input_name+'.fam')}")
            
            df_fam.columns = ["ID1", "ID2", "F_ID", "M_ID", "Sex", "Phenotype"]
            recode = {1:'Control', 2:'Patient'}
            df_fam["Phenotype"] = df_fam["Phenotype"].map(recode)
            df_fam = df_fam[['ID1', 'ID2', 'Phenotype']].copy()
            logger.info(f"Case-control markers file has {df_fam.shape[0]} rows and {df_fam.shape[1]} columns")

        if color_hue_file is not None and case_control_markers:
            # merge color hue file with case control markers
            df_metadata = df_color_hue.merge(
                df_fam,
                on=['ID1', 'ID2'],
                how='inner'
            )
            logger.info(f"Color hue file merged with case control markers file")
            logger.info(f"Merged file has {df_metadata.shape[0]} rows and {df_metadata.shape[1]} columns")
        elif color_hue_file is not None:
            df_metadata = df_color_hue.copy()
            logger.info(f"Color hue file used as metadata")
        elif case_control_markers:
            df_metadata = df_fam.copy()
            logger.info(f"Case control markers file used as metadata")
        else:
            df_metadata = None
            logger.info(f"No metadata file provided")

        count=1

        # create a dataframe to store parameters
        df_params = pd.DataFrame(
            columns=['n_neighbors', 'min_dist', 'metric', 'warnings']
        )

        for params in param_grid:

            # generate umap plot for data that passed QC
            warnings = self._umap_plots(
                plot_name   =f"umap_2d_{count}.pdf",
                n_neighbors =params['n_neighbors'],
                min_dist    =params['min_dist'],
                metric      =params['metric'],
                random_state=random_state,
                df_metadata =df_metadata,
                hue_col     =hue_col,
                umap_kwargs=umap_kwargs
            )

            self.files_to_keep.append(f"umap_2d_{count}.jpeg")

            df_params.loc[count, 'n_neighbors']= params['n_neighbors']
            df_params.loc[count, 'min_dist']   = params['min_dist']
            df_params.loc[count, 'metric']     = params['metric']
            df_params.loc[count, 'warnings']   = warnings

            count +=1

        # save parameters to a csv file
        df_params.to_csv(
            os.path.join(self.results_dir, 'plots_parameters.csv'),
            index=True,
            sep='\t'
        )

        self.files_to_keep.append('plots_parameters.csv')

        return
    
    def _umap_plots(self, plot_name: str, n_neighbors: int, min_dist: float, metric: str, random_state: Optional[int] = None, df_metadata: Optional[pd.DataFrame] = None, hue_col: Optional[str] = None, umap_kwargs: dict = dict()) -> list:
        """
        Generate and save UMAP (Uniform Manifold Approximation and Projection) plots from PCA data.
        This method reads eigenvector data from a file, performs UMAP dimensionality reduction,
        and creates a 2D scatter plot with optional metadata coloring/styling.
        
        Parameters
        ----------
        plot_name : str
            Name of the output plot file
        n_neighbors : int
            Number of neighbors to consider for manifold approximation
        min_dist : float
            Minimum distance between points in the low dimensional representation
        metric : str
            Distance metric to use for UMAP calculation
        random_state : int, optional
            Random seed for reproducibility
        df_metadata : pd.DataFrame, optional
            DataFrame containing metadata to merge with the eigenvector data
        hue_col : str, optional
            Column name in metadata to use for point coloring
            Additional keyword arguments to pass to UMAP
        umap_kwargs : dict, optional
            Additional keyword arguments to pass to UMAP constructor

        Returns
        -------
        list or None
            List of warning messages if any were generated during execution, None otherwise
        
        Notes
        -----
        The method expects an eigenvector file in the results directory with the naming pattern
        {input_name}.eigenvec. The plot will be saved in the plots directory with the provided
        plot_name.
        If metadata is provided and contains a 'Phenotype' column, it will be used for styling
        points unless hue_col is specified.
        """

        # load .eigenvec file
        df_eigenvec = pd.read_csv(
            self.results_dir / (self.input_name+'.eigenvec'),
            header=None,
            sep=' '
        )
        logger.info(f"Eigenvector file loaded from {self.results_dir / (self.input_name+'.eigenvec')}")
        logger.info(f"Eigenvector file has {df_eigenvec.shape[0]} rows and {df_eigenvec.shape[1]} columns")

        # rename columns
        num_pc = df_eigenvec.shape[1]-2
        new_cols = [f"pca_{k}" for k in range(1,num_pc+1)]
        df_eigenvec.columns = ['ID1', 'ID2'] + new_cols

        df_ids = df_eigenvec[['ID1', 'ID2']].copy()
        df_vals= df_eigenvec[new_cols].to_numpy()

        if df_metadata is not None:
            # merge metadata with eigenvector data
            df_ids = df_ids.merge(
                df_metadata,
                on=['ID1', 'ID2'],
                how='inner'
            )
            logger.info(f"Metadata file merged with eigenvector file")
            logger.info(f"Merged file has {df_ids.shape[0]} rows and {df_ids.shape[1]} columns")

            if 'Phenotype' in df_ids.columns:
                style_col = 'Phenotype'
            else:
                style_col = None
            
            if style_col and hue_col is None:
                hue_col, style_col = style_col, None

        del df_eigenvec

        # instantiate umap class
        D2_redux = umap.UMAP(
            n_components=2,
            n_neighbors =n_neighbors,
            min_dist    =min_dist,
            metric      =metric,
            random_state=random_state,
            **umap_kwargs
        )

        with warnings.catch_warnings(record=True) as w:

            warnings.simplefilter("always")
            
            # generates umap projection
            umap_2D_proj = D2_redux.fit_transform(df_vals)

            df_2D = pd.concat([df_ids, pd.DataFrame(data=umap_2D_proj, columns=['umap1', 'umap2'])], axis=1)

            del df_vals

            # generates and saves a 2D scatter plot
            # size given in inches
            sns.set_context(font_scale=0.9)
            fig, ax = plt.subplots(figsize=(5,5))

            scatter_kwargs = {}
            if style_col is not None:
                scatter_kwargs['style'] = style_col
                
            scatter_plot= sns.scatterplot(
                data=df_2D, 
                x='umap1', 
                y='umap2', 
                hue=hue_col,
                marker='.',
                s=5,
                alpha=0.5,
                ax=ax,
                **scatter_kwargs
            )
            if df_metadata is not None:
                plt.legend(
                    bbox_to_anchor=(0., 1.02, 1., .102), 
                    loc='lower left',
                    ncols=3, 
                    mode="expand", 
                    borderaxespad=0.,
                    fontsize=7,
                    markerscale=2
                )
                
            # Set tick label size
            ax.tick_params(axis='both', labelsize=7)

            # Set axis label and size
            ax.set_xlabel('UMAP1', fontsize=7)
            ax.set_ylabel('UMAP2', fontsize=7)
            ax.set_aspect('equal', adjustable='datalim')

            plt.tight_layout()

            fig.savefig(self.plots_dir / plot_name, dpi=500)
            plt.close(fig)


            warning = [warn.message.args[0] for warn in w] # type: ignore
            return warning

class FstSummary:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, high_ld_file: Path=Path(), built: str = '38', recompute_merge: bool = True, reference_files: dict = dict()) -> None:
        """
        Initialize FstSummary object for Fst analysis.
        
        Parameters
        ----------
        input_path : Path
            Path to the directory containing input files
        input_name : str
            Name of the input file
        output_path : Path
            Path to the directory where results will be saved
        
        Raises
        ------
        TypeError
            If input types are incorrect for any parameter
        FileNotFoundError
            If input_path or output_path do not exist
        """

        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(input_name, str): 
            raise TypeError("input_name should be a string")
        if not input_path.exists():
            raise FileNotFoundError("input_path does not exist")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")
        if not isinstance(built, str):
            raise TypeError("built should be a string")
        if built not in ['37', '38']:
            raise ValueError("built should be either '37' or '38'") 
        if not high_ld_file.is_file():
            logger.info(f"High LD file not found at {high_ld_file}")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions()
            ld_fetcher.get_ld_regions()

            if ld_fetcher.ld_regions is None:
                raise FileNotFoundError("Could not fetch LD regions file.")
                
            high_ld_file = ld_fetcher.ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_file}")

        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.recompute_merge = recompute_merge
        self.high_ld_regions = high_ld_file
        self.built = built

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

        self.results_dir = self.output_path / 'fst_results' 
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.merging_dir = self.results_dir / 'merging'
        self.merging_dir.mkdir(parents=True, exist_ok=True)

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
            output_name= 'cleaned-with-ref',
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

        for file in self.merging_dir.iterdir():
            if file.is_file() and '-merged' not in file.name and file.suffix != '.log':
                file.unlink()

        return
    
    def add_population_tags(self) -> None:
        """
        Add population tags to the merged dataset.
        This method adds population super-population tags from the reference dataset to
        the merged dataset. It reads population information from the reference PSAM file,
        merges it with the study dataset, and assigns 'StPop' (study population) to samples
        not present in the reference dataset.

        Requirements
        ------------
            - Merged dataset files (.bed, .bim, .fam) must exist in the merging directory
            - Reference files dictionary must contain a valid 'psam' Path
        
        Raises
        ------
            FileNotFoundError: If any of the required merged files are not found
            ValueError: If the reference files dictionary doesn't contain a valid 'psam' Path
        
        Side Effects
        ------------
            - Creates a new tab-separated file with population tags at 
              {merging_dir}/cleaned-with-ref-merged-pop-tags.csv
            - Sets self.population_tags to the path of the created file
        
        Returns
        -------
            None
        """
       
        merged_bed = self.merging_dir / 'cleaned-with-ref-merged.bed'
        merged_bim = self.merging_dir / 'cleaned-with-ref-merged.bim'
        merged_fam = self.merging_dir / 'cleaned-with-ref-merged.fam'
        
        if not merged_bed.is_file():
            raise FileNotFoundError(f"Merged data file not found at {merged_bed}")
        if not merged_bim.is_file():
            raise FileNotFoundError(f"Merged BIM file not found at {merged_bim}")
        if not merged_fam.is_file():
            raise FileNotFoundError(f"Merged FAM file not found at {merged_fam}")

        if 'psam' not in self.reference_files or not isinstance(self.reference_files['psam'], Path):
            raise ValueError("Reference files dictionary must contain a valid 'psam' Path")

        reference_tags = self.reference_files['psam']
        
        df_tags = pd.read_csv(reference_tags, sep=r"\s+", engine='python')
        df_tags['ID'] = '0'
        df_tags = df_tags[['ID', '#IID', 'SuperPop']]
        df_tags = df_tags.rename(columns={'ID': 'ID1', '#IID': 'ID2', 'SuperPop': 'SuperPop'})

        logger.info(f"Population tags loaded from {reference_tags}")
        logger.info(f'Population tags columns: {df_tags.columns.tolist()}')

        df_merged_fam = pd.read_csv(merged_fam, sep=r"\s+", header=None, engine='python')
        df_merged_fam = df_merged_fam.rename(columns={0: 'ID1', 1: 'ID2'})
        #df_merged_fam['SuperPop'] = 'StPop'

        logger.info(f"Merged BIM file loaded from {merged_fam}")
        logger.info(f'Merged BIM file columns: {df_merged_fam.columns.tolist()}')

        df = pd.merge(
            df_merged_fam[['ID1', 'ID2']],
            df_tags,
            on=['ID1', 'ID2'],
            how='left'
        )
        df['SuperPop'] = df['SuperPop'].fillna('StPop')
        logger.info(f"Added population tags to the merged dataset")

        self.population_tags = self.merging_dir / 'cleaned-with-ref-merged-pop-tags.csv'
        df.to_csv(
            self.population_tags,
            index=False,
            sep='\t'
        )

        return
    
    def compute_fst(self) -> None:
        """
        Compute FST (fixation index) statistics between populations.

        This method calculates FST statistics between each super-population in the dataset
        and a study population ('StPop'). The process involves:
        1. Reading population tags from the specified file
        2. For each unique super-population (except 'StPop'):
        - Creating population filter files (keep and within files)
        - Running PLINK commands to filter the dataset and compute FST statistics

        The method requires the following instance variables to be set:
            - population_tags: Path to a file containing population information
            - results_dir: Directory where results will be stored
            - merging_dir: Directory containing the merged genotype data

        Returns:
        --------
        None
        """

        df_tags = pd.read_csv(self.population_tags, sep=r"\s+", engine='python')

        files = dict()

        for pop in df_tags['SuperPop'].unique():
            if pop != 'StPop':
                df_temp = df_tags[(df_tags['SuperPop'] == pop) | (df_tags['SuperPop'] == 'StPop')].reset_index(drop=True)
                df_temp[['ID1', 'ID2']].to_csv(self.results_dir / f'keep-{pop}_StPop.txt', sep='\t', index=False, header=False)
                df_temp.to_csv(self.results_dir / f'within-{pop}_StPop.txt', index=False, header=False, sep='\t',)

                files[pop] = (self.results_dir / f'keep-{pop}_StPop.txt', self.results_dir / f'within-{pop}_StPop.txt')

                logger.info(f"Created keep and within files for population {pop}")

        input_file = self.merging_dir / 'cleaned-with-ref-merged'

        for key in files.keys():

            keep_file, within_file = files[key]
            output_file = self.results_dir / f'keep-{key}-StPop'

            plink_cmd1 = f"plink --bfile {input_file} --keep {keep_file} --make-bed --out {output_file}"
            plink_cmd2 = f"plink --bfile {output_file} --fst --within {within_file} --out {self.results_dir / f'fst-{key}-StPop'}"

            plink_cmds = [plink_cmd1, plink_cmd2]
            for cmd in plink_cmds:
                shell_do(cmd, log=True)
        logger.info("Fst computation completed for all populations.")

        return
    
    def report_fst(self) -> pd.DataFrame:
        """
        Generate a report of Fst results.
        This method reads the Fst results from the results directory and generates a summary report.
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing the Fst results summary
        
        Raises
        ------
        FileNotFoundError
            If no Fst result files are found in the results directory.
        """

        df_summary = pd.DataFrame(columns=['SuperPop', 'Fst', 'WeightedFst'])

        # Get a list of all log files in the results directory
        log_files = [f for f in self.results_dir.iterdir() if f.is_file() and f.suffix == '.log']

        if not log_files:
            raise FileNotFoundError(f"No log files found in {self.results_dir}")

        logger.info(f"Found {len(log_files)} log files in {self.results_dir}")

        # Extract the population names from log file names
        # Assuming log files follow the pattern 'fst-{population}-StPop.log'
        files = {}
        for log_file in log_files:
            if log_file.stem.startswith('fst-') and log_file.stem.endswith('-StPop'):
                pop = log_file.stem.split
                files[pop] = log_file

        if not files:
            raise FileNotFoundError(f"No Fst result files found in {self.results_dir}")

        for key in files.keys():
        
            log_file = files[key]
            with open(log_file, 'r') as f:

                lines = f.readlines()
                for line in lines:
                    if line.startswith('Mean Fst'):
                        fst = line.split(':')[1].strip()
                    if line.startswith('Weighted Fst'):
                        weighted_fst = line.split(':')[1].strip()
                df_summary = pd.concat([df_summary, pd.DataFrame({'SuperPop': [key], 'Fst': [fst], 'WeightedFst': [weighted_fst]})], ignore_index=True)

        df_summary.to_csv(
            self.results_dir / 'fst_summary.csv',
            index=False,
            sep='\t'
        )
        logger.info(f"Fst summary report generated at {self.results_dir / 'fst_summary.csv'}")

        for file in self.results_dir.iterdir():
            if file.is_file() and (file.suffix == '.bed' or file.suffix == '.bim' or file.suffix == '.fam'):
                file.unlink()
        
        return df_summary