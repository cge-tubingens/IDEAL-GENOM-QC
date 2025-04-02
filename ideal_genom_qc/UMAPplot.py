"""
Module to draw plots based on UMAP dimension reduction
"""

import os
import umap
import warnings
import logging
import psutil

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path

from ideal_genom_qc.Helpers import shell_do, delete_temp_files
from ideal_genom_qc.get_references import FetcherLDRegions
from sklearn.model_selection import ParameterGrid

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class UMAPplot:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name: str, high_ld_file: Path, built: str = '38', recompute_pca: bool = True) -> None:



        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(high_ld_file, Path):
            raise TypeError("high_ld_regions should be a Path object")
        if not isinstance(input_name, str): 
            raise TypeError("input_name should be a string")
        if not isinstance(output_name, str):
            raise TypeError("output_name should be a string")
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

            high_ld_file = ld_fetcher.ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_file}")

        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.output_name= output_name
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
        Compute Principal Component Analysis (PCA) to feed UMAP algorithm and generate plots.

        This method performs PCA analysis using PLINK1.9 on the input dataset specified by `input_name` and stores the results in the `results_dir`. The number of principal components to compute is defined in the configuration dictionary under the key 'umap_pca'. If `compute_all` is set to True, the PCA analysis is executed; otherwise, it assumes the principal components have already been computed.

        Raises:
        -------
            TypeError: If the 'umap_pca' value in the configuration dictionary is not an integer.

        Returns:
        --------
            dict: A dictionary containing the status of the process, the step name, and the output directory for the plots.
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
    
    def generate_plots(self, color_hue_file: Path = None, case_control_markers: bool = True, n_neighbors: list = [5], min_dist: list = [0.5], metric: list = ['euclidean'], random_state: int = None) -> None:


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
            # load color hue file
            df_color_hue = pd.read_csv(
                color_hue_file,
                sep=r'\s+',
                engine='python'
            )
            logger.info(f"Color hue file loaded from {color_hue_file}")
            logger.info(f"Column {df_color_hue.columns[2]} will be used for color hue")
            df_color_hue.columns = ["ID1", "ID2", df_color_hue.columns[2]]
            logger.info(f"Color hue file has {df_color_hue.shape[0]} rows and {df_color_hue.shape[1]} columns")
            hue_col = df_color_hue.columns[2]
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
                plot_name   =f"umap_2d_{count}.jpeg",
                n_neighbors =params['n_neighbors'],
                min_dist    =params['min_dist'],
                metric      =params['metric'],
                random_state=random_state,
                df_metadata =df_metadata,
                hue_col     =hue_col,
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
    
    def _umap_plots(self, plot_name: str, n_neighbors: int, min_dist: float, metric: str, random_state: int = None, df_metadata: pd.DataFrame = None, hue_col: str = None):
        
        """
        Generates UMAP plots from PCA eigenvector data and saves the plot to a .jpeg file.

        Parameters:
        -----------
        path_to_data : str
            Path to the .eigenvec file containing PCA eigenvector data.
        output_file : str
            Path to the output file where the UMAP plot will be saved.
        geo_path : str
            Path to the file containing geographic information. If the file does not exist, the plot will be generated without geographic information acting as hue.
        fam_path : str
            Path to the .fam (PLINK1.9 file format) file containing family information.
        n_neighbors : int
            The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
        min_dist : float
            The effective minimum distance between embedded points used for manifold approximation.
        metric : str
            The metric to use for distance computation during manifold approximation.

        Returns:
        --------
        list or None
            A list of warning messages if any warnings were raised during the UMAP projection, otherwise None.
        
        Notes:
        ------
        The function reads PCA eigenvector data (from .eigenvec file) and family information (from .fam file), performs UMAP dimensionality reduction, and generates a 2D scatter plot. To ensure reproducibility of the plot, the random state is set to 42.

        If geographic information is provided, it will be included in the plot. The plot is saved to the specified output file in a .jpeg file.
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
            random_state=random_state
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

            scatter_plot= sns.scatterplot(
                data=df_2D, 
                x='umap1', 
                y='umap2', 
                hue=hue_col,
                marker='.',
                s=5,
                alpha=0.6,
                ax=ax,
                style=style_col,
                edgecolor='none'
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

            plt.tight_layout()

            scatter_fig = scatter_plot.get_figure()
            scatter_fig.savefig(self.plots_dir / plot_name, dpi=400)
            plt.close()


            if isinstance(w, list):
                warning = [warn.message.args[0] for warn in w]
                return warning
            else:
                return None
