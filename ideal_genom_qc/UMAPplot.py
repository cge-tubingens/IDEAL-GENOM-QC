"""
Module to draw plots based on UMAP dimension reduction
"""

import os
import umap
import warnings

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ideal_genom_qc.Helpers import shell_do, delete_temp_files
from sklearn.model_selection import ParameterGrid

class UMAPplot:

    def __init__(self, input_path:str, input_name:str, dependables_path:str, output_path:str, compute_all:bool=True) -> None:

        """
        
        """

        # check if paths are set
        if input_path is None or dependables_path is None or output_path is None:
            raise ValueError("values for sampleQC_path, pcaQC_path, dependables_path and output_path must be set upon initialization.")

        # Check path validity of clean data
        bed_path = os.path.join(input_path, input_name + '.bed')
        fam_path = os.path.join(input_path, input_name + '.fam')
        bim_path = os.path.join(input_path, input_name + '.bim')

        bed_check = os.path.exists(bed_path)
        fam_check = os.path.exists(fam_path)
        bim_check = os.path.exists(bim_path)

        if not bed_check:
            raise FileNotFoundError(".bed file not found")
        if not fam_check:
            raise FileNotFoundError(".fam file not found")
        if not bim_check:
            raise FileNotFoundError(".bim file not found")
        
        # Check path validity of dependables
        if not os.path.exists(dependables_path):
            raise FileNotFoundError("dependables_path is not a valid path")
        # Check path validity of output_path
        if not os.path.exists(output_path):
            raise FileNotFoundError("output_path is not a valid path")

        self.input_path = input_path
        self.input_name = input_name
        self.dependables= dependables_path
        self.compute_all= compute_all

        self.files_to_keep= []

        # create results folder
        self.results_dir = os.path.join(output_path, 'umap_plots')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        pass

    def ld_pruning(self, maf:float=0.01, geno:float=0.1, mind:float=0.2, hwe:float=5e-8, ind_pair:list=[50, 5, 0.2])->dict:

        """
        Prune samples based on Linkage Disequilibrium (LD).

        This method performs LD-based sample pruning using PLINK1.9 commands. It filters samples based on Minor Allele Frequency (maf), genotype missingness (geno), individual missingness (mind), and Hardy-Weinberg Equilibrium (hwe). Additionally, it excludes SNPs located in high LD regions specified in the dependables path. The resulting pruned dataset is saved as a new binary file.

        Raises:
        -------
        - TypeError: If maf, geno, mind, or hwe is not of type float.
        - ValueError: If maf, geno, mind, or hwe is not within the specified range.
        - FileNotFoundError: If the file with high LD regions is not found.

        Returns:
        --------
        - dict: A dictionary containing information about the process completion status, the step performed, and the output files generated.
        """

        input_path      = self.input_path
        input_name      = self.input_name
        dependables_path= self.dependables
        results_dir     = self.results_dir
        compute_all     = self.compute_all


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
        if maf < 0.0 or maf > 0.5:
            raise ValueError("maf should be between 0 and 0.5")
        
        # Check if geno is in range
        if geno < 0.05 or geno > 0.1:
            raise ValueError("geno should be between 0.05 and 0.1")
        
        # Check if mind is in range
        if mind < 0 or mind > 1:
            raise ValueError("mind should be between 0 and 1")
        
        # Check if mind is around typical values
        if mind <= 0.02 and mind >= 0.1:
            warnings.warn(f"The 'mind' value {mind} is outside the recommended range of 0.02 to 0.1.", UserWarning)

        # Check if hwe is in range
        if hwe < 0 or hwe > 1:
            raise ValueError("hwe should be between 0 and 1")
        
        # check existence of high LD regions file
        high_ld_regions_file = os.path.join(dependables_path, 'high-LD-regions.txt')
        if not os.path.exists(high_ld_regions_file):
            raise FileNotFoundError(f"File with high LD region was not found: {high_ld_regions_file}")

        step = "ld_prune"

        if os.cpu_count() is not None:
            max_threads = os.cpu_count()-2
        else:
            max_threads = 10

        # generates prune.in and prune.out files
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --maf {maf} --geno {geno} --mind {mind} --hwe {hwe} --exclude {high_ld_regions_file} --range --indep-pairwise {ind_pair[0]} {ind_pair[1]} {ind_pair[2]} --threads {max_threads} --out {os.path.join(results_dir, input_name)}"

        # prune and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --extract {os.path.join(results_dir, input_name+'.prune.in')} --make-bed --threads {max_threads} --out {os.path.join(results_dir, input_name+'-LDpruned')}"

        self.files_to_keep.append(input_name+'-LDpruned.bed')
        self.files_to_keep.append(input_name+'-LDpruned.bim')
        self.files_to_keep.append(input_name+'-LDpruned.fam')

        if compute_all:
            # execute PLINK commands
            cmds = [plink_cmd1, plink_cmd2]
            for cmd in cmds:
                shell_do(cmd, log=True)
        else:
            print(f"\033[1m LD prunning already performed.\033[0m")

        # report
        process_complete = True

        outfiles_dict = {
            'plink_out': results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def compute_pcas(self, pca:int=10)->None:
    
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

        input_name = self.input_name
        results_dir= self.results_dir
        compute_all= self.compute_all

        # Check type of pca
        if not isinstance(pca, int):
            raise TypeError("pca should be of type int.")

        step= "compute_pca_for_umap_plots"

        # runs pca analysis
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'-LDpruned')} --keep-allele-order --maf 0.01 --out {os.path.join(results_dir, input_name)} --pca {pca}"

        if compute_all:
            # execute plink command
            cmd = plink_cmd1
            shell_do(cmd, log=True)
        else:
            print(f"\033[1m Principal components already computed.\033[0m")

        self.files_to_keep.append(input_name+'.eigenvec')

        # report
        process_complete = True

        outfiles_dict = {
            'plots_out': self.results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def generate_plots(self, n_neighbors:list=[5], min_dist:list=[0.5], metric:list=['euclidean'])->None:

        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir
        dependables= self.dependables

        step = "draw_umap_plots"

        # generate a parameter grid
        params_dict = {
            'n_neighbors': n_neighbors,
            'min_dist'   : min_dist,
            'metric'     : metric
        }
        param_grid = ParameterGrid(params_dict)

        count=1

        # path to geographic information
        geo_info_path = os.path.join(dependables, 'geographic_info.txt')

        # create a dataframe to store parameters
        df_params = pd.DataFrame(
            columns=['n_neighbors', 'min_dist', 'metric', 'warnings']
        )

        for params in param_grid:

            # generate umap plot for data that passed QC
            warnings = self.umap_plots(
                path_to_data=os.path.join(results_dir, input_name+'.eigenvec'),
                output_file =os.path.join(results_dir, f"umap_2d_{count}.jpeg"),
                geo_path    =geo_info_path,
                fam_path    =os.path.join(input_path, input_name+".fam"),
                n_neighbors =params['n_neighbors'],
                min_dist    =params['min_dist'],
                metric      =params['metric'],
            )

            self.files_to_keep.append(f"umap_2d_{count}.jpeg")

            df_params.loc[count, 'n_neighbors']= params['n_neighbors']
            df_params.loc[count, 'min_dist']   = params['min_dist']
            df_params.loc[count, 'metric']     = params['metric']
            df_params.loc[count, 'warnings']   = warnings

            count +=1

        # save parameters to a csv file
        df_params.to_csv(
            os.path.join(results_dir, 'plots_parameters.csv'),
            index=True,
            sep='\t'
        )

        self.files_to_keep.append('plots_parameters.csv')

        # delete temporary files
        delete_temp_files(self.files_to_keep, results_dir)

        # report
        process_complete = True

        outfiles_dict = {
            'plots_out': self.results_dir
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    @staticmethod
    def umap_plots(path_to_data:str, output_file:str, geo_path:str, fam_path:str, n_neighbors:int, min_dist:float, metric:str):
        
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

        import warnings

        # load .eigenvec file
        df_eigenvec = pd.read_csv(
            path_to_data,
            header=None,
            sep=' '
        )

        # load .fam file
        df_fam = pd.read_csv(
            fam_path,
            header=None,
            sep=' '
        )
        df_fam.columns = ["ID1", "ID2", "F_ID", "M_ID", "Sex", "Phenotype"]
        recode = {1:'control', 2:'case'}
        df_fam["Phenotype"] = df_fam["Phenotype"].map(recode)

        # rename columns
        num_pc = df_eigenvec.shape[1]-2
        new_cols = [f"pca_{k}" for k in range(1,num_pc+1)]
        df_eigenvec.columns = ['ID1', 'ID2'] + new_cols

        df_ids = df_eigenvec[['ID1', 'ID2']].copy()
        df_vals= df_eigenvec[new_cols].to_numpy()

        df_ids = df_ids.merge(
            df_fam[["ID1", "ID2", "Phenotype"]], on=['ID1', 'ID2']
        )

        del df_eigenvec, df_fam

        # instantiate umap class
        D2_redux = umap.UMAP(
            n_components=2,
            n_neighbors =n_neighbors,
            min_dist    =min_dist,
            metric      =metric,
            random_state=42
        )

        with warnings.catch_warnings(record=True) as w:

            warnings.simplefilter("always")
            
            # generates umap projection
            umap_2D_proj = D2_redux.fit_transform(df_vals)

            del df_vals

            if os.path.isfile(geo_path):

                # load file with geographic info
                df_geo = pd.read_csv(
                    geo_path,
                    sep=' ',
                    index_col=False 
                )

                # prepares data for plotting
                df_2D = pd.concat([df_ids, pd.DataFrame(data=umap_2D_proj, columns=['umap1', 'umap2'])], axis=1)
                df_2D = pd.merge(
                    df_2D,
                    df_geo,
                    left_on='ID2',
                    right_on=df_geo.columns[0]
                ).drop(columns=[df_geo.columns[0]])

                # generates and saves a 2D scatter plot
                # size given in inches
                sns.set_context(font_scale=0.9)
                fig, ax = plt.subplots(figsize=(5,5))
                scatter_plot= sns.scatterplot(
                    data=df_2D, 
                    x='umap1', 
                    y='umap2', 
                    hue=df_geo.columns[1],
                    marker='.',
                    s=10,
                    alpha=0.6,
                    ax=ax,
                    style="Phenotype",
                    edgecolor='none'
                )
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
                scatter_fig.savefig(output_file, dpi=600)
                plt.close()
            else:
                # prepares data for plotting
                df_2D = pd.concat(
                    [df_ids, pd.DataFrame(data=umap_2D_proj, columns=['umap1', 'umap2'])], axis=1
                )

                # generates and saves a 2D scatter plot
                # size given in inches
                fig, ax = plt.subplots(figsize=(5,5))
                scatter_plot= sns.scatterplot(
                    data=df_2D, 
                    x='umap1',
                    y='umap2',
                    marker='.',
                    alpha=0.6,
                    s=10,
                    ax=ax,
                    style="Phenotype",
                    edgecolor='none'
                )

                # Set tick label size
                ax.tick_params(axis='both', labelsize=7)

                # Set axis label and size
                ax.set_xlabel('UMAP1', fontsize=7)
                ax.set_ylabel('UMAP2', fontsize=7)

                scatter_fig = scatter_plot.get_figure()
                plt.tight_layout()
                scatter_fig.savefig(output_file, dpi=600)
                plt.close()

            if isinstance(w, list):
                warning = [warn.message.args[0] for warn in w]
                return warning
            else:
                return None


