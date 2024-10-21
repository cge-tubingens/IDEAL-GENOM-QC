"""
Module to draw plots based on UMAP dimension reduction
"""

import os
import umap

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from cge_comrare_pipeline.Helpers import shell_do, delete_temp_files
from sklearn.model_selection import ParameterGrid

class UMAPplot:

    def __init__(self, input_path:str, input_name:str, dependables_path:str, config_dict:dict, output_path:str) -> None:

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
        self.config_dict= config_dict

        self.files_to_keep= []

        # create results folder
        self.results_dir = os.path.join(output_path, 'umap_plots')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        pass

    def ld_pruning(self)->dict:

        """
        Prune samples based on Linkage Disequilibrium (LD).

        This method performs LD-based sample pruning using PLINK commands. It filters samples based on Minor Allele Frequency (maf), genotype missingness (geno), individual missingness (mind), and Hardy-Weinberg Equilibrium (hwe). Additionally, it excludes SNPs located in high LD regions specified in the dependables path. The resulting pruned dataset is saved as a new binary file.

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

        maf      = self.config_dict['maf']
        geno     = self.config_dict['geno']
        mind     = self.config_dict['mind']
        hwe      = self.config_dict['hwe']
        ind_pair = self.config_dict['indep-pairwise']

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
        #if maf < 0.05 or maf > 0.1:
        #    raise ValueError("maf should be between 0.05 and 0.1")
        
        # Check if geno is in range
        if geno < 0.05 or geno > 0.1:
            raise ValueError("geno should be between 0.05 and 0.1")
        
        # Check if mind is in range
        #if mind < 0.1 or mind > 0.15:
        #    raise ValueError("mind should be between 0.1 and 0.15")
        
        # Check if hwe is in range
        if hwe < 0.00000001 or hwe > 0.001:
            raise ValueError("hwe should be between 0.00000001 and 0.001")
        
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
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --extract {os.path.join(results_dir, input_name+'.prune.in')} --make-bed --threads {max_threads} --out {os.path.join(results_dir, input_name+'.pruned')}"

        # execute PLINK commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

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

    def compute_pcas(self)->None:

        input_name = self.input_name
        results_dir= self.results_dir

        pca = self.config_dict['umap_pca']

        step= "compute_pca_for_umap_plots"

        # runs pca analysis
        plink_cmd1 = f"plink --bfile {os.path.join(results_dir, input_name+'.pruned')} --keep-allele-order --maf 0.01 --out {os.path.join(results_dir, 'cleaned_samples.pca')} --pca {pca}"

        # execute plink command
        cmd = plink_cmd1
        shell_do(cmd, log=True)

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
    
    def generate_plots(self)->None:

        input_path = self.input_path
        input_name = self.input_name
        results_dir= self.results_dir
        dependables= self.dependables

        n_neighbors = self.config_dict['umap_n_neighbors']
        min_dist    = self.config_dict['umap_min_dist']
        metric      = self.config_dict['umap_metric']

        step = "draw_umap_plots"

        params_dict = {
            'n_neighbors': n_neighbors,
            'min_dist'   : min_dist,
            'metric'     : metric
        }
        param_grid = ParameterGrid(params_dict)

        count=1

        geo_info_path = os.path.join(dependables, 'geographic_info.txt')

        df_params = pd.DataFrame(columns=['n_neighbors', 'min_dist', 'metric', 'warnings'])

        for params in param_grid:

            # generate umap plot for data that passed Sample QC
            warnings = self.umap_plots(
                path_to_data=os.path.join(results_dir, 'cleaned_samples.pca.eigenvec'),
                output_file =os.path.join(results_dir, f"umap_2d_{count}.png"),
                geo_path=geo_info_path,
                fam_path=os.path.join(input_path, input_name+".fam"),
                n_neighbors =params['n_neighbors'],
                min_dist    =params['min_dist'],
                metric      =params['metric'],
                fig_num=count
            )

            df_params.loc[count, 'n_neighbors']= params['n_neighbors']
            df_params.loc[count, 'min_dist']   = params['min_dist']
            df_params.loc[count, 'metric']     = params['metric']
            df_params.loc[count, 'warnings']   = warnings

            count +=1

        df_params.to_csv(
            os.path.join(results_dir, 'plots_parameters.csv'),
            index=True
        )

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
    def umap_plots(path_to_data:str, output_file:str, geo_path:str, fam_path:str, n_neighbors:int, min_dist:float, metric:str, fig_num:int=1):

        """
        Generates a 2D UMAP projection plot with geographic information and saves it to a file.

        Parameters:
        -----------
        path_to_data (str): 
            Path to the .eigenvec file containing the PCA data.
        output_file (str): 
            Path to the output file where the generated plot will be saved.
        geo_path (str): 
            Path to the file containing geographic information.
        n_neighbors (int): 
            The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
        min_dist (float): 
            The minimum distance between points in the low-dimensional space.
        metric (str): 
            The metric to use for the UMAP algorithm.

        Returns:
        --------
        pandas.DataFrame
        """

        import warnings

        # load .eigenvec file
        df_eigenvec = pd.read_csv(
            path_to_data,
            header=None,
            sep=' '
        )

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

        df_ids = df_ids.merge(df_fam[["ID1", "ID2", "Phenotype"]], on=['ID1', 'ID2'])

        del df_eigenvec, df_fam

        # instantiate umap class
        D2_redux = umap.UMAP(
            n_components=2,
            n_neighbors =n_neighbors,
            min_dist    =min_dist,
            metric      =metric
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
                sns.set_context(font_scale=0.9)
                fig, ax = plt.subplots(figsize=(10,10))
                scatter_plot= sns.scatterplot(
                    data=df_2D, 
                    x='umap1', 
                    y='umap2', 
                    hue=df_geo.columns[1],
                    marker='.',
                    alpha=0.6,
                    ax=ax,
                    style="Phenotype"
                )
                plt.legend(fontsize='10', markerscale=2)
                
                # caption = f"Figure {fig_num}: min_dis={min_dist}, n_neighbors={n_neighbors}, metric={metric}."
                # plt.figtext(0.5, 0.05, wrap=True, horizontalalignment='center', fontsize=12)

                scatter_fig = scatter_plot.get_figure()
                scatter_fig.savefig(output_file)
                plt.close()
            else:
                # prepares data for plotting
                df_2D = pd.concat([df_ids, pd.DataFrame(data=umap_2D_proj, columns=['umap1', 'umap2'])], axis=1)

                # generates and saves a 2D scatter plot
                fig, ax = plt.subplots(figsize=(10,10))
                scatter_plot= sns.scatterplot(
                    data=df_2D, 
                    x='umap1',
                    y='umap2',
                    marker='.',
                    alpha=0.6,
                    ax=ax,
                    style="Phenotype"
                )
                
                # caption = f"Figure {fig_num}: min_dis={min_dist}, n_neighbors={n_neighbors}."
                # plt.figtext(0.5, -0.05, wrap=True, horizontalalignment='center', fontsize=12)

                scatter_fig = scatter_plot.get_figure()
                scatter_fig.savefig(output_file)
                plt.close()

            if isinstance(w, list):
                warning = [warn.message.args[0] for warn in w]
                return warning
            else:
                return None


