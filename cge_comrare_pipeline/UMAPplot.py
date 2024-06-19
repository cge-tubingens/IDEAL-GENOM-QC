"""
Module to draw plots based on UMAP dimension reduction
"""

import os
import subprocess
import shutil
import umap

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from cge_comrare_pipeline.Helpers import shell_do, delete_temp_files
from sklearn.model_selection import ParameterGrid

class UMAPplot:

    def __init__(self, sampleQC_path:str, sampleQC_name:str, pcaQC_path:str, pcaQC_name:str, dependables_path:str, config_dict:dict, output_path:str) -> None:

        """
        
        """

        # check if paths are set
        if sampleQC_path is None or pcaQC_path is None or dependables_path is None or output_path is None:
            raise ValueError("values for sampleQC_path, pcaQC_path, dependables_path and output_path must be set upon initialization.")

        # Check path validity of past sample QC data
        bed_path = os.path.join(sampleQC_path, sampleQC_name + '.bed')
        fam_path = os.path.join(sampleQC_path, sampleQC_name + '.fam')
        bim_path = os.path.join(sampleQC_path, sampleQC_name + '.bim')

        bed_check = os.path.exists(bed_path)
        fam_check = os.path.exists(fam_path)
        bim_check = os.path.exists(bim_path)

        if not bed_check:
            raise FileNotFoundError(".bed file not found")
        if not fam_check:
            raise FileNotFoundError(".fam file not found")
        if not bim_check:
            raise FileNotFoundError(".bim file not found")
        
        # Check path validity of past sample QC data
        bed_path = os.path.join(pcaQC_path, pcaQC_name + '.bed')
        fam_path = os.path.join(pcaQC_path, pcaQC_name + '.fam')
        bim_path = os.path.join(pcaQC_path, pcaQC_name + '.bim')

        bed_check = os.path.exists(bed_path)
        fam_check = os.path.exists(fam_path)
        bim_check = os.path.exists(bim_path)

        if not bed_check:
            raise FileNotFoundError(".bed file not found")
        if not fam_check:
            raise FileNotFoundError(".fam file not found")
        if not bim_check:
            raise FileNotFoundError(".bim file not found")

        if not os.path.exists(dependables_path):
            raise FileNotFoundError("dependables_path is not a valid path")
        if not os.path.exists(output_path):
            raise FileNotFoundError("output_path is not a valid path")

        self.sampleQC_path= sampleQC_path
        self.pcaQC_path   = pcaQC_path
        self.sampleQC_name= sampleQC_name
        self.pcaQC_name   = pcaQC_name
        self.dependables  = dependables_path
        self.config_dict  = config_dict

        self.files_to_keep= []

        # create results folder
        self.results_dir = os.path.join(output_path, 'umap_plots')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

    def compute_pcas(self)->None:

        sampleQC_path= self.sampleQC_path
        pcaQC_path   = self.pcaQC_path
        sampleQC_name= self.sampleQC_name
        pcaQC_name   = self.pcaQC_name
        results_dir  = self.results_dir

        pca = self.config_dict['pca']

        step= "compute_pca_for_umap_plots"

        # runs pca analysis
        plink_cmd1 = f"plink --bfile {os.path.join(sampleQC_path, sampleQC_name)} --keep-allele-order --maf 0.01 --out {os.path.join(results_dir, '.sampleQC.pca')} --pca {pca}"

        plink_cmd2 = f"plink --bfile {os.path.join(pcaQC_path, pcaQC_name)} --keep-allele-order --maf 0.01 --out {os.path.join(results_dir, '.nooutliers.pca')} --pca {pca}"

        # execute plink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
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

        sampleQC_path= self.sampleQC_path
        pcaQC_path   = self.pcaQC_path
        sampleQC_name= self.sampleQC_name
        pcaQC_name   = self.pcaQC_name
        results_dir  = self.results_dir
        dependables  = self.dependables

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

        df_plots = pd.DataFrame(columns=['n_neighbors', 'min_dist', 'metric'])

        for params in param_grid:

            # generate umap plot for data that passed Sample QC
            df_samp = self.umap_plots(
                path_to_data=os.path.join(results_dir, '.sampleQC.pca.eigenvec'),
                output_file =os.path.join(results_dir, f"sampleQC_umap_2d_{count}.pdf"),
                geo_path=geo_info_path,
                fam_path=os.path.join(sampleQC_path, sampleQC_name+".fam"),
                n_neighbors =params['n_neighbors'],
                min_dist    =params['min_dist'],
                metric      =params['metric']
            )

            # generate umap plot for data that passed Sample QC and Ethnicity check
            df_out = self.umap_plots(
                path_to_data=os.path.join(results_dir, '.nooutliers.pca.eigenvec'),
                output_file =os.path.join(results_dir, f"nooutliers_umap_2d_{count}.pdf"),
                geo_path=geo_info_path,
                fam_path=os.path.join(pcaQC_path, pcaQC_name+".fam"),
                n_neighbors =params['n_neighbors'],
                min_dist    =params['min_dist'],
                metric      =params['metric']
            )

            df_plots.loc[count, 'n_neighbors'] = params['n_neighbors']
            df_plots.loc[count, 'min_dist'] = params['min_dist']
            df_plots.loc[count, 'metric'] = params['metric']

            count +=1

        df_plots.to_csv(
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
    def umap_plots(path_to_data:str, output_file:str, geo_path:str, fam_path:str, n_neighbors:int, min_dist:float, metric:str)->pd.DataFrame:

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
            df_2D = pd.concat([df_ids, pd.DataFrame(data=umap_2D_proj, columns=['umap_1', 'umap_2'])], axis=1)
            df_2D = pd.merge(
                df_2D,
                df_geo,
                left_on='ID2',
                right_on=df_geo.columns[0]
            ).drop(columns=[df_geo.columns[0]])

            # generates and saves a 2D scatter plot
            fig, ax = plt.subplots(figsize=(10,10))
            scatter_plot= sns.scatterplot(
                data=df_2D, 
                x='umap_1', 
                y='umap_2', 
                hue=df_geo.columns[1],
                marker='.',
                alpha=0.6,
                ax=ax,
                style="Phenotype"
            )
            scatter_fig = scatter_plot.get_figure()
            scatter_fig.savefig(output_file)
            plt.close()
        else:
            # prepares data for plotting
            df_2D = pd.concat([df_ids, pd.DataFrame(data=umap_2D_proj, columns=['umap_1', 'umap_2'])], axis=1)

            # generates and saves a 2D scatter plot
            fig, ax = plt.subplots(figsize=(10,10))
            scatter_plot= sns.scatterplot(
                data=df_2D, 
                x='umap_1',
                y='umap_2',
                marker='.',
                alpha=0.6,
                ax=ax,
                style="Phenotype"
            )
            scatter_fig = scatter_plot.get_figure()
            scatter_fig.savefig(output_file)
            plt.close()

        return df_2D

