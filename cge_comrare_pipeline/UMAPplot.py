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

        self.sampleQC_path = sampleQC_path
        self.pcaQC_path= pcaQC_path
        self.sampleQC_name = sampleQC_name
        self.pcaQC_name= pcaQC_name
        self.dependables= dependables_path
        self.config_dict= config_dict

        self.files_to_keep = []
        
        # create figures folder
        self.plots_dir = os.path.join(output_path, 'plots')
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)
