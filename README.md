# Genotype Quality Control Pipeline

[![Documentation Status](https://readthedocs.org/projects/ideal-genom-qc/badge/?version=latest)](https://ideal-genom-qc.readthedocs.io/en/latest/)
[![PyPI version](https://img.shields.io/pypi/v/ideal-genom-qc.svg)](https://pypi.org/project/ideal-genom-qc/)


This Python package is designed to execute a genotype quality control pipeline, encapsulating several years of research at CGE Tübingen.

## Basic requirements

The quality control pipeline is built on `PLINK` as the main tool. The `ideal_genom_qc` serves as a wrapper for the various QC pipeline steps. To run the pipeline, `PLINK1.9` and `PLINK2` must be installed on the system.

The pipeline is designed to seamlessly run with minimal input and produce cleaned binary files as a result as well as several plots along the way. To accomplish this, the following folder structure is expected:

```
projectFolder
    |
    |---inputData
    |
    |---outputData
    |
    |---configFiles
    |
    |---dependables
```
1. The `inputData` folder should contain the binary files with the genotype data to be analyzed in `PLINK` format (`.bed`, `.bim`, `.fam` files).

2. The `outputData` folder will contain the resultant files of the quality control pipeline. Below, the pipeline output will be detailed.

3. The `dependables` folder is designed to contain complemenatry files for the quality control pipeline. This folder is optional.

4. The `configFiles` folder is essential for the correct functioning of the pipeline. It should contain three configuration files: `parameters.JSON`, `paths.JSON` and `steps.JSON`.

## Configuration Files

These three files contain all the information necessary to run the pipeline.

### Quality Control Pipeline Parameters

The `parameters.JSON` file contains values for `PLINK` commands that will be used in the pipeline as well as other parameters to tailor other steps. The parameters for the CLI (command line interface) must be provided in a `.JSON` file with the following structure:

```
{
    "sample_qc": {
        "rename_snp"   : true,
        "hh_to_missing": true,
        "use_kingship" : true,
        "ind_pair"     : [50, 5, 0.2],
        "mind"         : 0.2,
        "sex_check"    : [0.2, 0.8],
        "maf"          : 0.01,
        "het_deviation": 3,
        "kingship"     : 0.354,
        "ibd_threshold": 0.185
    },
    "ancestry_qc": {
        "ind_pair"     : [50, 5, 0.2],
        "pca"          : 10,
        "maf"          : 0.01,
        "ref_threshold": 4,
        "stu_threshold": 4,
        "reference_pop": "SAS",
        "num_pcs"      : 10,
    },
    "variant_qc": {
        "chr_y": 24,
        "miss_data_rate": 0.2,
        "diff_genotype_rate": 1e-5,
        "geno": 0.1,
        "maf": 5e-8,
        "hwe": 5e-8,
    },
    "umap_plot": {
        "umap_maf": 0.01,
        "umap_mind": 0.2,
        "umap_geno": 0.1,
        "umap_hwe": 5e-8,
        "umap_ind_pair": [50, 5, 0.2],
        "umap_pca": 10,
        "n_neighbors": [5, 10, 15],
        "metric": ["euclidean", "chebyshev"],
        "min_dist": [0.01, 0.1, 0.2],
        "random_state": 42,
        "case_control_marker": true,
        "color_hue_file": "path/to/color_hue_file.txt",
        "umap_kwargs": {}
    }
}
```

The values that come with each parameter are the default values used in our research group. If the user wishes to change at least one of them, please provide the full information in the configuration file.

### Paths to Project Folders

The `paths.JSON` file contains the addresses to the project folder as well as the prefix of the input and output data. The file must contain the following fields:

```
{
    "input_directory"      : "<path to folder with project input data>",
    "input_prefix"         : "<prefix of the input data>",
    "output_directory"     : "<path to folder where the output data will go>",
    "output_prefix"        : "<prefix for the output data>",
    "high_ld_file"         : "<path to file with high LD regions>"
}
```

If the CLI is run locally the user should provide the full path to file and directories. If no high LD file is provided or if the path is wrong, the library will use the one it has by default.

### Pipeline Steps

The `steps.JSON` file has the following structure:

```
{
    "ancestry": true,
    "sample"  : true,
    "variant" : true,
    "umap"    : true
}
```

With the above configuration, all three steps will run seamlessly, which is the recommended initial configuration. If you want to skip some steps, change the value to `false`. For example,

```
{
    "sample"   : false,
    "ancestry" : false,
    "variant"  : true,
    "umap"     : true
}
```

allows you to run only the variant QC and generate the UMAP plot(s). Note that an exception will be raised if the ancestry cehck step has not been run, as the necessary files for the variant step would not be available.

## Dependable Files

This folder should contain additional files to run the quality control pipeline. For example, the user might use this directory to store the high LD regions files in case it wants to use a different one from the library's default. Moreover, if the user wants to explore the population structure with respect to some category, the corresponding file should be located in this folder.

```
dependables
    |
    |---high-LD_regions.txt
    |
    |---population_categories.txt
```

Regarding the `population_structure.txt`, we expect a file with three colums, the first two are the ones for the `IID` and `FID` from **PLINK** `.fam` file, and the third one with the category that wants to be explored.

The other external files needed to perform the QC pipeline are the reference genome files. The library has the facility of fetch and process the reference genome automatically.

## Output Data

This folder has the following structure:
```
outputData
    |
    |---ancestry_results
    |
    |---umap_plots
    |
    |---sample_qc_results
    |
    |---variant_qc_results
```

### Results of ancestry outliers analysis

This folder contains the results from the ancestry analysis. Once the process is finished the folder will contain three folders and several files (we intend to reduce the files at a leter step). The three folders are
1. `fail_samples`: it contains a `.txt` file with the samples that failed the ancestry check; 
2. `clean_files`: it contains the cleaned files after the ancestry check in `PLINK` format;
3. `ancestryQC_plots`: it contains two plots showing the PCA decomposition of the study population against the reference panel.
The files are those resulting from the several steps of the ancestry check.

Recall that the cleaned binary files will feed the next steps.

### UMAP Plots

In this folder one can find the plot(s) generated after the UMAP dimensionality reduction, in order to explore the structure of the study population.

### Results of Sample Quality Control

This folder contains the results from the Sample Quality Control. Once the process is done the folder will contain three folders and multiple files. The three folders are
1. `fail_samples`: it contains `.txt` files with the samples that failed the different stages of the sample QC; 
2. `clean_files`: it contains the cleaned files after the sample quality control in `PLINK` format;
3. `sampleQC_plots`: it contains different plots that serve as a report of the different stages and might suggest a different selection of parameters.
The files are those resulting from the several steps of the sample QC.

Recall that the cleaned binary files will feed the next steps.

### Results of Variant Quality Control

This folder contains the results from the Variant Quality Control. Once the process is done the folder will contain three folders and several files. The three folders are
1. `fail_samples`: it contains `.txt` files with the samples that failed the different stages of the variant QC; 
2. `clean_files`: it contains the cleaned files after the variant quality control in `PLINK` format;
3. `variantQC_plots`: it contains different plots that serve as a report of the different stages and might suggest a different selection of parameters.
The files are those resulting from the steps of the varaint QC.

These cleaned binary files are ready for the next steps of the GWAS analysis.

## Installation and usage

The library can be installed by cloning the GitHub repository:

```
git clone https://github.com/cge-tubingens/IDEAL-GENOM-QC.git
```

or directly from PyPI:

```
pip install ideal_genom_qc
```

It is important to remark that the version in PyPI is the stable one, while the one on GitHub is on development.

### Setting up the environment

The virtual environment can be created using either `Poetry` or `pip`. Since this is a `Poetry`-based project, we recommend using `Poetry`. Once `Poetry` is installed on your system (refer to [`Poetry` documentation](https://python-poetry.org/docs/) for installation details), navigate to the cloned repository folder and run the following command:

```
poetry install
```

It is important to remark that currently the project has been updated to use `Poetry 2.0`.

### Pipeline usage options

#### 1. Inside a virtual environment

After running the `poetry install` activate the virtual environment with 

```
poetry shell
```

 Once the environment is active, you can execute the pipeline with the following command:

```
python3 ideal_genom_qc --path_params <path to parameters.JSON> 
                             --file_folders <path to paths.JSON> 
                             --steps <path to steps.JSON>
                             --recompute-merge true
                             --built 38
```

The first three parameters are the path to the three configuration files. The fourth is used to control the pipeline behavior.

#### 2. Using `Poetry` directly

One of the benefits of using `Poetry` s that it eliminates the need to activate a virtual environment. Run the pipeline directly with:

```
poetry run python3 ideal_genom_qc --path_params <path to parameters.JSON> 
                             --file_folders <path to paths.JSON> 
                             --steps <path to steps.JSON>
                             --recompute-merge true
                             --built 38
```
#### 3. Jupyter Notebooks

The package includes Jupyter notebooks located in the notebooks folder. Each notebook corresponds to a specific step of the pipeline. Simply provide the required parameters to execute the steps interactively.

Using the notebooks is a great way to gain a deeper understanding of how the pipeline operates.

#### 4. Docker Container

A `Dockerfile` is provided to build a container for the pipeline. Since the container interacts with physical files, it is recommended to use the following command:

```
docker run -v <path to project folder>:/data <docker_image_name>:<tag> --path_params <relative path to parameters.JSON> --file_folders <relative path to paths.JSON> --steps <relative path to steps.JSON> ---recompute-merge true --built 38
```

It is important to remark that the path to the files in `paths.JSON` must be relative to their location inside `data` folder in the Docker container.



