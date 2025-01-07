# Genotype Quality Control Pipeline

This Python package is designed to execute a genotype quality control pipeline, encapsulating several years of research at CGE Tübingen.

## Basic requirements

The quality control pipeline is built on `PLINK 1.9` as the main tool. The `ideal_genom_qc` serves as a wrapper for the various QC pipeline steps. To run the pipeline, `PLINK 1.9` must be installed on the system.

The pipeline is designed to seamlessly run with minimal input and produce cleaned binary files as a result as well as several plots along the way. To accomplish this, the following folder structure is expected:

```
projectFolder
    |
    |---inputData
    |
    |---outputData
    |
    |---dependables
    |
    |---configFiles
```
1. The `inputData` folder should contain the binary files with the genotype data to be analyzed in `PLINK` format (`.bed`, `.bim`, `.fam` files).

2. The `outputData` folder will contain the resultant files of the quality control pipeline. Below, the pipeline output will be detailed.

3. The `dependables` folder is designed to contain necessary files for the pipeline.

4. The `configFiles` folder is essential for the correct functioning of the pipeline. It should contain three configuration files: `parameters.JSON`, `paths.JSON` and `steps.JSON`.

## Configuration Files

These three files contain all the information necessary to run the pipeline.

### Quality Control Pipeline Parameters

The `parameters.JSON` file contains values for `PLINK` commands that will be used in the pipeline as well as other parameters to tailor other steps. The parameters for the CLI (command line interface) must be provided in a `.JSON` file with the following structure:

```
{
    "sample_qc": {
        "ind_pair": [
            50,
            5,
            0.2
        ],
        "mind": 0.2,
        "sex_check": [
            0.2,
            0.8
        ],
        "het_deviation": 3,
        "maf": 0.01,
        "kingship": 0.354,
        "ibd_threshold": 0.185
    },
    "ancestry_qc": {
        "maf": 0.01,
        "ind_pair": [
            50,
            5,
            0.2
        ],
        "pca": 10,
        "ref_threshold": 4,
        "stu_threshold": 4,
        "num_pca": 10
    },
    "variant_qc": {
        "chr-y": 24,
        "miss_data_rate": 0.2,
        "diff_genotype_rate": 1e-5,
        "geno": 0.1,
        "hwe": 5e-8,
        "maf": 5e-8
    },
    "umap_plot": {
    	"maf": 0.01,
    	"geno": 0.1,
    	"mind": 0.2,
    	"hwe": 5e-8,
    	"ind_pair": [
            50,
            5,
            0.2
        ],
        "pca": 10,
        "n_neighbors": [5],
        "min_dist": [0.5],
        "metric": ["euclidean"]
    }
}
```

If you wish to change at least one of the default values, please provide the full information in the configuration file.

### Paths to Project Folders

The `paths.JSON` file contains the addresses to the project folder as well as the prefix of the input and output data. The file must contain the following fields:

```
{
    "input_directory"      : "<path to folder with project input data>",
    "input_prefix"         : "<prefix of the input data>",
    "output_directory"     : "<path to folder where the output data will go>",
    "output_prefix"        : "<prefix for the output data>",
    "dependables_directory": "<path to folder with dependables files>"
}
```

If the CLI is run locally you should provide the full path to the directories.

### Pipeline Steps

The `steps.JSON` file has the following structure:

```
{
    "pca"    : true,
    "sample" : true,
    "variant": true,
    "umap"   : true
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

This folder should contain additional files to run the quality control pipeline. The structure inside the directory should be as follows:

```
dependables
    |
    |---all_phase3.bed
    |
    |---all_phase3.bim
    |
    |---all_phase3.fam
    |
    |---all_phase3.psam
    |
    |---high-LD-regions.txt
    |
    |---geographic_info.txt
```

Notice that the files `all_phase3.bed`, `all_phase3.bim`, `all_phase3.fam` and `all_phase3.psam` correspond to the 1000 Genomes phase 3. In addition, the file `high-LD-regions.txt` corresponds to the build 38, in order to be consistent with 1000 Genomes phase 3 build.

The file `geographic_info.txt` is optional. It should contain two columns: the first one for individual ids, and the second one with geographical (or another categorical information that allows to split the study population) information. With this file it can be built a plot to see the genes distribution by categories. If the file is not present the plot will not be generated. The separator between the two columns must be a blank space or tab.

In a future realease is intended to add this step in the pipeline.

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

### Setting up the environment

The virtual environment can be created using either `Poetry` or `pip`. Since this is a `Poetry`-based project, we recommend using `Poetry`. Once `Poetry` is installed on your system (refer to [`Poetry` documentation](https://python-poetry.org/docs/) for installation details), navigate to the cloned repository folder and run the following command:

```
poetry install
```

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
                             --use-kingship <true or false>
```

The first three parameters are the path to the three configuration files. The fourth is used to control the pipeline behavior.

#### 2. Using `Poetry` directly

One of the benefits of using `Poetry` s that it eliminates the need to activate a virtual environment. Run the pipeline directly with:

```
poetry run python3 ideal_genom_qc --path_params <path to parameters.JSON> 
                             --file_folders <path to paths.JSON> 
                             --steps <path to steps.JSON>
                             --use-kingship <true or false>
```
#### 3. Jupyter Notebooks

The package includes Jupyter notebooks located in the notebooks folder. Each notebook corresponds to a specific step of the pipeline. Simply provide the required parameters to execute the steps interactively.

Using the notebooks is a great way to gain a deeper understanding of how the pipeline operates.

#### 4. Docker Container

A `Dockerfile` is provided to build a container for the pipeline. Since the container interacts with physical files, it is recommended to use the following command:

```
docker run -v <path to project folder>:/data <docker_image_name>:<tag> --path_params <relative path to parameters.JSON> --file_folders <relative path to paths.JSON> --steps <relative path to steps.JSON> --use-kingship false
```

It is important to remark that the path to the files in `paths.JSON` must be relative to their location inside `data` folder in the Docker container.

