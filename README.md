# Genotype Quality Control Pipeline

This Python package is designed to execute a genotype quality control pipeline, encapsulating several years of research at CGE Tübingen.

## Basic requirements

The quality control pipeline is built on `PLINK 1.9` as the main tool. The `cge-comrare-pipeline` serves as a wrapper for the various pipeline steps. To run the pipeline, `PLINK 1.9` must be installed on the system.

The pipeline is designed to seamlessly run with minimal input and produce cleaned binary files as a result. To accomplish this, the following folder structure is expected:

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

4. The `configFiles` older is essential for the correct functioning of the pipeline. It should contain three configuration files: `parameters.JSON`, `paths.JSON` and `steps.JSON`.

## Configuration Files

These two files contain all the information necessary to run the pipeline.

### Quality Control Pipeline Parameters

The `parameters.JSON` file contains values for `PLINK` commands that will be used in the pipeline. If this file is not provided, the default values of the pipeline will be taken into account. These are

```
{
    "maf" : 0.05,
    "geno": 0.1,
    "mind": 0.1,
    "hwe" : 0.00000005,
    "sex_check": [0.2, 0.8],
    "indep-pairwise": [50, 5, 0.2],
    "chr": 24,
    "outlier_threshold": 6,
    "pca": 10
}
```

If you wish to change at least one of the default values, please provide the full information in the configuration file. The corresponding `.JSON` file for the `cge-comrare-pipeline` default parameters can be found in the repository.

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

### Pipeline Steps

The `steps.JSON` file has the following structure:

```
{
    "pca"    : true,
    "sample" : true,
    "variant": true
}
```

With the above configuration, all three steps will run seamlessly, which is the recommended initial configuration. If you want to skip some steps, change the value to `false`. For example,

```
{
    "pca"    : false,
    "sample" : true,
    "variant": true
}
```

allows you to run only the sample and variant quality control. Note that an exception will be raised if the PCA step has not been run, as the necessary files for the sample steps would not be available.

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
```

Notice that the files `all_phase3.bed`, `all_phase3.bim`, `all_phase3.fam` and `all_phase3.psam` correspond to the 1000 Genomes phase 3. In addition, the file `high-LD-regions.txt` corresponds to the build 38, in order to be consistent with 1000 Genomes phase 3 build.

## Output Data

## Usage

The pipeline is easy to use. Once installed in the system or in a virtual enviroment one needs to run the following command:

```
python3 cge_comrare_pipeline --path_params <path to parameters.JSON> 
                             --file_folders <path to paths.JSON> 
                             --steps <path to steps.JSON>
```

## Docker Container

The Docker container for the pipeline is under construction.
