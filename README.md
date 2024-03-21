# Genotype Quality Control Pipeline

This is a Python package designed to perform a genotype quality control pipeline. It encapsulates several years of research at CGE Tübingen.

## Basic requirements

The quality control pipeline is build on `PLINK 1.9` as main tool. The `cge-comrare-pipeline` works as a wrapper for the different pipeline steps. Then, to run the pipeline, `PLINK 1.9` must be installed in the system. 

The pipeline is designed to run seamlessly with a "minimum" input and get cleaned binary files as result. In order to accomplish this, it is expected the following folder structure:

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
1. The folder `inputData` should contain the binary files with the genotype data to analyze in `PLINK` format (`.bed`, `.bim`, `.fam` files).

2. The folder `outputData` will contain the resultant files of the quality control pipeline. Bellow it will be treated in detail the pipeline output.

3. The folder `dependables` is designed to contain necessary files for the pipeline. It should contain the file `high_LD_regions.txt`. In future versions it might contain additional files.

4. The folder `configFiles` is essential for the pipeline correct functioning. It should contain two configuration files: `parameters.JSON` and `paths.JSON`.

## Configuration files

These two files contain all the information necessary to run the pipeline.

### Quality control pipeline parameters

The file `parameters.JSON` contains values for `PLINK` commands that will be used in the pipeline. If this file is not provided, the default values of the pipeline will be taken into account. These are

```
    "maf" : 0.05,
    "geno": 0.1,
    "mind": 0.1,
    "hwe" : 0.00000005,
    "sex_check": [0.2, 0.8],
    "indep-pairwise": [50, 5, 0.2],
    "chr": 24,
    "outlier_threshold": 6,
    "pca": 10.
```

If one wants to change at least one of the default values, please provide the full information in the configuration file. In repository can be found the `.JSON` file corresponding to the `cge-comrare-pipeline` default parameters.

### Paths to project folders

The file `paths.JSON` contain the addresses to the project folder as well as the prefix of the input and output data. The file must contain the following fields:

```
{
    "input_directory"      : "<path to folder with project input data>",
    "input_prefix"         : "<prefix of the input data>",
    "output_directory"     : "<path to folder where the output data will go>",
    "output_prefix"        : "<prefix for the output data>",
    "dependables_directory": "<path to folder with dependables files>"
}
```

## Output data

## Usage

The pipeline is easy to use. Once installed in the system or in a virtual enviroment one needs to run the following command:

```
python3 cge_comrare_pipeline --path_params <path to parameters.JSON> --file_folders <path to paths.JSON>
```

