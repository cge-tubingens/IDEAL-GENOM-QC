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


### Paths to project folders

## Outpu data

## Usage
