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
