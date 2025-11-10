Configuration Guide
===================

This guide provides detailed information about configuring IDEAL-GENOM-QC for your specific needs. Understanding these parameters will help you optimize the quality control process for your genomic data.

Configuration Files Overview
-----------------------------

IDEAL-GENOM-QC uses three JSON configuration files:

1. **parameters.json**: QC thresholds and algorithm parameters
2. **paths.json**: File and directory paths
3. **steps.json**: Pipeline steps to execute

Parameters Configuration
------------------------

The ``parameters.json`` file controls the behavior of each QC module. Each section corresponds to a specific QC step.

Sample QC Parameters
^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
        "sample_qc": {
            "rename_snp": true,
            "hh_to_missing": true,
            "use_kinship": true,
            "ind_pair": [50, 5, 0.2],
            "mind": 0.2,
            "sex_check": [0.2, 0.8],
            "maf": 0.01,
            "het_deviation": 3,
            "kinship": 0.354,
            "ibd_threshold": 0.185
        }
    }

**Parameter Descriptions:**

- ``rename_snp`` (bool): Whether to rename SNPs for consistency
- ``hh_to_missing`` (bool): Convert heterozygous haploid calls to missing
- ``use_kinship`` (bool): Use KING kinship estimation for relatedness
- ``ind_pair`` (list): LD pruning parameters [window_size, step_size, r2_threshold]
- ``mind`` (float): Maximum missing data rate per individual (0.0-1.0)
- ``sex_check`` (list): F coefficient thresholds [female_max, male_min]
- ``maf`` (float): Minor allele frequency threshold
- ``het_deviation`` (float): Standard deviations from mean heterozygosity
- ``kinship`` (float): Kinship coefficient threshold for relatedness
- ``ibd_threshold`` (float): IBD threshold for identifying duplicates

**Recommended Values:**

We have recommended the default values based on our experience with case-control studies. Adjust these parameters based on your study design and data quality. It is advisable to run the pipeline with different thresholds and visually inspect the results.


Ancestry QC Parameters
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
        "ancestry_qc": {
            "ind_pair": [50, 5, 0.2],
            "pca": 10,
            "maf": 0.01,
            "ref_threshold": 4,
            "stu_threshold": 4,
            "reference_pop": "SAS",
            "num_pcs": 10
        }
    }

**Parameter Descriptions:**

- ``ind_pair``: LD pruning parameters for PCA
- ``pca``: Number of principal components to compute
- ``maf``: MAF threshold for variants used in PCA
- ``ref_threshold``: SD threshold for reference population outliers
- ``stu_threshold``: SD threshold for study population outliers
- ``reference_pop``: Reference population code (EUR, AFR, AMR, EAS, SAS)
- ``num_pcs``: Number of PCs to use for ancestry assignment

**Reference Population Codes:**

- ``EUR``: European
- ``AFR``: African
- ``AMR``: Admixed American
- ``EAS``: East Asian
- ``SAS``: South Asian

Variant QC Parameters
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
        "variant_qc": {
            "chr_y": 24,
            "miss_data_rate": 0.2,
            "diff_genotype_rate": 1e-5,
            "geno": 0.1,
            "maf": 5e-8,
            "hwe": 5e-8
        }
    }

**Parameter Descriptions:**

- ``chr_y``: Chromosome Y identifier (23 or 24)
- ``miss_data_rate``: Maximum missing data rate per variant
- ``diff_genotype_rate``: Differential missing rate between cases/controls
- ``geno``: Maximum missing genotype rate
- ``maf``: Minor allele frequency threshold (very low for rare variants)
- ``hwe``: Hardy-Weinberg equilibrium p-value threshold


UMAP Plotting Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
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
            "color_hue_file": null,
            "umap_kwargs": {}
        }
    }

**Parameter Descriptions:**

- ``n_neighbors``: List of neighbor values for UMAP
- ``metric``: Distance metrics to try
- ``min_dist``: Minimum distances for UMAP embedding
- ``random_state``: Random seed for reproducibility
- ``case_control_marker``: Whether to color by case/control status
- ``color_hue_file``: Optional file for custom sample coloring
- ``umap_kwargs``: Additional UMAP parameters

Paths Configuration
-------------------

The ``paths.json`` file specifies all file and directory locations:

.. code-block:: json

    {
        "input_directory": "/full/path/to/inputData",
        "input_prefix": "mydata",
        "output_directory": "/full/path/to/outputData", 
        "output_prefix": "clean_data",
        "high_ld_file": "/full/path/to/high-LD-regions.txt"
    }

**Path Guidelines:**

- Use **absolute paths** for reliability
- Ensure directories exist before running
- Use consistent naming conventions
- Avoid spaces and special characters in paths

**Docker Paths:**

When using Docker, paths should be relative to the container's ``/data`` directory:

.. code-block:: json

    {
        "input_directory": "/data/inputData",
        "input_prefix": "mydata",
        "output_directory": "/data/outputData",
        "output_prefix": "clean_data",
        "high_ld_file": "/data/dependables/high-LD-regions.txt"
    }

Steps Configuration
-------------------

The ``steps.json`` file controls which pipeline steps to execute:

.. code-block:: json

    {
        "ancestry": true,
        "sample": true,
        "variant": true,
        "umap": true,
        "fst": true
    }

**Step Dependencies:**

- ``sample`` → ``ancestry`` → ``variant`` → ``umap`` → ``fst``
- You can skip steps, but maintain dependencies
- Results from previous steps are required for subsequent steps

**Common Configurations:**

*Full Pipeline:*

.. code-block:: json

    {
        "ancestry": true,
        "sample": true,
        "variant": true,
        "umap": true,
        "fst": true
    }

*Sample QC Only:*

.. code-block:: json

    {
        "ancestry": false,
        "sample": true,
        "variant": false,
        "umap": false,
        "fst": false
    }

Advanced Configuration
----------------------

Custom LD Regions
^^^^^^^^^^^^^^^^^^

Provide your own high-LD regions file:

.. code-block:: text

    # high-LD-regions.txt format
    1   48000000    52000000    # Chromosome, start, end
    2   85000000    100000000
    6   25000000    35000000

Population Categories
^^^^^^^^^^^^^^^^^^^^^

For custom population visualization, create a population file:

.. code-block:: text

    # population_categories.txt format
    FID1    IID1    EUR
    FID2    IID2    AFR
    FID3    IID3    SAS

Performance Tuning
-------------------

**Memory Optimization:**

- Increase ``ind_pair`` window size for large datasets
- Reduce ``pca`` components if memory is limited
- Process chromosomes separately for very large datasets

**Speed Optimization:**

- Use SSD storage for temporary files
- Increase available CPU cores
- Consider splitting large datasets

**Disk Space Management:**

- Monitor intermediate file sizes
- Clean up temporary files regularly
- Use compression for archival storage

Best Practices
--------------

1. **Version Control**: Keep configuration files under version control
2. **Documentation**: Document parameter choices and rationale
3. **Validation**: Always validate results visually
4. **Backup**: Keep copies of successful configurations
5. **Testing**: Test parameter changes on small datasets first

Troubleshooting
---------------

**Common Configuration Issues:**

- **Path not found**: Check absolute paths and permissions
- **Parameter out of range**: Verify threshold values are reasonable
- **JSON syntax errors**: Validate JSON format
- **Memory errors**: Reduce dataset size or adjust parameters

See the :doc:`troubleshooting` guide for more detailed solutions.