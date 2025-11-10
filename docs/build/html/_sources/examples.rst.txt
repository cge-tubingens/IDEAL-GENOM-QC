Examples
========

This page provides practical examples of using IDEAL-GENOM-QC for different types of genomic studies. Each example includes complete configuration files and step-by-step instructions.

Example 1: Basic Case-Control Study
-----------------------------------

This example demonstrates a standard GWAS quality control pipeline for a case-control study.

**Study Setup:**
- 2,000 samples (1,000 cases, 1,000 controls)
- 500,000 SNPs
- European population
- Standard QC thresholds

**Configuration Files:**

*parameters.json:*

.. code-block:: json

    {
        "sample_qc": {
            "rename_snp": true,
            "hh_to_missing": true,
            "use_kinship": true,
            "ind_pair": [50, 5, 0.2],
            "mind": 0.1,
            "sex_check": [0.2, 0.8],
            "maf": 0.01,
            "het_deviation": 3,
            "kinship": 0.354,
            "ibd_threshold": 0.185
        },
        "ancestry_qc": {
            "ind_pair": [50, 5, 0.2],
            "pca": 10,
            "maf": 0.05,
            "ref_threshold": 3,
            "stu_threshold": 3,
            "reference_pop": "EUR",
            "num_pcs": 10
        },
        "variant_qc": {
            "chr_y": 24,
            "miss_data_rate": 0.1,
            "diff_genotype_rate": 1e-4,
            "geno": 0.05,
            "maf": 0.01,
            "hwe": 1e-6
        },
        "umap_plot": {
            "umap_maf": 0.05,
            "umap_mind": 0.1,
            "umap_geno": 0.05,
            "umap_hwe": 1e-6,
            "umap_ind_pair": [50, 5, 0.2],
            "umap_pca": 10,
            "n_neighbors": [10, 15, 20],
            "metric": ["euclidean"],
            "min_dist": [0.1, 0.3, 0.5],
            "random_state": 42,
            "case_control_marker": true
        }
    }

*paths.json:*

.. code-block:: json

    {
        "input_directory": "/data/gwas_study/inputData",
        "input_prefix": "gwas_data",
        "output_directory": "/data/gwas_study/outputData",
        "output_prefix": "gwas_clean",
        "high_ld_file": "/data/gwas_study/dependables/high-LD-regions_GRCH38.txt"
    }

*steps.json:*

.. code-block:: json

    {
        "ancestry": true,
        "sample": true,
        "variant": true,
        "umap": true,
        "fst": true
    }

**Execution:**

.. code-block:: bash

    python -m ideal_genom_qc \\
        --path_params /data/gwas_study/configFiles/parameters.json \\
        --file_folders /data/gwas_study/configFiles/paths.json \\
        --steps /data/gwas_study/configFiles/steps.json \\
        --recompute-merge true \\
        --built 38

**Expected Results:**
- ~1,800 samples passing all QC steps
- ~450,000 SNPs after variant QC
- Clear population clustering in UMAP plots
- Distinct case/control visualization

Example 2: Multi-Ethnic Population Study
-----------------------------------------

This example shows QC for a diverse population study with multiple ethnicities.

**Study Setup:**
- 5,000 samples from multiple populations
- Mixed ancestry backgrounds
- Population structure analysis focus

**Configuration Highlights:**

.. code-block:: json

    {
        "ancestry_qc": {
            "ref_threshold": 6,
            "stu_threshold": 6,
            "reference_pop": "ALL",
            "num_pcs": 20
        },
        "umap_plot": {
            "n_neighbors": [5, 10, 15, 30],
            "metric": ["euclidean", "manhattan", "cosine"],
            "min_dist": [0.01, 0.1, 0.3],
            "case_control_marker": false,
            "color_hue_file": "/data/population_categories.txt"
        }
    }

**Population Categories File:**

.. code-block:: text

    # population_categories.txt
    FAM001  SAMPLE001  African
    FAM002  SAMPLE002  European  
    FAM003  SAMPLE003  East_Asian
    FAM004  SAMPLE004  Hispanic
    FAM005  SAMPLE005  Mixed

**Python API Usage:**

.. code-block:: python

    from ideal_genom_qc import AncestryQC, UMAPplot
    from pathlib import Path
    
    # Multi-ethnic ancestry analysis
    ancestry_qc = AncestryQC(
        input_path=Path("inputData"),
        input_name="multi_ethnic_data",
        output_path=Path("outputData"),
        output_name="ancestry_clean",
        built="38"
    )
    
    # Use more lenient thresholds for diverse populations
    ancestry_qc.run_ancestry_qc(
        ref_threshold=6,
        stu_threshold=6,
        reference_pop="ALL",
        num_pcs=20
    )
    
    # Create detailed UMAP plots
    umap_plotter = UMAPplot(
        input_path=Path("outputData/ancestry_results/clean_files"),
        input_name="ancestry_clean",
        output_path=Path("outputData/umap_plots"),
        color_hue_file=Path("dependables/population_categories.txt")
    )

Jupyter Notebook Examples
--------------------------

The package includes interactive Jupyter notebooks in the ``notebooks/`` folder:

- ``sample_qc.ipynb``: Interactive sample QC with live plotting
- ``ancestry_qc.ipynb``: Population structure analysis
- ``variant_qc.ipynb``: Variant-level quality control
- ``umap_plot.ipynb``: Advanced visualization techniques

**Running Notebooks:**

.. code-block:: bash

    # Install Jupyter
    pip install jupyter
    
    # Launch notebook server
    cd notebooks/
    jupyter notebook
    
    # Open desired notebook and follow instructions

**Hardware Recommendations:**

- **CPU**: 8+ cores for large datasets
- **RAM**: 16GB minimum, 64GB for biobank-scale
- **Storage**: SSD preferred, 100GB+ free space
- **Network**: Fast connection for reference data download

Next Steps
----------

- Explore the :doc:`configuration` guide for parameter tuning
- Check the :doc:`troubleshooting` guide for common issues
- Review the API documentation for advanced usage
- Join our community discussions on GitHub