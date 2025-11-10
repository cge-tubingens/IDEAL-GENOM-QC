Getting Started
===============

This guide will help you get up and running with IDEAL-GENOM-QC quickly. We'll walk through setting up your first quality control pipeline step by step.

Overview
--------

IDEAL-GENOM-QC follows a structured workflow:

1. **Project Setup**: Organize your data and configuration files
2. **Configuration**: Set parameters, paths, and pipeline steps
3. **Execution**: Run the QC pipeline
4. **Results**: Review outputs and visualizations

Project Structure
-----------------

IDEAL-GENOM-QC expects a specific folder structure:

.. code-block:: text

    my_project/
    ├── inputData/
    │   ├── mydata.bed
    │   ├── mydata.bim
    │   └── mydata.fam
    ├── outputData/
    ├── configFiles/
    │   ├── parameters.json
    │   ├── paths.json
    │   └── steps.json
    └── dependables/
        └── high-LD-regions.txt

**Folder Descriptions:**

- ``inputData/``: Contains your PLINK binary files (.bed, .bim, .fam)
- ``outputData/``: Will contain all pipeline results
- ``configFiles/``: Contains configuration files (required)
- ``dependables/``: Optional files like custom LD regions

Step 1: Prepare Your Data
-------------------------

Ensure your genetic data is in PLINK binary format:

- ``.bed``: Binary genotype file
- ``.bim``: Variant information file  
- ``.fam``: Sample information file

If your data is in a different format, convert it first:

.. code-block:: bash

    # Convert VCF to PLINK
    plink --vcf mydata.vcf --make-bed --out inputData/mydata

Step 2: Configuration Files
---------------------------

Create three JSON configuration files in the ``configFiles/`` directory:

**parameters.json** - QC Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
        "sample_qc": {
            "rename_snp": true,
            "hh_to_missing": true,
            "use_kingship": true,
            "ind_pair": [50, 5, 0.2],
            "mind": 0.2,
            "sex_check": [0.2, 0.8],
            "maf": 0.01,
            "het_deviation": 3,
            "kingship": 0.354,
            "ibd_threshold": 0.185
        },
        "ancestry_qc": {
            "ind_pair": [50, 5, 0.2],
            "pca": 10,
            "maf": 0.01,
            "ref_threshold": 4,
            "stu_threshold": 4,
            "reference_pop": "SAS",
            "num_pcs": 10
        },
        "variant_qc": {
            "chr_y": 24,
            "miss_data_rate": 0.2,
            "diff_genotype_rate": 1e-5,
            "geno": 0.1,
            "maf": 5e-8,
            "hwe": 5e-8
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
            "case_control_marker": true
        }
    }

**paths.json** - File Paths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
        "input_directory": "/full/path/to/my_project/inputData",
        "input_prefix": "mydata",
        "output_directory": "/full/path/to/my_project/outputData",
        "output_prefix": "clean_data",
        "high_ld_file": "/full/path/to/my_project/dependables/high-LD-regions.txt"
    }

The paths should be absolute paths to ensure proper file access. In addition, the high-LD-regions.txt is optional; if not provided (or if the path is not correct), the default file included with IDEAL-GENOM-QC will be used.

**steps.json** - Pipeline Steps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
        "ancestry": true,
        "sample": true,
        "variant": true,
        "umap": true,
        "fst": true
    }

In the current version, all steps must be included in the JSON file, even if set to false. It is important to remark that the steps will always be executed in the following order: Sample QC, Ancestry QC, Variant QC, UMAP Plotting, and FST Calculation. Moreover, each step depends on the previous one, meaning that if a step is set to true, all preceding steps must also be set to true or, at least, the resulting files from the previous step must be available in the output folder.

Step 3: Running the Pipeline
-----------------------------

**Command Line Interface:**

.. code-block:: bash

    python -m ideal_genom_qc \\
        --path_params /path/to/configFiles/parameters.json \\
        --file_folders /path/to/configFiles/paths.json \\
        --steps /path/to/configFiles/steps.json \\
        --recompute-merge true \\
        --built 38

**Using Poetry:**

.. code-block:: bash

    poetry run python -m ideal_genom_qc \\
        --path_params configFiles/parameters.json \\
        --file_folders configFiles/paths.json \\
        --steps configFiles/steps.json \\
        --recompute-merge true \\
        --built 38

**Python API:**

A shorter example using the Python API to run Sample QC:

.. code-block:: python

    from ideal_genom_qc import SampleQC, AncestryQC, VariantQC
    from pathlib import Path

    # Initialize components
    sample_qc = SampleQC(
        input_path=Path("inputData"),
        input_name="mydata",
        output_path=Path("outputData"),
        output_name="clean_data",
        high_ld_file=Path("dependables/high-LD-regions.txt")
    )

    # Run sample QC
    sample_qc.execute_sample_qc_pipeline(
        mind=0.2,
        sex_check=[0.2, 0.8],
        maf=0.01,
        het_deviation=3
    )

Step 4: Understanding Results
-----------------------------

After running the pipeline, your ``outputData/`` folder will contain:

.. code-block:: text

    outputData/
    ├── ancestry_results/
    │   ├── clean_files/
    │   ├── fail_samples/
    │   └── ancestryQC_plots/
    ├── sample_qc_results/
    │   ├── clean_files/
    │   ├── fail_samples/
    │   └── sampleQC_plots/
    ├── variant_qc_results/
    │   ├── clean_files/
    │   ├── fail_samples/
    │   └── variantQC_plots/
    └── umap_plots/

**Key Output Files:**

- ``clean_files/``: Final cleaned PLINK files ready for GWAS
- ``fail_samples/``: Lists of samples/variants that failed QC
- ``*QC_plots/``: Visualization reports for each QC step

**Quality Control Reports:**

Each QC step generates informative plots:

- **Sample QC**: Heterozygosity, missing data, kinship plots
- **Ancestry QC**: PCA plots showing population structure
- **Variant QC**: MAF, HWE, missing data distributions
- **UMAP Plots**: Population structure visualization

Example Workflow
----------------

Here's a complete example using the Python API:

.. code-block:: python

    #!/usr/bin/env python3
    from ideal_genom_qc import SampleQC, AncestryQC, VariantQC
    from pathlib import Path
    import json

    # Load configuration
    with open('configFiles/parameters.json') as f:
        params = json.load(f)

    # Initialize QC components
    base_path = Path(".")
    
    sample_qc = SampleQC(
        input_path=base_path / "inputData",
        input_name="mydata",
        output_path=base_path / "outputData",
        output_name="clean_data",
        high_ld_file=base_path / "dependables" / "high-LD-regions.txt"
    )

    # Run complete QC pipeline
    print("Running Sample QC...")
    sample_qc.run_sample_qc(**params['sample_qc'])

    print("Running Ancestry QC...")
    ancestry_qc = AncestryQC(
        input_path=base_path / "outputData" / "sample_qc_results" / "clean_files",
        input_name="clean_data_after_sample_qc",
        output_path=base_path / "outputData",
        output_name="clean_data",
        built="38"
    )
    ancestry_qc.run_ancestry_qc(**params['ancestry_qc'])

    print("QC Pipeline completed!")

Next Steps
----------

- Explore the :doc:`examples` for more advanced usage
- Check the :doc:`configuration` guide for parameter tuning
- Review the API documentation for each module
- See the :doc:`troubleshooting` guide if you encounter issues

**Need Help?**

- Check our `FAQ <faq.html>`_
- Browse the examples in the ``notebooks/`` folder
- Report issues on `GitHub <https://github.com/cge-tubingens/IDEAL-GENOM-QC/issues>`_