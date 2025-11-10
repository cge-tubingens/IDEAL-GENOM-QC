.. ideal-genome-qc documentation master file, created by
   sphinx-quickstart on Fri Apr 11 15:12:15 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

IDEAL-GENOM-QC Documentation
============================

.. image:: https://readthedocs.org/projects/ideal-genom-qc/badge/?version=latest
   :target: https://ideal-genom-qc.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/ideal-genom-qc.svg
   :target: https://pypi.org/project/ideal-genom-qc/
   :alt: PyPI version

**IDEAL-GENOM-QC** is a comprehensive Python library designed for performing automated and reproducible quality control (QC) on human genomic data, specifically tailored for case-control studies. This package serves as a sophisticated wrapper around PLINK tools, encapsulating years of research expertise from CGE Tübingen.

🎯 **Key Features**
------------------

✨ **Comprehensive QC Pipeline**
   - **Sample QC**: Automated sample quality assessment and filtering
   - **Ancestry QC**: Population structure analysis and outlier detection
   - **Variant QC**: SNP-level quality control and filtering
   - **Visualization**: Rich plotting capabilities for QC reporting

🔬 **Advanced Analytics**
   - **PCA Analysis**: Principal Component Analysis for population structure
   - **UMAP Plotting**: Modern dimensionality reduction for data exploration
   - **Kinship Analysis**: Relationship detection and handling
   - **Sex Check**: Automated sex verification from genetic data

🛠️ **User-Friendly Design**
   - **Minimal Configuration**: Ready-to-use with sensible defaults
   - **Flexible Parameters**: Easily customizable for different studies
   - **Multiple Interfaces**: Command-line, Python API, and Jupyter notebooks
   - **Docker Support**: Containerized deployment option

📊 **Rich Reporting**
   - **Automated Plots**: Publication-ready visualizations
   - **Quality Metrics**: Comprehensive QC statistics
   - **Clean Data Output**: PLINK-formatted results ready for GWAS

Quick Start
-----------

**Installation**

Install from PyPI (recommended for stable version):

.. code-block:: bash

    pip install ideal-genom-qc

Or install the development version:

.. code-block:: bash

    git clone https://github.com/cge-tubingens/IDEAL-GENOM-QC.git
    cd IDEAL-GENOM-QC
    poetry install

**Basic Usage**

.. code-block:: bash

    # Run the complete QC pipeline
    python -m ideal_genom_qc --path_params config/parameters.json \
                             --file_folders config/paths.json \
                             --steps config/steps.json \
                             --recompute-merge true \
                             --built 38

**Python API**

.. code-block:: python

    from ideal_genom_qc import SampleQC, AncestryQC, VariantQC
    
    # Initialize sample QC
    sample_qc = SampleQC(
        input_path="data/input",
        input_name="mydata",
        output_path="data/output",
        output_name="clean_data",
        high_ld_file="data/high_ld_regions.txt"
    )
    
    # Run QC pipeline
    sample_qc.run_sample_qc()

Documentation Contents
======================

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   getting_started
   configuration
   examples

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   SampleQC
   AncestryQC
   PopStructure
   VariantQC
   Helpers
   get_references

.. toctree::
   :maxdepth: 1
   :caption: Additional Resources

   faq
   troubleshooting
   contributing

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

