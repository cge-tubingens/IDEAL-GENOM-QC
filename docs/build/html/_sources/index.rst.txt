.. ideal-genome-qc documentation master file, created by
   sphinx-quickstart on Fri Apr 11 15:12:15 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ideal-genome-qc documentation
=============================

**ideal-genom-qc** is a Python library for performing automated and reproducible quality control (QC) on genomic data, with a focus on ancestry, variant filtering, and sample-level QC. It is designed to integrate with modern bioinformatics tools and supports visualizations and report generation.

Features
--------
- Sample, Ancestry and Variant QC pipelines
- Additional population structure analysis using PCA and UMAP
- Full automated pipeline for QC
- Many plots to report the results
- Clean API and extensible modules

Installation
------------
You can install the latest version from PyPI:

.. code-block:: bash

    pip install ideal-genom-qc

Or clone the repository and install it locally: (Disclaimer: this version could have some errors because contains the latest changes)

.. code-block:: bash

   git clone https://github.com/cge-tubingens/IDEAL-GENOM-QC.git


Documentation
-------------


.. toctree::
   :maxdepth: 2

   SampleQC
   AncestryQC
   UMAPplot
   VariantQC
   Helpers
   get_references

