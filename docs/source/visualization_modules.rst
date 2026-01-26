Visualization Modules
=====================

The ``ideal_genom.visualizations`` package provides functions for creating publication-ready plots for GWAS and genomic analysis.

Module Overview
---------------

.. contents:: Modules
   :local:
   :depth: 1

manhattan_type
--------------

Generate Manhattan and Miami plots for genome-wide association studies (GWAS).

**Module:** ``ideal_genom.visualizations.manhattan_type``

Features:
^^^^^^^^^

- Data processing and visualization of GWAS summary statistics
- Annotation of SNPs with gene information from various sources
- Highlighting and labeling of specific SNPs of interest
- Support for both Manhattan (single study) and Miami (two studies) plots

Key Functions:
^^^^^^^^^^^^^^

.. py:function:: compute_relative_pos(data, chr_col='CHR', pos_col='POS', p_col='p')

   Compute the relative position of probes/SNPs across chromosomes and add a -log10(p-value) column.

   :param data: Input DataFrame containing genomic data
   :type data: pandas.DataFrame
   :param chr_col: Column name for chromosome identifiers
   :type chr_col: str
   :param pos_col: Column name for base pair positions
   :type pos_col: str
   :param p_col: Column name for p-values
   :type p_col: str
   :return: DataFrame with added columns for relative positions and -log10(p-values)
   :rtype: pandas.DataFrame

.. py:function:: manhattan(df_gwas, plots_dir, pval_col='P', chr_col='CHR', pos_col='POS', snp_col='SNP', p_threshold=5e-8, annotate=None, annotation_type='ensembl', genome_build='38', api_request=True, highlight_snps=None, alpha=0.7, save_name='manhattan.jpeg', colors=None, chr_text_shift=None, fig_size=(18, 6), dpi=500)

   Generate a Manhattan plot from GWAS summary statistics.

   :param df_gwas: DataFrame containing GWAS results
   :type df_gwas: pandas.DataFrame
   :param plots_dir: Directory path where the plot will be saved
   :type plots_dir: str
   :param pval_col: Column name for p-values
   :type pval_col: str
   :param chr_col: Column name for chromosome
   :type chr_col: str
   :param pos_col: Column name for base pair position
   :type pos_col: str
   :param snp_col: Column name for SNP identifiers
   :type snp_col: str
   :param p_threshold: Genome-wide significance threshold
   :type p_threshold: float
   :param annotate: List of SNP IDs to annotate with gene names
   :type annotate: Optional[list]
   :param annotation_type: Source for gene annotation ('ensembl', 'refseq', or 'both')
   :type annotation_type: str
   :param genome_build: Genome build version ('37' or '38')
   :type genome_build: str
   :param api_request: Whether to use API for annotation (if False, uses local GTF)
   :type api_request: bool
   :param highlight_snps: List of SNP IDs to highlight in different color
   :type highlight_snps: Optional[list]
   :param alpha: Transparency level for points
   :type alpha: float
   :param save_name: Filename for saving the plot
   :type save_name: str
   :param colors: Custom colors for alternating chromosomes
   :type colors: Optional[list]
   :param chr_text_shift: Shift amount for chromosome labels
   :type chr_text_shift: Optional[float]
   :param fig_size: Figure size (width, height) in inches
   :type fig_size: tuple
   :param dpi: Resolution for saved figure
   :type dpi: int
   :return: True if successful
   :rtype: bool

.. py:function:: miami(df_gwas1, df_gwas2, plots_dir, pval_col='P', chr_col='CHR', pos_col='POS', snp_col='SNP', p_threshold=5e-8, annotate1=None, annotate2=None, annotation_type='ensembl', genome_build='38', api_request=True, highlight_snps1=None, highlight_snps2=None, alpha=0.7, save_name='miami.jpeg', colors=None, chr_text_shift=None, fig_size=(18, 12), dpi=500, plot1_label='Study 1', plot2_label='Study 2')

   Generate a Miami plot (back-to-back Manhattan plots) comparing two GWAS studies.

   :param df_gwas1: DataFrame containing GWAS results for first study
   :type df_gwas1: pandas.DataFrame
   :param df_gwas2: DataFrame containing GWAS results for second study
   :type df_gwas2: pandas.DataFrame
   :param plots_dir: Directory path where the plot will be saved
   :type plots_dir: str
   :param plot1_label: Label for the top plot
   :type plot1_label: str
   :param plot2_label: Label for the bottom plot
   :type plot2_label: str
   :return: True if successful
   :rtype: bool

   *Other parameters are the same as manhattan() function*

Usage Example:
^^^^^^^^^^^^^^

.. code-block:: python

   import pandas as pd
   from ideal_genom.visualizations.manhattan_type import manhattan, miami
   
   # Load GWAS summary statistics
   gwas_df = pd.read_csv("gwas_results.txt", sep="\t")
   
   # Generate Manhattan plot
   manhattan(
       df_gwas=gwas_df,
       plots_dir="./plots",
       pval_col='P',
       chr_col='CHR',
       pos_col='BP',
       snp_col='SNP',
       p_threshold=5e-8,
       annotate=['rs12345', 'rs67890'],  # Annotate specific SNPs
       annotation_type='ensembl',
       genome_build='38',
       save_name='my_manhattan.jpeg'
   )
   
   # Generate Miami plot comparing two studies
   gwas_df2 = pd.read_csv("gwas_results2.txt", sep="\t")
   miami(
       df_gwas1=gwas_df,
       df_gwas2=gwas_df2,
       plots_dir="./plots",
       plot1_label='Discovery cohort',
       plot2_label='Replication cohort',
       save_name='my_miami.jpeg'
   )

plots
-----

Functions for generating various plots for GWAS data analysis.

**Module:** ``ideal_genom.visualizations.plots``

Features:
^^^^^^^^^

- QQ plots for visualizing the distribution of p-values
- Beta-beta scatter plots for comparing effect sizes between studies
- Trumpet plots for visualizing power and effect sizes
- Support for both binary and quantitative traits

Key Functions:
^^^^^^^^^^^^^^

.. py:function:: qqplot_draw(df_gwas, plots_dir, lambda_val=None, pval_col='P', conf_color='lightgray', save_name='qq_plot.jpeg', fig_size=(10, 10), dpi=500)

   Create a Q-Q (Quantile-Quantile) plot from GWAS results.

   This function generates a Q-Q plot comparing observed vs expected -log10(p-values)
   from GWAS results, including confidence intervals and genomic inflation factor (λ).

   :param df_gwas: DataFrame containing GWAS results with p-values
   :type df_gwas: pandas.DataFrame
   :param plots_dir: Directory path where the plot will be saved
   :type plots_dir: str
   :param lambda_val: Genomic inflation factor (calculated if None)
   :type lambda_val: Optional[float]
   :param pval_col: Column name for p-values
   :type pval_col: str
   :param conf_color: Color for confidence interval bands
   :type conf_color: str
   :param save_name: Filename for saving the plot
   :type save_name: str
   :param fig_size: Figure size (width, height) in inches
   :type fig_size: tuple
   :param dpi: Resolution for saved figure
   :type dpi: int
   :return: True if successful
   :rtype: bool

.. py:function:: beta_beta_plot(df1, df2, plots_dir, beta_col='BETA', se_col='SE', snp_col='SNP', pval_col='P', p_threshold=5e-8, annotate=None, annotation_type='ensembl', genome_build='38', api_request=True, save_name='beta_beta.jpeg', fig_size=(10, 10), dpi=500, x_label='Study 1', y_label='Study 2')

   Create a beta-beta scatter plot comparing effect sizes between two GWAS studies.

   :param df1: DataFrame containing GWAS results for first study
   :type df1: pandas.DataFrame
   :param df2: DataFrame containing GWAS results for second study
   :type df2: pandas.DataFrame
   :param plots_dir: Directory path where the plot will be saved
   :type plots_dir: str
   :param beta_col: Column name for effect sizes (beta)
   :type beta_col: str
   :param se_col: Column name for standard errors
   :type se_col: str
   :param snp_col: Column name for SNP identifiers
   :type snp_col: str
   :param pval_col: Column name for p-values
   :type pval_col: str
   :param p_threshold: Significance threshold for highlighting SNPs
   :type p_threshold: float
   :param annotate: List of SNP IDs to annotate
   :type annotate: Optional[list]
   :param annotation_type: Source for gene annotation
   :type annotation_type: str
   :param genome_build: Genome build version
   :type genome_build: str
   :param api_request: Whether to use API for annotation
   :type api_request: bool
   :param save_name: Filename for saving the plot
   :type save_name: str
   :param fig_size: Figure size (width, height) in inches
   :type fig_size: tuple
   :param dpi: Resolution for saved figure
   :type dpi: int
   :param x_label: Label for x-axis
   :type x_label: str
   :param y_label: Label for y-axis
   :type y_label: str
   :return: True if successful
   :rtype: bool

.. py:function:: trumpet_plot_binary(df_gwas, plots_dir, beta_col='BETA', se_col='SE', snp_col='SNP', pval_col='P', p_threshold=5e-8, annotate=None, annotation_type='ensembl', genome_build='38', api_request=True, maf_col='MAF', n_cases=None, n_controls=None, prevalence=0.5, alpha_val=0.05, save_name='trumpet_binary.jpeg', fig_size=(10, 10), dpi=500)

   Create a trumpet plot for binary traits, showing power curves and effect sizes.

   :param df_gwas: DataFrame containing GWAS results
   :type df_gwas: pandas.DataFrame
   :param plots_dir: Directory path where the plot will be saved
   :type plots_dir: str
   :param beta_col: Column name for effect sizes
   :type beta_col: str
   :param se_col: Column name for standard errors
   :type se_col: str
   :param snp_col: Column name for SNP identifiers
   :type snp_col: str
   :param pval_col: Column name for p-values
   :type pval_col: str
   :param p_threshold: Significance threshold
   :type p_threshold: float
   :param annotate: List of SNP IDs to annotate
   :type annotate: Optional[list]
   :param maf_col: Column name for minor allele frequency
   :type maf_col: str
   :param n_cases: Number of cases in the study
   :type n_cases: Optional[int]
   :param n_controls: Number of controls in the study
   :type n_controls: Optional[int]
   :param prevalence: Disease prevalence
   :type prevalence: float
   :param alpha_val: Significance level for power calculation
   :type alpha_val: float
   :param save_name: Filename for saving the plot
   :type save_name: str
   :return: True if successful
   :rtype: bool

.. py:function:: trumpet_plot_quantitative(df_gwas, plots_dir, beta_col='BETA', se_col='SE', snp_col='SNP', pval_col='P', p_threshold=5e-8, annotate=None, annotation_type='ensembl', genome_build='38', api_request=True, maf_col='MAF', n_samples=None, alpha_val=0.05, save_name='trumpet_quantitative.jpeg', fig_size=(10, 10), dpi=500)

   Create a trumpet plot for quantitative traits, showing power curves and effect sizes.

   :param df_gwas: DataFrame containing GWAS results
   :type df_gwas: pandas.DataFrame
   :param plots_dir: Directory path where the plot will be saved
   :type plots_dir: str
   :param n_samples: Total number of samples in the study
   :type n_samples: Optional[int]
   :return: True if successful
   :rtype: bool

   *Other parameters are the same as trumpet_plot_binary() function*

Usage Example:
^^^^^^^^^^^^^^

.. code-block:: python

   import pandas as pd
   from ideal_genom.visualizations.plots import (
       qqplot_draw, beta_beta_plot, 
       trumpet_plot_binary, trumpet_plot_quantitative
   )
   
   # Load GWAS results
   gwas_df = pd.read_csv("gwas_results.txt", sep="\t")
   
   # Generate QQ plot
   qqplot_draw(
       df_gwas=gwas_df,
       plots_dir="./plots",
       pval_col='P',
       save_name='my_qq_plot.jpeg'
   )
   
   # Beta-beta plot comparing two studies
   gwas_df2 = pd.read_csv("gwas_results2.txt", sep="\t")
   beta_beta_plot(
       df1=gwas_df,
       df2=gwas_df2,
       plots_dir="./plots",
       beta_col='BETA',
       se_col='SE',
       x_label='Discovery',
       y_label='Replication',
       save_name='my_beta_beta.jpeg'
   )
   
   # Trumpet plot for binary trait
   trumpet_plot_binary(
       df_gwas=gwas_df,
       plots_dir="./plots",
       n_cases=1000,
       n_controls=1000,
       prevalence=0.1,
       save_name='my_trumpet_binary.jpeg'
   )
   
   # Trumpet plot for quantitative trait
   trumpet_plot_quantitative(
       df_gwas=gwas_df,
       plots_dir="./plots",
       n_samples=2000,
       save_name='my_trumpet_quant.jpeg'
   )

zoom_heatmap
------------

Create zoomed heatmap visualizations of SNP associations, gene annotations, and linkage disequilibrium (LD) patterns.

**Module:** ``ideal_genom.visualizations.zoom_heatmap``

Features:
^^^^^^^^^

- Filter and annotate SNP data in a genomic region
- Calculate LD matrices using PLINK
- Generate three-panel plots with:
  
  1. Association plot with SNPs colored by functional consequences
  2. Gene track showing gene locations and orientations
  3. LD heatmap showing correlation patterns between SNPs

Key Functions:
^^^^^^^^^^^^^^

.. py:function:: filter_sumstats(data_df, lead_snp, snp_col, p_col, pos_col, chr_col, pval_threshold=5e-8, radius=10e6)

   Filter GWAS summary statistics based on a lead SNP, p-value threshold and genomic region.

   :param data_df: DataFrame containing GWAS summary statistics
   :type data_df: pandas.DataFrame
   :param lead_snp: Lead SNP identifier to center the region around
   :type lead_snp: str
   :param snp_col: Column name for SNP identifiers
   :type snp_col: str
   :param p_col: Column name for p-values
   :type p_col: str
   :param pos_col: Column name for base pair positions
   :type pos_col: str
   :param chr_col: Column name for chromosome
   :type chr_col: str
   :param pval_threshold: P-value threshold for filtering
   :type pval_threshold: float
   :param radius: Genomic radius around lead SNP (in base pairs)
   :type radius: Union[float, int]
   :return: Filtered DataFrame
   :rtype: pandas.DataFrame

.. py:function:: compute_ld(plink_file, snp_list, output_dir, lead_snp=None, ld_window_kb=10000, ld_window_snps=10000, threads=1)

   Compute linkage disequilibrium matrix for a list of SNPs using PLINK.

   :param plink_file: Path to PLINK binary file prefix (without .bed/.bim/.fam)
   :type plink_file: Union[str, Path]
   :param snp_list: List of SNP IDs for LD calculation
   :type snp_list: list
   :param output_dir: Directory to save output files
   :type output_dir: Union[str, Path]
   :param lead_snp: Lead SNP for coloring (optional)
   :type lead_snp: Optional[str]
   :param ld_window_kb: LD window size in kilobases
   :type ld_window_kb: int
   :param ld_window_snps: LD window size in number of SNPs
   :type ld_window_snps: int
   :param threads: Number of threads for PLINK
   :type threads: int
   :return: LD matrix as DataFrame
   :rtype: pandas.DataFrame

.. py:function:: create_zoom_heatmap(sumstats_df, plink_file, lead_snp, output_dir, snp_col='SNP', chr_col='CHR', pos_col='BP', p_col='P', beta_col='BETA', pval_threshold=5e-8, radius=500000, ld_window_kb=1000, genome_build='38', annotation_type='ensembl', api_request=True, fig_size=(14, 12), dpi=300, threads=1, save_name='zoom_heatmap.png')

   Create a comprehensive zoom heatmap plot with association, gene track, and LD panels.

   :param sumstats_df: DataFrame containing GWAS summary statistics
   :type sumstats_df: pandas.DataFrame
   :param plink_file: Path to PLINK binary file prefix
   :type plink_file: Union[str, Path]
   :param lead_snp: Lead SNP identifier to center the plot
   :type lead_snp: str
   :param output_dir: Directory to save output files
   :type output_dir: Union[str, Path]
   :param snp_col: Column name for SNP identifiers
   :type snp_col: str
   :param chr_col: Column name for chromosome
   :type chr_col: str
   :param pos_col: Column name for base pair position
   :type pos_col: str
   :param p_col: Column name for p-values
   :type p_col: str
   :param beta_col: Column name for effect sizes
   :type beta_col: str
   :param pval_threshold: P-value threshold for filtering
   :type pval_threshold: float
   :param radius: Genomic radius around lead SNP (in base pairs)
   :type radius: Union[float, int]
   :param ld_window_kb: LD window size in kilobases
   :type ld_window_kb: int
   :param genome_build: Genome build version ('37' or '38')
   :type genome_build: str
   :param annotation_type: Source for gene annotation ('ensembl', 'refseq', or 'both')
   :type annotation_type: str
   :param api_request: Whether to use API for functional annotation
   :type api_request: bool
   :param fig_size: Figure size (width, height) in inches
   :type fig_size: tuple
   :param dpi: Resolution for saved figure
   :type dpi: int
   :param threads: Number of threads for PLINK
   :type threads: int
   :param save_name: Filename for saving the plot
   :type save_name: str
   :return: Path to saved figure
   :rtype: Path

Usage Example:
^^^^^^^^^^^^^^

.. code-block:: python

   import pandas as pd
   from pathlib import Path
   from ideal_genom.visualizations.zoom_heatmap import create_zoom_heatmap
   
   # Load GWAS summary statistics
   sumstats = pd.read_csv("gwas_results.txt", sep="\t")
   
   # Create zoom heatmap around a lead SNP
   create_zoom_heatmap(
       sumstats_df=sumstats,
       plink_file=Path("data/genotypes"),  # Without .bed/.bim/.fam extension
       lead_snp='rs12345',
       output_dir=Path("./plots"),
       snp_col='SNP',
       chr_col='CHR',
       pos_col='BP',
       p_col='P',
       beta_col='BETA',
       pval_threshold=5e-8,
       radius=500000,  # 500kb window
       genome_build='38',
       annotation_type='ensembl',
       api_request=True,
       save_name='rs12345_zoom.png'
   )

Notes
-----

**Dependencies:**
   - matplotlib
   - seaborn
   - pandas
   - numpy
   - textalloc (for label positioning)
   - pyensembl (for gene annotations)
   - PLINK 1.9 or 2.0 (for LD calculations)

**Annotation Sources:**
   All plotting functions support gene annotation from:
   
   - **Ensembl**: Via REST API or local GTF files
   - **RefSeq**: Via local GTF files
   - **Both**: Combined annotations from both sources

**Genome Builds:**
   Supported genome builds are GRCh37/hg19 ('37') and GRCh38/hg38 ('38')

**Output Formats:**
   - JPEG format for Manhattan, Miami, QQ, beta-beta, and trumpet plots
   - PNG format for zoom heatmaps (recommended for better quality with complex graphics)
   - All plots are publication-ready with customizable DPI

See Also
--------

- :doc:`gwas_modules` - GWAS analysis modules that generate data for visualization
- :doc:`Helpers` - Annotation utilities used by visualization functions
- :doc:`api_overview` - Complete API reference
