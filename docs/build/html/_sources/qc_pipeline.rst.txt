Quality Control Pipeline
=========================

The Quality Control (QC) pipeline in IDEAL-GENOM provides comprehensive tools for ensuring the quality and integrity of genotype data in case-control studies. The pipeline systematically addresses individual-level quality issues, population stratification, and variant-level quality metrics to produce analysis-ready data.

Overview
--------

The QC pipeline consists of three sequential components:

1. **Sample QC**: Individual-level quality control and relatedness analysis
2. **Ancestry QC**: Population structure detection and ancestry outlier removal
3. **Variant QC**: Marker-level quality filtering and Hardy-Weinberg equilibrium testing

**Key Features:**

- Automated multi-stage quality control workflow
- Sex verification and heterozygosity analysis
- Kinship and IBD-based relatedness detection
- Ancestry projection with 1000 Genomes reference
- Hardy-Weinberg equilibrium testing by case-control status
- Differential missingness and batch effect detection
- Resource-aware parallel processing
- Publication-ready QC plots and comprehensive reports

Use Cases
---------

This pipeline is ideal for:

- **GWAS Preparation**: Ensure data quality before association testing
- **Case-Control Studies**: Remove biased samples and variants

Prerequisites
-------------

**Required Software:**

- PLINK 1.9 (for functionalities not yet in PLINK 2.0)
- PLINK 2.0 (for most QC operations)
- Python 3.11+ with required packages (matplotlib, seaborn, pandas, numpy, scipy, umap-learn)

**Required Data:**

- PLINK binary files (.bed, .bim, .fam)
- Phenotype data in .fam file (case=2, control=1, missing=0/-9)
- 1000 Genomes reference data (auto-downloaded if not provided)

**Input Data Requirements:**

- Genome-wide genotype data (not pre-filtered)
- At least 50,000 SNPs for reliable QC metrics
- Both cases and controls with phenotype labels
- Autosomal chromosomes (1-22)

Quick Start
-----------

**1. Get the QC Configuration Template**

.. code-block:: bash

    # Copy the QC pipeline template
    cp yaml_configs/qc_pipeline_config_template.yaml my_qc_pipeline.yaml

**2. Edit Configuration**

Update paths to your data:

.. code-block:: yaml

    pipeline:
      name: "my_qc_pipeline"
      base_output_dir: "/data/qc_output"
      steps:
        - name: "sample_qc"
          enabled: true
          init_params:
            input_path: "/data/raw_data"
            input_name: "mydata"
            output_path: "${base_output_dir}/sample_qc"
            build: "38"

**3. Run the Pipeline**

.. code-block:: bash

    # Validate configuration
    ideal-genom validate --config my_qc_pipeline.yaml
    
    # Run the full QC pipeline
    ideal-genom run --config my_qc_pipeline.yaml

Pipeline Components
-------------------

1. Sample QC
~~~~~~~~~~~~

The first stage performs individual-level quality control to identify and remove problematic samples.

**Quality Metrics:**

- **Genotyping Call Rate**: Remove samples with high missingness (default: >2%)
- **Sex Verification**: Check X chromosome homozygosity against reported sex (F<0.2=female, F>0.8=male)
- **Heterozygosity Rate**: Identify samples with unusual heterozygosity (±3 SD from mean)
- **Relatedness**: Detect duplicates, parent-offspring, siblings using kinship coefficients

**Outputs:**

- List of samples passing QC filters
- Quality control plots (missingness, heterozygosity, sex check, kinship heatmap)
- Cleaned PLINK files ready for ancestry analysis
- Detailed QC metrics log

**Configuration Example:**

.. code-block:: yaml

    - name: "sample_qc"
      enabled: true
      init_params:
        input_path: "/data/raw_data"
        input_name: "mydata"
        output_path: "/data/qc_output/sample_qc"
        build: "38"
      run_params:
        geno_mind: 0.02          # Max individual missingness: 2%
        geno_snp: 0.05           # Max SNP missingness: 5% (for QC only)
        het_threshold: 3         # Heterozygosity SD threshold
        pihat_threshold: 0.1875  # Relatedness threshold (2nd degree)
        maf: 0.01                # MAF threshold for LD pruning

**Typical Sample Exclusions:**

- Genotyping call rate <98% (~0.5-2% of samples)
- Sex discordance (XX as male or XY as female, ~0.1-0.5%)
- Extreme heterozygosity (possible contamination, ~0.1-0.5%)
- Cryptic relatedness (duplicates/relatives, ~1-5%)

2. Ancestry QC
~~~~~~~~~~~~~~

The second stage projects samples onto reference populations to identify ancestry outliers and ensure population homogeneity.

**Methodology:**

1. **Reference Integration**: Merge study data with 1000 Genomes Phase 3 reference
2. **LD Pruning**: Prune to independent SNPs for ancestry analysis (r²<0.2)
3. **PCA Projection**: Compute principal components with reference populations
4. **Outlier Detection**: Identify samples >6 SD from population centroid

**Reference Populations:**

- **EUR**: European (British, Finnish, Spanish, Italian, Iberian)
- **AFR**: African (Yoruba, Luhya, Gambian, Mende, Esan)
- **EAS**: East Asian (Han Chinese, Japanese, Southern Han Chinese, Dai, Kinh)
- **SAS**: South Asian (Gujarati, Punjabi, Bengali, Tamil, Telugu)
- **AMR**: Admixed American (Mexican, Puerto Rican, Colombian, Peruvian)

**Outputs:**

- PCA plots showing study samples vs. reference populations
- List of ancestry outliers
- Cleaned PLINK files ready for variant QC
- Population assignment for each sample

**Configuration Example:**

.. code-block:: yaml

    - name: "ancestry_qc"
      enabled: true
      init_params:
        input_path: "${base_output_dir}/sample_qc"
        input_name: "mydata_SampleQCed"
        output_path: "${base_output_dir}/ancestry_qc"
        build: "38"
        reference_1kG: "/data/references/1000genomes_build_38"
        reference_1kG_name: "1kG_phase3_GRCh38"
      run_params:
        expected_population: "EUR"
        ancestry_filter: true       # Remove outliers
        sd_threshold: 6             # SD for outlier detection

**Typical Ancestry Exclusions:**

- Population outliers 3 SD from centroid (~0.5-2% for homogeneous cohorts) using the Chebyshev metric as default.
- Visual inspection of PCA plots recommended

3. Variant QC
~~~~~~~~~~~~~

The final stage performs marker-level quality control to ensure variant reliability.

**Quality Filters:**

- **SNP Call Rate**: Remove variants with high missingness (default: >2%)
- **Minor Allele Frequency (MAF)**: Filter rare variants (default: <1%)
- **Hardy-Weinberg Equilibrium**: Separate HWE testing in cases (p<1e-10) and controls (p<1e-6)
- **Differential Missingness**: Remove variants with case-control missingness differences (p<1e-5)

**Outputs:**

- Analysis-ready PLINK binary files
- Variant-level QC statistics
- Manhattan plots showing removed variants
- Summary report with variant counts

**Configuration Example:**

.. code-block:: yaml

    - name: "variant_qc"
      enabled: true
      init_params:
        input_path: "${base_output_dir}/ancestry_qc"
        input_name: "mydata_AncestryQCed"
        output_path: "${base_output_dir}/variant_qc"
        build: "38"
      run_params:
        geno: 0.02               # Max SNP missingness: 2%
        maf: 0.01                # Min MAF: 1%
        hwe_control: 1e-6        # HWE p-value in controls
        hwe_case: 1e-10          # HWE p-value in cases
        diff_miss: 1e-5          # Differential missingness p-value

**Typical Variant Exclusions:**

- Genotyping call rate <98% (~1-3% of variants)
- MAF <1% (~20-40% of variants, depending on array)
- HWE violations (~0.1-0.5% in controls, ~0.01-0.1% in cases)
- Differential missingness (~0.1-0.5%)
- Final retention: typically 60-80% of initial markers

Advanced Configuration
----------------------

Data Flow
~~~~~~~~~

The pipeline handles data flow automatically using variable substitution:

.. code-block:: yaml

    pipeline:
      name: "comprehensive_qc"
      base_output_dir: "/data/qc_results"
      steps:
        - name: "sample_qc"
          init_params:
            output_path: "${base_output_dir}/step1_sample"
        
        - name: "ancestry_qc"
          init_params:
            # Automatically uses output from sample_qc
            input_path: "${base_output_dir}/step1_sample"
            output_path: "${base_output_dir}/step2_ancestry"
        
        - name: "variant_qc"
          init_params:
            # Automatically uses output from ancestry_qc
            input_path: "${base_output_dir}/step2_ancestry"
            output_path: "${base_output_dir}/step3_variant"

Resource Management
~~~~~~~~~~~~~~~~~~~

Control computational resources for each step:

.. code-block:: yaml

    - name: "sample_qc"
      init_params:
        threads: 8                  # CPU cores
        memory_gb: 16               # RAM allocation
      run_params:
        run_het: true               # Enable heterozygosity check
        run_ibd: true               # Enable relatedness analysis

Population-Specific Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adjust parameters for different populations:

.. code-block:: yaml

    # European cohort (default settings work well)
    - name: "ancestry_qc"
      run_params:
        expected_population: "EUR"
        sd_threshold: 6

    # Admixed cohort (more lenient thresholds)
    - name: "ancestry_qc"
      run_params:
        expected_population: "AMR"
        sd_threshold: 8           # More permissive
        ancestry_filter: false    # Keep all samples

    # African cohort
    - name: "ancestry_qc"
      run_params:
        expected_population: "AFR"
        sd_threshold: 6

Interpretation Guide
--------------------

Sample QC Metrics
~~~~~~~~~~~~~~~~~

**Call Rate**

- Normal: >98% (most samples should have >99%)
- Problematic: <98% indicates genotyping failure
- Action: Remove samples <98% call rate

**Heterozygosity**

- Normal: Within ±3 SD of cohort mean
- High: May indicate DNA contamination
- Low: May indicate inbreeding or DNA quality issues
- Action: Remove outliers beyond ±3 SD

**Sex Discordance**

- F statistic: <0.2 (XX/female), >0.8 (XY/male), 0.2-0.8 (ambiguous)
- Discordance: F statistic doesn't match reported sex
- Action: Remove discordant samples or verify phenotype file

**Relatedness**

- PI_HAT: 0.5 (parent-offspring/full sibling), 0.25 (half-sibling), 0.125 (first cousin)
- Threshold: 0.1875 (between 2nd and 3rd degree)
- Action: Remove one sample from each related pair (keep higher call rate)

Ancestry QC Metrics
~~~~~~~~~~~~~~~~~~~

**PCA Interpretation**

- PC1 typically separates African vs. non-African
- PC2 typically separates European vs. East Asian
- Outliers: >6 SD from expected population centroid
- Visual inspection recommended for borderline cases

**Population Stratification**

- Presence of clusters suggests population substructure
- May need to include PCs as covariates in GWAS
- Consider ancestry-stratified analysis for multi-ethnic cohorts

Variant QC Metrics
~~~~~~~~~~~~~~~~~~

**Hardy-Weinberg Equilibrium**

- Controls: Strict threshold (p<1e-6) - genotyping errors
- Cases: Lenient threshold (p<1e-10) - allow disease-associated variants
- Rationale: Disease SNPs may violate HWE in cases due to selection

**Differential Missingness**

- Tests: Is missingness rate different between cases and controls?
- Significance: p<1e-5 indicates batch effects or technical artifacts
- Action: Remove variants with significant differential missingness

**MAF Filtering**

- Threshold: Typically 1% for GWAS, higher (5%) for rare variant analysis
- Rationale: Rare variants have low statistical power and higher error rates
- Consider: Keep all variants for imputation reference panels

Best Practices
--------------

Workflow Recommendations
~~~~~~~~~~~~~~~~~~~~~~~~

1. **Pre-QC Checks**
   
   - Verify input file integrity (check .bim, .bed, .fam consistency)
   - Ensure phenotypes are correctly coded (1=control, 2=case)
   - Review number of samples and variants before QC

2. **QC Order Matters**
   
   - Always run: Sample QC → Ancestry QC → Variant QC
   - Rationale: Poor samples bias variant statistics; ancestry outliers bias association tests

3. **Iterative QC**
   
   - Run QC, review plots, adjust thresholds if needed
   - Document all threshold changes and rationale
   - Rerun with adjusted parameters if necessary

4. **Documentation**
   
   - Save all QC plots and summary statistics
   - Record number of samples/variants removed at each step
   - Document rationale for any non-standard thresholds

Troubleshooting
~~~~~~~~~~~~~~~

**High Sample Failure Rate (>10%)**

- Check: DNA quality, genotyping platform issues
- Review: Call rate distributions, batch effects
- Action: May need to adjust thresholds or exclude problematic batches

**Many Ancestry Outliers (>5%)**

- Check: Is expected_population correct?
- Review: PCA plots - is there substructure?
- Action: Consider ancestry-stratified QC or adjust sd_threshold

**Excessive Variant Removal (>50%)**

- Check: Are MAF and HWE thresholds too strict?
- Review: Variant call rate distributions
- Action: May indicate genotyping quality issues

**No Relatedness Detected (Unexpected)**

- Check: Are kinship calculations enabled (run_ibd: true)?
- Review: Are there family structures in study design?
- Action: Verify input data includes expected relatives

Output Files
------------

Each QC step produces standardized outputs:

**Sample QC Outputs:**

- ``mydata_SampleQCed.bed/bim/fam`` - Cleaned PLINK files
- ``sample_qc_report.txt`` - Summary statistics
- ``missingness_plot.png`` - Sample call rate distribution
- ``heterozygosity_plot.png`` - F coefficient distribution
- ``sex_check_plot.png`` - X chromosome homozygosity
- ``kinship_heatmap.png`` - Relatedness matrix
- ``removed_samples.txt`` - List of excluded samples with reasons

**Ancestry QC Outputs:**

- ``mydata_AncestryQCed.bed/bim/fam`` - Cleaned PLINK files
- ``pca_plot_PC1_PC2.png`` - Primary PCs with reference
- ``pca_plot_PC3_PC4.png`` - Secondary PCs
- ``ancestry_outliers.txt`` - List of outlier samples
- ``pca_coordinates.txt`` - PC scores for all samples

**Variant QC Outputs:**

- ``mydata_variantQCed.bed/bim/fam`` - Analysis-ready PLINK files
- ``variant_qc_summary.txt`` - Filtering statistics
- ``hwe_case_plot.png`` - HWE distribution in cases
- ``hwe_control_plot.png`` - HWE distribution in controls
- ``maf_distribution.png`` - Allele frequency spectrum
- ``removed_variants.txt`` - List of excluded variants with reasons

Next Steps
----------

After completing QC, the clean data is ready for:

1. **GWAS Analysis**
   
   .. code-block:: bash
   
       # Use QC output as GWAS input
       ideal-genom run --config gwas_config.yaml

2. **Imputation**
   
   - Upload to Michigan/TopMed imputation servers
   - Use ``mydata_variantQCed`` files for pre-imputation QC

3. **Population Structure Analysis**
   
   .. code-block:: bash
   
       # Compute FST statistics, run PCA/UMAP
       ideal-genom run --config population_config.yaml


See Also
--------

- :doc:`configuration` - Detailed YAML configuration reference
- :doc:`SampleQC` - Sample QC API documentation
- :doc:`AncestryQC` - Ancestry QC API documentation  
- :doc:`VariantQC` - Variant QC API documentation
- :doc:`examples` - Example configurations and workflows
- :doc:`gwas_pipeline` - GWAS pipeline documentation
