GWAS Pipeline
=============

The GWAS (Genome-Wide Association Study) pipeline in IDEAL-GENOM provides comprehensive tools for performing association analysis between genetic variants and phenotypes. The pipeline includes preparatory steps, statistical analysis using both fixed and mixed effects models, and gene annotation of significant findings.

Overview
--------

The GWAS pipeline consists of three main components:

1. **Preparatory Analysis**: LD pruning and PCA decomposition to prepare data
2. **GLM Analysis**: Fixed effects association testing using generalized linear models
3. **GLMM Analysis**: Mixed model analysis accounting for relatedness and population structure

**Key Features:**

- Automated LD pruning with high-LD region filtering
- Principal component analysis for population stratification
- Fixed effects (GLM) and random effects (GLMM) association testing
- Independent signal identification using GCTA COJO
- Automatic gene annotation with Ensembl or RefSeq
- Support for binary (case-control) phenotypes
- Resource-aware parallel processing

Prerequisites
-------------

**Required Software:**

- PLINK 2.0 (for GLM analysis)
- GCTA (for GLMM analysis and COJO)
- Python 3.11+ with required packages

**Required Data:**

- Quality-controlled PLINK binary files (.bed, .bim, .fam)
- Phenotype data in .fam file (case=2, control=1, missing=0/-9)
- High-LD regions file (auto-downloaded if not provided)

**Recommended Preprocessing:**

Run the QC pipeline first to ensure data quality:

.. code-block:: bash

    # Run QC pipeline first
    ideal-genom run --config qc_pipeline.yaml
    
    # Then run GWAS on clean data
    ideal-genom run --config gwas_pipeline.yaml

Quick Start
-----------

**1. Get the GWAS Configuration Template**

.. code-block:: bash

    # Copy the GWAS template
    cp yaml_configs/gwas_config_template.yaml my_gwas.yaml

**2. Edit Configuration**

Update paths to your QC-cleaned data:

.. code-block:: yaml

    pipeline:
      name: "my_gwas"
      base_output_dir: "/data/gwas_output"
      steps:
        - name: "preparatory"
          enabled: true
          init_params:
            input_path: "/data/qc_output/variant_qc/clean_files"
            input_name: "mydata_variantQCed"
            output_path: "${base_output_dir}"
            build: "38"

**3. Run the Pipeline**

.. code-block:: bash

    ideal-genom run --config my_gwas.yaml

Pipeline Steps
--------------

Step 1: Preparatory Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Prepares genetic data for association testing by performing LD pruning and PCA decomposition.

**What it does:**

1. Filters high-LD regions
2. Applies QC filters (MAF, missingness, HWE)
3. Performs LD-based pruning to identify independent SNPs
4. Computes principal components for population structure

**Configuration:**

.. code-block:: yaml

    - name: "preparatory"
      enabled: true
      module: "ideal_genom.gwas.preparatory"
      class: "Preparatory"
      init_params:
        input_path: "/data/qc/clean_files"
        input_name: "mydata_clean"
        output_path: "${base_output_dir}"
        output_name: "gwas_prep"
        high_ld_regions_file: "auto"
        build: "38"
      execute_params:
        # QC filters
        mind: 0.1                   # Max missing per individual
        maf: 0.01                   # Minor allele frequency
        geno: 0.1                   # Max missing per SNP
        hwe: 5.0e-6                 # Hardy-Weinberg p-value
        
        # LD pruning
        ind_pair: [50, 5, 0.2]      # Window, step, r² threshold
        
        # PCA
        pca: 10                     # Number of PCs to compute
        
        # Resources
        memory: null                # Auto-detect
        threads: null               # Auto-detect

**Parameters Explained:**

- **mind** (float, 0-1): Maximum missing genotype rate per individual
- **maf** (float, 0-0.5): Minor allele frequency threshold for filtering
- **geno** (float, 0-1): Maximum missing genotype rate per SNP
- **hwe** (float, 0-1): Hardy-Weinberg equilibrium p-value threshold
- **ind_pair** (list): LD pruning [window_size, step_size, r²_threshold]
  
  - window_size: Number of variants in sliding window (default: 50)
  - step_size: Window shift size in variants (default: 5)
  - r² threshold: Prune variants with r² > threshold (default: 0.2)

- **pca** (int): Number of principal components to compute (typically 10-20)
- **memory** (int or null): Memory in MB (null = auto-detect 2/3 available RAM)
- **threads** (int or null): CPU threads (null = auto-detect cores - 2)

**Output Files:**

.. code-block:: text

    preparatory/
    ├── gwas_prep-prunning.bed/bim/fam    # After QC filters
    ├── gwas_prep-prunning.prune.in       # SNPs passing LD pruning
    ├── gwas_prep-prunning.prune.out      # SNPs removed by LD pruning
    ├── gwas_prep-pruned.bed/bim/fam      # Pruned dataset
    ├── mydata_clean.eigenvec             # PC scores
    └── mydata_clean.eigenval             # Eigenvalues (variance explained)

Step 2: GLM Analysis
^^^^^^^^^^^^^^^^^^^^

Performs fixed effects association testing using generalized linear models with PCA covariates.

**What it does:**

1. Runs logistic regression with PCA covariates
2. Identifies genome-wide significant variants (p < 5×10⁻⁸)
3. Uses GCTA COJO to find independent signals
4. Annotates significant variants with gene information

**Configuration:**

.. code-block:: yaml

    - name: "gwas_glm"
      enabled: true
      module: "ideal_genom.gwas.gen_linear_model"
      class: "GWAS_GLM"
      init_params:
        input_path: "${steps.preparatory.input_path}"
        input_name: "${steps.preparatory.input_name}"
        output_path: "${base_output_dir}"
        output_name: "gwas_glm_results"
        recompute: true
      execute_params:
        # Association testing filters
        maf: 0.01                   # MAF threshold
        mind: 0.1                   # Missing per individual
        hwe: 5.0e-6                 # HWE p-value
        ci: 0.95                    # Confidence interval
        
        # Annotation
        gtf_path: null              # Custom GTF (null = auto-download)
        build: "38"                 # Genome build
        anno_source: "ensembl"      # "ensembl" or "refseq"

**Parameters Explained:**

- **maf** (float): MAF threshold for variants included in analysis
- **mind** (float): Maximum missing rate per individual
- **hwe** (float): Hardy-Weinberg equilibrium p-value threshold
- **ci** (float): Confidence interval for effect size estimates (0-1)
- **gtf_path** (string or null): Path to custom GTF annotation file
- **build** (string): "37" for GRCh37/hg19, "38" for GRCh38/hg38
- **anno_source** (string): Annotation source - "ensembl" or "refseq"
- **recompute** (bool): If false, skip if results already exist

**Output Files:**

.. code-block:: text

    gwas_glm/
    ├── gwas_glm_results_glm.PHENO1.glm.logistic.hybrid
    │   # Full GWAS results with all tested variants
    ├── gwas_glm_results_glm.PHENO1.glm.logistic.hybrid.adjusted
    │   # Results with adjusted p-values
    ├── cojo_file.ma
    │   # Prepared for COJO analysis
    ├── gwas_glm_results-cojo.jma.cojo
    │   # Independent genome-wide significant variants
    └── top_hits_annotated.tsv
        # Annotated top hits with gene names, positions, effects

Step 3: GLMM Analysis
^^^^^^^^^^^^^^^^^^^^^

Performs mixed model association testing accounting for population structure and cryptic relatedness.

**What it does:**

1. Computes genetic relationship matrix (GRM) from pruned SNPs
2. Creates sparse GRM for computational efficiency
3. Runs fastGWA mixed model with PCA and sex covariates
4. Identifies independent signals using GCTA COJO
5. Annotates significant variants with gene information

**Configuration:**

.. code-block:: yaml

    - name: "gwas_glmm"
      enabled: true
      module: "ideal_genom.gwas.gen_linear_mix_model"
      class: "GWAS_GLMM"
      init_params:
        input_path: "${steps.preparatory.input_path}"
        input_name: "${steps.preparatory.input_name}"
        output_path: "${base_output_dir}"
        output_name: "gwas_glmm_results"
        recompute: true
      execute_params:
        # Association parameters
        maf: 0.01
        
        # GRM computation
        pruned_file: "${steps.preparatory.pruned_file}"
        max_threads: null           # Auto-detect
        
        # Annotation
        gtf_path: null
        build: "38"
        anno_source: "ensembl"

**Parameters Explained:**

- **maf** (float): MAF threshold for variants in analysis
- **pruned_file** (string): Path to pruned PLINK files (without extension)
- **max_threads** (int or null): Threads for GRM computation
- Other parameters same as GLM

**Output Files:**

.. code-block:: text

    gwas_glmm/
    ├── mydata_clean_pheno.phen
    │   # Phenotype file (FID, IID, phenotype)
    ├── mydata_clean_sex.covar
    │   # Sex covariate file
    ├── mydata_clean_grm.grm.bin
    │   # Genetic relationship matrix (binary)
    ├── mydata_clean_grm.grm.id
    │   # Sample IDs for GRM
    ├── mydata_clean_grm.grm.N.bin
    │   # Number of variants used per pair
    ├── mydata_clean_sparse.grm.sp
    │   # Sparse GRM (efficient storage)
    ├── mydata_clean_sparse.grm.id
    │   # Sample IDs for sparse GRM
    ├── gwas_glmm_results_assocSparseCovar_pca_sex-mlm-binary.fastGWA
    │   # Full GWAS results
    ├── cojo_file.ma
    │   # Prepared for COJO analysis
    ├── gwas_glmm_results-cojo.jma.cojo
    │   # Independent genome-wide significant variants
    └── top_hits_annotated.tsv
        # Annotated top hits

Complete Workflow Example
--------------------------

**Full GWAS Configuration:**

.. code-block:: yaml

    pipeline:
      name: "my_gwas_study"
      base_output_dir: "/data/gwas_results"
      
      steps:
        # Step 1: Prepare data
        - name: "preparatory"
          enabled: true
          module: "ideal_genom.gwas.preparatory"
          class: "Preparatory"
          init_params:
            input_path: "/data/qc/clean_files"
            input_name: "study_clean"
            output_path: "${base_output_dir}"
            output_name: "gwas_prep"
            high_ld_regions_file: "auto"
            build: "38"
          execute_params:
            mind: 0.02
            maf: 0.01
            geno: 0.02
            hwe: 1.0e-6
            ind_pair: [50, 5, 0.2]
            pca: 10
        
        # Step 2: Fixed effects analysis
        - name: "gwas_glm"
          enabled: true
          module: "ideal_genom.gwas.gen_linear_model"
          class: "GWAS_GLM"
          init_params:
            input_path: "${steps.preparatory.input_path}"
            input_name: "${steps.preparatory.input_name}"
            output_path: "${base_output_dir}"
            output_name: "gwas_glm"
            recompute: false
          execute_params:
            maf: 0.01
            mind: 0.02
            hwe: 1.0e-6
            ci: 0.95
            build: "38"
            anno_source: "ensembl"
        
        # Step 3: Mixed model analysis (optional)
        - name: "gwas_glmm"
          enabled: false  # Enable if needed
          module: "ideal_genom.gwas.gen_linear_mix_model"
          class: "GWAS_GLMM"
          init_params:
            input_path: "${steps.preparatory.input_path}"
            input_name: "${steps.preparatory.input_name}"
            output_path: "${base_output_dir}"
            output_name: "gwas_glmm"
            recompute: false
          execute_params:
            maf: 0.01
            pruned_file: "${steps.preparatory.pruned_file}"
            build: "38"
            anno_source: "ensembl"
    
    settings:
      logging:
        level: "INFO"
        file_logging: true
      resources:
        max_memory: null
        max_threads: null
      files:
        keep_intermediate: true

**Running the Analysis:**

.. code-block:: bash

    # Validate configuration
    ideal-genom validate --config my_gwas.yaml
    
    # Dry run to preview
    ideal-genom run --config my_gwas.yaml --dry-run
    
    # Execute pipeline
    ideal-genom run --config my_gwas.yaml

Parameter Recommendations
^^^^^^^^^^^^^^^^^^^^^^^^^

**Conservative Analysis (reduce false positives):**

.. code-block:: yaml

    execute_params:
      maf: 0.05      # Higher MAF
      hwe: 1.0e-10   # Stricter HWE
      mind: 0.01     # Stricter missingness
      geno: 0.01

**Liberal Analysis (increase power):**

.. code-block:: yaml

    execute_params:
      maf: 0.01      # Lower MAF
      hwe: 1.0e-4    # Relaxed HWE
      mind: 0.05
      geno: 0.05

**Rare Variant Analysis:**

.. code-block:: yaml

    execute_params:
      maf: 0.001     # Include rare variants
      geno: 0.02     # Stricter missingness for rare variants
      hwe: 1.0e-6

Resource Management
^^^^^^^^^^^^^^^^^^^

**For Large Datasets (> 500K SNPs, > 10K samples):**

.. code-block:: yaml

    settings:
      resources:
        max_memory: 64000    # 64GB
        max_threads: 16      # Use many cores
      files:
        keep_intermediate: false  # Save disk space

**For Small Datasets:**

.. code-block:: yaml

    settings:
      resources:
        max_memory: 16000    # 16GB sufficient
        max_threads: 4
      files:
        keep_intermediate: true  # Keep for inspection

Troubleshooting
---------------

Common Issues and Solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Issue: "No genome-wide significant hits"**

Solutions:

1. Check sample size (need sufficient power)
2. Verify phenotype coding (cases=2, controls=1)
3. Review QC stringency (may be too strict)
4. Check for population stratification
5. Consider suggestive hits (p < 1×10⁻⁵)

**Issue: "Inflation of test statistics (λ > 1.1)"**

Solutions:

1. Check for population stratification
2. Increase number of PCs used as covariates
3. Consider using GLMM instead of GLM
4. Review sample quality (duplicates, relatedness)

**Issue: "GRM computation fails or runs out of memory"**

Solutions:

1. Reduce max_threads (more memory per thread)
2. Increase max_memory setting
3. Ensure pruned file has reasonable number of SNPs (50-100K)
4. Check available system memory

**Issue: "COJO analysis produces no results"**

Solutions:

1. Verify you have genome-wide significant hits
2. Check reference LD panel is appropriate
3. Ensure sufficient sample size
4. Review MAF threshold (not too high)

**Issue: "Annotation fails"**

Solutions:

1. Check internet connection (for auto-download)
2. Provide custom GTF file with gtf_path
3. Verify genome build matches your data
4. Check gene database is accessible

Performance Optimization
^^^^^^^^^^^^^^^^^^^^^^^^

**Speed up analysis:**

1. Use GLM instead of GLMM when possible
2. Set ``recompute: false`` for completed steps
3. Reduce number of PCs if population homogeneous
4. Use ``keep_intermediate: false`` to save I/O
5. Allocate more threads for parallel processing

**Reduce memory usage:**

1. Reduce number of threads (more memory per thread)
2. Use sparse GRM for GLMM
3. Process chromosomes separately if needed
4. Close other applications during analysis

See Also
--------

- :doc:`configuration` - Full parameter reference
- :doc:`qc_pipeline` - Run QC before GWAS
- :doc:`examples` - Complete workflow examples
- :doc:`troubleshooting` - Detailed problem-solving guide

Additional Resources
--------------------

**PLINK Documentation:**
- PLINK 2.0: https://www.cog-genomics.org/plink/2.0/
- Logistic regression: https://www.cog-genomics.org/plink/2.0/assoc

**GCTA Documentation:**
- GCTA overview: https://yanglab.westlake.edu.cn/software/gcta/
- fastGWA: https://yanglab.westlake.edu.cn/software/gcta/#fastGWA
- COJO: https://yanglab.westlake.edu.cn/software/gcta/#COJO
