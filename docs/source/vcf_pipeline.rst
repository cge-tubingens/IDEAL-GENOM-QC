VCF Processing Pipeline
========================

The VCF (Variant Call Format) processing pipeline in IDEAL-GENOM provides tools for handling post-imputation VCF files. It automates the complete workflow from imputed VCF files to analysis-ready PLINK binary format, including quality filtering, normalization, annotation, and format conversion.

Overview
--------

The VCF pipeline consists of two main components:

1. **VCF Processing**: Filter, normalize, annotate, and concatenate imputed VCF files
2. **PLINK Conversion**: Convert processed VCF files to PLINK binary format

**Key Features:**

- Parallel processing of multiple chromosome VCF files
- Quality filtering based on imputation R² score
- Reference genome normalization
- Variant ID annotation with dbSNP
- Automatic handling of multiallelic variants
- Efficient conversion to PLINK format
- Resource-aware parallel execution

Use Cases
---------

This pipeline is ideal for:

- **Post-Imputation Processing**: Clean and prepare imputed genotype data
- **TOPMed/HRC Imputation**: Process results from imputation servers
- **Data Harmonization**: Standardize variant representation and IDs
- **GWAS Preparation**: Convert VCF to PLINK for association analysis

Prerequisites
-------------

**Required Software:**

- BCFtools (v1.10+) for VCF manipulation
- PLINK 2.0 for format conversion
- Python 3.11+ with required packages

**Required Data:**

- Imputed VCF files (can be gzipped)
- Reference genome FASTA (auto-downloaded if not provided)
- dbSNP VCF for annotation (optional but recommended)

**Input File Formats:**

Supported input formats:

- ``.vcf.gz`` - Gzipped VCF files (recommended)
- ``.vcf`` - Uncompressed VCF files
- ``.zip`` - Zipped VCF files (automatically extracted)

Quick Start
-----------

**1. Get the VCF Configuration Template**

.. code-block:: bash

    # Copy the VCF processing template
    cp yaml_configs/vcf_config_template.yaml my_vcf_pipeline.yaml

**2. Edit Configuration**

Update paths to your imputed data:

.. code-block:: yaml

    pipeline:
      name: "vcf_processing"
      base_output_dir: "/data/output"
      steps:
        - name: "process_vcf"
          enabled: true
          init_params:
            input_path: "/data/imputed_vcfs"
            output_path: "${base_output_dir}"
            output_name: "processed_data.vcf.gz"
          execute_params:
            r2_threshold: 0.3
            build: "38"

**3. Run the Pipeline**

.. code-block:: bash

    ideal-genom run --config my_vcf_pipeline.yaml

Pipeline Steps
--------------

Step 1: VCF Processing
^^^^^^^^^^^^^^^^^^^^^^^

Processes imputed VCF files through multiple quality control and normalization steps.

**What it does:**

1. **Unzip Files**: Extract compressed VCF files if needed
2. **Filter by R²**: Remove poorly imputed variants (R² < threshold)
3. **Normalize**: Split multiallelic variants and left-align indels
4. **Reference Normalization**: Ensure variants match reference genome
5. **Annotate IDs**: Add or update variant IDs with dbSNP rs numbers
6. **Concatenate**: Merge all chromosome VCFs into a single file
7. **Index**: Create tabix index for the final VCF

**Configuration:**

.. code-block:: yaml

    - name: "process_vcf"
      enabled: true
      module: "ideal_genom.post_imputation.vcf_process"
      class: "ProcessVCF"
      init_params:
        input_path: "/data/imputed_vcfs"
        output_path: "${base_output_dir}"
        input_name: "placeholder"
        output_name: "processed_data.vcf.gz"
      execute_params:
        password: null              # For password-protected zip files
        r2_threshold: 0.3           # Minimum imputation quality
        build: "38"                 # Genome build ("37" or "38")
        ref_genome: null            # Custom reference (auto-downloaded if null)
        ref_annotation: null        # dbSNP VCF for annotation
        max_threads: null           # Thread count (auto-detect if null)

**Parameters Explained:**

**init_params:**

- **input_path** (string): Directory containing imputed VCF files

  - Can contain ``.vcf.gz``, ``.vcf``, or ``.zip`` files
  - Automatically detects chromosome-specific files

- **output_path** (string): Where to save processed files
- **input_name** (string): Placeholder for compatibility (can be any value)
- **output_name** (string): Name for final concatenated VCF file

**execute_params:**

- **password** (string or null): Password for encrypted zip files
- **r2_threshold** (float, 0-1): Minimum R² (imputation quality) score

  - 0.3: Liberal (includes more variants, some lower quality)
  - 0.5: Moderate (balance quality and quantity)
  - 0.8: Conservative (high quality, fewer variants)
  - Variants with R² < threshold are removed

- **build** (string): Genome build version

  - "37": GRCh37/hg19
  - "38": GRCh38/hg38

- **ref_genome** (string or null): Path to reference genome FASTA

  - ``null``: Auto-download from Ensembl based on build
  - Custom path: Use your own reference FASTA file

- **ref_annotation** (string or null): Path to dbSNP VCF for variant annotation

  - ``null``: Skip annotation step
  - Path to dbSNP VCF: Annotate variants with rs IDs

- **max_threads** (int or null): Maximum threads for parallel processing

  - ``null``: Auto-detect (uses available cores)
  - Specific number: Limit thread usage

**Processing Pipeline Details:**

.. code-block:: text

    Input VCF files (per chromosome)
           ↓
    [1. Unzip if needed]
           ↓
    [2. Filter: R² >= threshold]
           ↓
    [3. Normalize: split multiallelic, left-align]
           ↓
    [4. Reference normalization]
           ↓
    [5. Annotate with dbSNP IDs]
           ↓
    [6. Concatenate all chromosomes]
           ↓
    [7. Index final VCF]
           ↓
    Output: processed_data.vcf.gz + .tbi

**Output Files:**

.. code-block:: text

    process_vcf/
    ├── unzipped/                   # Extracted VCF files
    │   ├── chr1.vcf.gz
    │   ├── chr2.vcf.gz
    │   └── ...
    ├── filtered/                   # After R² filtering
    │   ├── chr1_filtered.vcf.gz
    │   ├── chr2_filtered.vcf.gz
    │   └── ...
    ├── normalized/                 # After normalization
    │   ├── chr1_norm.vcf.gz
    │   ├── chr2_norm.vcf.gz
    │   └── ...
    ├── ref_normalized/             # After reference normalization
    │   ├── chr1_refnorm.vcf.gz
    │   └── ...
    ├── annotated/                  # After dbSNP annotation
    │   ├── chr1_annotated.vcf.gz
    │   └── ...
    └── processed_data.vcf.gz       # Final concatenated file
        └── processed_data.vcf.gz.tbi  # Tabix index

Step 2: PLINK Conversion
^^^^^^^^^^^^^^^^^^^^^^^^^

Converts processed VCF files to PLINK binary format for GWAS and other analyses.

**What it does:**

1. Converts VCF to PLINK1.9 binary format (.bed/.bim/.fam)
2. Handles sample ID formatting
3. Updates family information if provided
4. Optimizes for downstream analysis

**Configuration:**

.. code-block:: yaml

    - name: "plink_conversion"
      enabled: true
      module: "ideal_genom.post_imputation.vcf_to_plink"
      class: "GetPLINK"
      init_params:
        input_path: "${steps.process_vcf.output_path}"
        input_name: "${steps.process_vcf.concatenated_file}"
        output_path: "${base_output_dir}"
        output_name: "imputed_data_plink"
      execute_params:
        double_id: true             # Use sample ID for both FID and IID
        for_fam_update_file: null   # FAM file for updating family info
        threads: null               # Thread count (auto-detect)
        memory: null                # Memory in MB (auto-detect)

**Parameters Explained:**

**init_params:**

- **input_path** (string): Directory containing processed VCF
- **input_name** (string): Name of input VCF file (without extension)
- **output_path** (string): Where to save PLINK files
- **output_name** (string): Prefix for output PLINK files

**execute_params:**

- **double_id** (bool): Sample ID handling

  - ``true``: Set FID = IID (recommended for most cases)
  - ``false``: Use VCF sample IDs as-is

- **for_fam_update_file** (string or null): Path to FAM file for updating

  - ``null``: No update, use VCF phenotype information
  - Path: Update .fam with information from this file (FID, IID, phenotype, sex)

- **threads** (int or null): Number of threads for conversion
- **memory** (int or null): Memory allocation in MB

**Output Files:**

.. code-block:: text

    plink_conversion/
    ├── imputed_data_plink.bed      # Binary genotype file
    ├── imputed_data_plink.bim      # Variant information
    ├── imputed_data_plink.fam      # Sample information
    └── imputed_data_plink.log      # Conversion log

Complete Workflow Example
--------------------------

**Full VCF Processing Configuration:**

.. code-block:: yaml

    pipeline:
      name: "imputed_data_processing"
      base_output_dir: "/data/processed_imputation"
      
      steps:
        # Step 1: Process VCF files
        - name: "process_vcf"
          enabled: true
          module: "ideal_genom.post_imputation.vcf_process"
          class: "ProcessVCF"
          init_params:
            input_path: "/data/imputation_results"
            output_path: "${base_output_dir}"
            input_name: "placeholder"
            output_name: "imputed_clean.vcf.gz"
          execute_params:
            password: null
            r2_threshold: 0.3
            build: "38"
            ref_genome: null
            ref_annotation: "/data/references/dbSNP156_GRCh38.vcf.gz"
            max_threads: null
        
        # Step 2: Convert to PLINK
        - name: "plink_conversion"
          enabled: true
          module: "ideal_genom.post_imputation.vcf_to_plink"
          class: "GetPLINK"
          init_params:
            input_path: "${steps.process_vcf.output_path}"
            input_name: "${steps.process_vcf.concatenated_file}"
            output_path: "${base_output_dir}"
            output_name: "imputed_plink"
          execute_params:
            double_id: true
            for_fam_update_file: "/data/original_study.fam"
            threads: null
            memory: null
    
    settings:
      logging:
        level: "INFO"
        file_logging: true
      resources:
        max_memory: null
        max_threads: null
      files:
        keep_intermediate: true

**Running the Pipeline:**

.. code-block:: bash

    # Validate
    ideal-genom validate --config vcf_pipeline.yaml
    
    # Execute
    ideal-genom run --config vcf_pipeline.yaml

Input Data Organization
-----------------------

Expected Directory Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For chromosome-specific VCF files:

.. code-block:: text

    /data/imputation_results/
    ├── chr1.dose.vcf.gz
    ├── chr2.dose.vcf.gz
    ├── chr3.dose.vcf.gz
    ├── ...
    ├── chr22.dose.vcf.gz
    └── chrX.dose.vcf.gz

Or as zip files:

.. code-block:: text

    /data/imputation_results/
    ├── chr1.zip
    ├── chr2.zip
    ├── ...
    └── chr22.zip

**File Naming Patterns:**

The pipeline automatically detects files with these patterns:

- ``chr*.vcf.gz``
- ``chr*.vcf``
- ``*.dose.vcf.gz``
- ``chromosome_*.vcf.gz``
- ``chr*.zip``

VCF File Requirements
^^^^^^^^^^^^^^^^^^^^^

**Required INFO Fields:**

- **R2** or **INFO**: Imputation quality score (R²)

  - Used for quality filtering
  - Typically ranges from 0 to 1

Quality Control Considerations
-------------------------------

R² Threshold Selection
^^^^^^^^^^^^^^^^^^^^^^

**Recommended Thresholds:**

.. code-block:: text

    R² Threshold  |  Quality     |  Use Case
    --------------|--------------|--------------------------------
    0.3           |  Liberal     |  Discovery, maximize variants
    0.5           |  Moderate    |  Standard analysis
    0.8           |  Conservative|  High-confidence variants only
    0.9           |  Very strict |  Fine-mapping, functional studies

**Trade-offs:**

- **Lower threshold (0.3)**: More variants, some lower quality
- **Higher threshold (0.8)**: Fewer variants, higher confidence

**Check Imputation Quality:**

.. code-block:: bash

    # Before filtering
    bcftools query -f '%INFO/R2\n' chr1.vcf.gz | \\
        awk '{sum+=$1; count++} END {print "Mean R2:", sum/count}'
    
    # Distribution
    bcftools query -f '%INFO/R2\n' chr1.vcf.gz | \\
        awk '{if($1<0.3) low++; else if($1<0.8) med++; else high++} 
             END {print "R2<0.3:", low, "0.3-0.8:", med, "R2>=0.8:", high}'

Reference Genome Selection
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Auto-Download (Recommended):**

Set ``ref_genome: null`` and specify build:

- Build "37": Downloads GRCh37 from Ensembl
- Build "38": Downloads GRCh38 from Ensembl

**Custom Reference:**

Use your own reference if:

- Using non-standard reference assembly
- Working offline
- Need specific reference version

Ensure reference matches imputation reference!

dbSNP Annotation
^^^^^^^^^^^^^^^^

**Benefits of Annotation:**

- Standardized variant IDs (rs numbers)
- Easier result interpretation
- Facilitates meta-analysis
- Improves result sharing

**Download dbSNP:**

.. code-block:: bash

    # GRCh38
    wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
    wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi
    
    # GRCh37
    wget https://ftp.ncbi.nih.gov/snp/archive/b151/VCF/GCF_000001405.25.gz

Best Practices
--------------

Pre-Processing Checklist
^^^^^^^^^^^^^^^^^^^^^^^^^

Before running the VCF pipeline:

1. **Verify File Integrity:**

   .. code-block:: bash

       # Check VCF files
       for vcf in chr*.vcf.gz; do
           bcftools view -H $vcf | head -1 || echo "Error in $vcf"
       done

2. **Check Sample Counts:**

   .. code-block:: bash

       # Ensure consistent sample numbers
       for vcf in chr*.vcf.gz; do
           echo "$vcf: $(bcftools query -l $vcf | wc -l) samples"
       done

3. **Inspect Imputation Quality:**

   .. code-block:: bash

       # Check R² distribution
       bcftools query -f '%INFO/R2\n' chr1.vcf.gz | \\
           python -c "import sys; import numpy as np; \\
           r2 = [float(x) for x in sys.stdin]; \\
           print(f'Mean: {np.mean(r2):.3f}'); \\
           print(f'Median: {np.median(r2):.3f}')"

Workflow Optimization
^^^^^^^^^^^^^^^^^^^^^

**For Large Datasets (>500K variants per chromosome):**

.. code-block:: yaml

    execute_params:
      r2_threshold: 0.5          # More stringent
      max_threads: 16            # Use more cores
    
    settings:
      resources:
        max_memory: 64000        # 64GB
        max_threads: 16
      files:
        keep_intermediate: false  # Save disk space

**For Small Datasets:**

.. code-block:: yaml

    execute_params:
      r2_threshold: 0.3          # More lenient
      max_threads: 4
    
    settings:
      files:
        keep_intermediate: true  # Keep for inspection

**For Repeated Analysis:**

Set ``keep_intermediate: true`` and ``recompute: false`` in init_params to avoid reprocessing.

Common Workflows
----------------

Workflow 1: Basic VCF Processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Process and convert imputed VCF files:

.. code-block:: yaml

    pipeline:
      steps:
        - name: "process_vcf"
          enabled: true
          execute_params:
            r2_threshold: 0.3
            build: "38"
        
        - name: "plink_conversion"
          enabled: true
          execute_params:
            double_id: true

Workflow 2: High-Quality Variant Subset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Extract only high-confidence variants:

.. code-block:: yaml

    pipeline:
      steps:
        - name: "process_vcf"
          enabled: true
          execute_params:
            r2_threshold: 0.8      # Strict filtering
            build: "38"
            ref_annotation: "/data/dbSNP.vcf.gz"

Workflow 3: Update Phenotype Information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Merge with existing phenotype data:

.. code-block:: yaml

    pipeline:
      steps:
        - name: "process_vcf"
          enabled: true
        
        - name: "plink_conversion"
          enabled: true
          execute_params:
            for_fam_update_file: "/data/phenotypes.fam"
            double_id: true

Troubleshooting
---------------

Common Issues
^^^^^^^^^^^^^

**Issue: "R2 field not found in VCF"**

Solutions:

1. Check INFO field name (might be INFO, IMP, IMPINFO)
2. Verify VCF format with ``bcftools view -h file.vcf.gz``
3. Update code to match your VCF's R² field name

**Issue: "Reference genome download fails"**

Solutions:

1. Check internet connection
2. Download manually and provide path with ``ref_genome``
3. Use alternative Ensembl mirror

**Issue: "Memory error during concatenation"**

Solutions:

1. Set ``max_memory`` explicitly in settings
2. Reduce ``max_threads`` (more memory per thread)
3. Process subsets of chromosomes separately
4. Close other applications

**Issue: "Sample IDs don't match phenotype file"**

Solutions:

1. Verify sample ID format in VCF
2. Check FAM file sample ID format (FID and IID)
3. Use ``double_id: true`` if VCF has only one ID per sample
4. Manually create matching FAM file

**Issue: "Missing variants after filtering"**

Solutions:

1. Lower ``r2_threshold`` (e.g., from 0.8 to 0.3)
2. Check imputation quality of input files
3. Verify R² field is correctly parsed

Performance Issues
^^^^^^^^^^^^^^^^^^

**Slow processing:**

1. Increase ``max_threads`` in execute_params
2. Use SSD storage for temporary files
3. Reduce ``r2_threshold`` to filter earlier
4. Process chromosomes in parallel manually

**High memory usage:**

1. Reduce ``max_threads`` (each thread uses memory)
2. Process chromosomes separately
3. Set explicit ``max_memory`` limit
4. Enable ``keep_intermediate: false``

Output Validation
-----------------

Verify Pipeline Success
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # Check file integrity
    plink2 --bfile imputed_plink --validate
    
    # Verify sample count
    wc -l imputed_plink.fam
    
    # Check variant info
    head imputed_plink.bim

Integration with QC Pipeline
-----------------------------

After VCF processing, run the QC pipeline:

.. code-block:: bash

    # Step 1: Process VCF
    ideal-genom run --config vcf_pipeline.yaml
    
    # Step 2: Run QC on converted PLINK files
    # Update qc_pipeline.yaml to use processed data
    ideal-genom run --config qc_pipeline.yaml

**QC Configuration Update:**

.. code-block:: yaml

    pipeline:
      steps:
        - name: "sample_qc"
          init_params:
            input_path: "/data/processed_imputation/plink_conversion"
            input_name: "imputed_plink"

See Also
--------

- :doc:`configuration` - Full parameter reference
- :doc:`qc_pipeline` - Quality control after VCF processing
- :doc:`gwas_pipeline` - Association analysis workflow
- :doc:`troubleshooting` - Detailed problem-solving

External Resources
------------------

**BCFtools Documentation:**
- Main page: http://samtools.github.io/bcftools/
- VCF filtering: http://samtools.github.io/bcftools/bcftools.html#view
- Normalization: http://samtools.github.io/bcftools/bcftools.html#norm

**PLINK Documentation:**
- PLINK 2.0: https://www.cog-genomics.org/plink/2.0/
- VCF import: https://www.cog-genomics.org/plink/2.0/input#vcf

**Imputation Services:**
- TOPMed: https://imputation.biodatacatalyst.nhlbi.nih.gov/
- Michigan Imputation Server: https://imputationserver.sph.umich.edu/
- HRC: http://www.haplotype-reference-consortium.org/

**VCF Format:**
- Specification: https://samtools.github.io/hts-specs/VCFv4.2.pdf
- Best practices: https://www.internationalgenome.org/wiki/Analysis/vcf4.0/
