Installation Guide
==================

System Requirements
-------------------

Before installing IDEAL-GENOM-QC, ensure you have the following prerequisites:

**Software Dependencies**
   - Python 3.8 or higher
   - PLINK 1.9 (required for QC operations)
   - PLINK 2.0 (required for QC operations)

**Hardware Requirements**
   - Minimum 8GB RAM (16GB recommended for large datasets)
   - At least 10GB free disk space
   - Multi-core processor recommended

Installing PLINK
-----------------

IDEAL-GENOM-QC requires both PLINK 1.9 and PLINK 2.0 to be installed and accessible in your system PATH.

**On Ubuntu/Debian:**

.. code-block:: bash

    # Install PLINK 1.9
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
    unzip plink_linux_x86_64_20231211.zip
    sudo mv plink /usr/local/bin/

    # Install PLINK 2.0
    wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20231212.zip
    unzip plink2_linux_x86_64_20231212.zip
    sudo mv plink2 /usr/local/bin/

**On macOS and Windows**

For PLINK1.9 and PLINK2.0 installation instructions, please refer to the official PLINK website: `https://www.cog-genomics.org/plink/1.9/` and `https://www.cog-genomics.org/plink/2.0/` respectively.

**Verify Installation:**

.. code-block:: bash

    plink --version
    plink2 --version

Installing IDEAL-GENOM-QC
--------------------------

Option 1: PyPI Installation (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the stable version from PyPI:

.. code-block:: bash

    pip install ideal-genom-qc

Option 2: Development Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the latest features and development version:

.. code-block:: bash

    git clone https://github.com/cge-tubingens/IDEAL-GENOM-QC.git
    cd IDEAL-GENOM-QC

Please, be aware that the development version may be unstable and it may have some bugs.

**Using Poetry (Recommended):**

.. code-block:: bash

    # Install Poetry if not already installed
    curl -sSL https://install.python-poetry.org | python3 -

    # Install dependencies
    poetry install

    # Activate the virtual environment
    poetry shell

**Using pip:**

.. code-block:: bash

    pip install -r requirements.txt
    pip install -e .

Option 3: Conda Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # Create a new conda environment
    conda create -n ideal-genom-qc python=3.9
    conda activate ideal-genom-qc
    
    # Install from PyPI
    pip install ideal-genom-qc

Docker Installation
-------------------

For containerized deployment:

.. code-block:: bash

    # Pull the Docker image
    docker pull cge-tubingens/ideal-genom-qc:latest

    # Or build from source
    git clone https://github.com/cge-tubingens/IDEAL-GENOM-QC.git
    cd IDEAL-GENOM-QC
    docker build -t ideal-genom-qc .

Verification
------------

Test your installation:

.. code-block:: python

    import ideal_genom_qc
    print(ideal_genom_qc.__version__)


Troubleshooting
---------------

**Common Issues:**

1. **PLINK not found**: Ensure PLINK is installed and in your PATH
2. **Permission errors**: Use `pip install --user` for user-only installation
3. **Poetry issues**: Update Poetry to the latest version (2.0+)

**Getting Help:**

- Check the :doc:`troubleshooting` guide
- Report issues on `GitHub <https://github.com/cge-tubingens/IDEAL-GENOM-QC/issues>`_
- Contact the development team

Next Steps
----------

After installation, proceed to the :doc:`getting_started` guide to learn how to use IDEAL-GENOM-QC.