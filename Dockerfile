# Use an official Python runtime as a parent image
FROM python:3.12-slim

# Set the working directory in the container
WORKDIR /app

# install system dependencies in the container
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Clone the specific repository
RUN git clone https://github.com/cge-tubingens/IDEAL-GENOM-QC .

# Install Python dependencies
RUN pip install --no-cache-dir .

# Copy the local code to the container
COPY . .

# Download and extract PLINK 1.9 and PLINK 2.0
RUN wget -O /tmp/plink19.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip \
    && unzip /tmp/plink19.zip -d /usr/local/bin/ \
    && wget -O /tmp/plink2.zip https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240105.zip \
    && unzip /tmp/plink2.zip -d /usr/local/bin/ \
    && chmod +x /usr/local/bin/plink*

# Download and install GCTA
RUN wget -O /tmp/gcta.zip https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-linux-kernel-3-x86_64.zip \
    && unzip /tmp/gcta.zip -d /usr/local/bin/ \
    && chmod +x /usr/local/bin/gcta*

# Download and install BCFtools
RUN wget -O /tmp/bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.23/bcftools-1.23.tar.bz2 \
    && tar -xjf /tmp/bcftools.tar.bz2 -C /tmp/ \
    && cd /tmp/bcftools-1.23 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd /app \
    && rm -rf /tmp/bcftools*

# Set the entrypoint to run the CLI as defined in pyproject.toml
ENTRYPOINT ["ideal-genom"]
