# Use an official Python runtime as a parent image
FROM python:3.11-slim

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
RUN git clone https://github.com/cge-tubingens/cge-comrare-pipeline .

# Install Python dependencies
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Copy the local code to the container
COPY . .

# Download and extract the necessary apps
RUN wget -O /tmp/plink19.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip \
    && unzip /tmp/plink19.zip -d /usr/local/bin/ \
    && wget -O /tmp/plink2.zip https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240105.zip \
    && unzip /tmp/plink2.zip -d /usr/local/bin/

# Set the entrypoint to run the Python CLI tool
ENTRYPOINT ["python", "-m", "cge_comrare_pipeline"]
