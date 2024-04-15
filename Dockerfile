# Use an official Python runtime as a parent image
FROM python:3.11-slim

# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY cge_comrare_pipeline /app

# Install any needed dependencies specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install wget to download the external application
RUN apt-get update && apt-get install -y wget

# Download the necessary app from the web
RUN wget -O /tmp/plink19.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip

# Extract the downloaded app
RUN unzip /tmp/plink19.zip -d /usr/local/bin/

# Make PLINK1.9 executable
# RUN chmod +x /usr/local/bin/plink/plink

# Make the entry point script executable
# RUN chmod +x entrypoint.sh

# Define the entry point command
CMD ["bash"]
