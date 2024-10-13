#!/bin/bash

# Update package lists and upgrade existing packages
echo "Updating package lists..."
sudo apt update && sudo apt upgrade -y

# Install essential tools
echo "Installing essential tools..."
sudo apt install -y build-essential wget curl

# Install Cutadapt
echo "Installing Cutadapt..."
pip install cutadapt

# Install FastX Toolkit
echo "Installing FastX Toolkit..."
sudo apt install -y fastx-toolkit

# Install BWA (Burrows-Wheeler Aligner)
echo "Installing BWA..."
sudo apt install -y bwa

# Install Samtools
echo "Installing Samtools..."
sudo apt install -y samtools

# Install CD-HIT
echo "Installing CD-HIT..."
sudo apt install -y cd-hit

# Install Python libraries if not already installed
echo "Installing necessary Python libraries..."
pip install biopython

# (Optional) Install FastQC for quality control
echo "Installing FastQC..."
sudo apt install -y fastqc

# Create directories for your project
echo "Creating project directories..."
mkdir -p Diploma-Thesis/split Diploma-Thesis/mapping Diploma-Thesis/clustering

# Print completion message
echo "Setup completed successfully! You can now run your script."